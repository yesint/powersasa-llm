//! Triangle-mesh extraction of the solvent-accessible surface (SAS).
//!
//! The SAS is the boundary of the union of balls of radius `R_i = vdW_i + probe`.
//! A point `c_i + R_i·dir` on atom `i`'s sphere is *exposed* iff it lies outside
//! every other ball, i.e. `dir·e_k ≤ cosθ_k` for all neighbors `k`, where `e_k` is
//! the unit direction to neighbor `k` and
//! `cosθ_k = (d² + R_i² − R_k²) / (2·R_i·d)` (`d` = center distance). This is
//! exactly the per-atom burial test PowerSASA already uses for its area integral
//! (see `PowerSasa::calc_single`), so the mesh is the same exact SAS the analytic
//! area is computed over.
//!
//! We tessellate each exposed atom by clipping a unit **icosphere** against its
//! neighbor caps with spherical Sutherland–Hodgman: every kept vertex lies on the
//! SAS sphere and every cut edge is snapped onto the analytic intersection circle.
//! Robust by construction — holes, multiple caps, fully-exposed and fully-buried
//! atoms all fall out of the same path — and embarrassingly parallel per atom.

use std::collections::HashMap;

use nalgebra::{RealField, Vector3};

/// An indexed triangle mesh of a molecular surface, in the **original** input
/// coordinate frame (the power diagram's center shift is added back).
pub struct SurfaceMesh<Scalar>
where
    Scalar: RealField + Copy,
{
    /// World-space vertex positions.
    pub vertices: Vec<Vector3<Scalar>>,
    /// Outward unit normals (the sphere direction at each vertex).
    pub normals: Vec<Vector3<Scalar>>,
    /// Triangle indices into `vertices`.
    pub indices: Vec<u32>,
    /// Per-vertex source atom index (the generator the patch belongs to).
    pub atom_ids: Vec<u32>,
}

impl<Scalar> Default for SurfaceMesh<Scalar>
where
    Scalar: RealField + Copy,
{
    fn default() -> Self {
        Self {
            vertices: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
            atom_ids: Vec::new(),
        }
    }
}

/// A neighbor cap on the unit sphere: directions with `dir·axis > cos_theta` are
/// buried by that neighbor and must be clipped away.
pub(crate) struct Cap<Scalar>
where
    Scalar: RealField + Copy,
{
    pub axis: Vector3<Scalar>,
    pub cos_theta: Scalar,
}

/// Build a unit icosphere by recursively subdividing an icosahedron `subdiv`
/// times (0 → 20 triangles, 1 → 80, 2 → 320, 3 → 1280, …). Returns unit-length
/// vertex directions and triangle index triples. Winding is not made consistent —
/// the caller derives outward normals from the directions and renders two-sided.
pub(crate) fn unit_icosphere<Scalar>(subdiv: usize) -> (Vec<Vector3<Scalar>>, Vec<[usize; 3]>)
where
    Scalar: RealField + Copy,
{
    let t = (1.0 + 5.0_f64.sqrt()) * 0.5;
    let raw: [[f64; 3]; 12] = [
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut verts: Vec<Vector3<Scalar>> = raw
        .iter()
        .map(|v| {
            let x = Vector3::new(
                Scalar::from_f64(v[0]).unwrap(),
                Scalar::from_f64(v[1]).unwrap(),
                Scalar::from_f64(v[2]).unwrap(),
            );
            x.normalize()
        })
        .collect();

    let mut faces: Vec<[usize; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];

    for _ in 0..subdiv {
        // Midpoint cache so shared edges produce a single shared vertex.
        let mut cache: HashMap<(usize, usize), usize> = HashMap::new();
        let mut new_faces: Vec<[usize; 3]> = Vec::with_capacity(faces.len() * 4);
        let mut midpoint = |a: usize, b: usize, verts: &mut Vec<Vector3<Scalar>>| -> usize {
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(&m) = cache.get(&key) {
                return m;
            }
            let m = (verts[a] + verts[b]).normalize();
            let idx = verts.len();
            verts.push(m);
            cache.insert(key, idx);
            idx
        };
        for f in &faces {
            let a = midpoint(f[0], f[1], &mut verts);
            let b = midpoint(f[1], f[2], &mut verts);
            let c = midpoint(f[2], f[0], &mut verts);
            new_faces.push([f[0], a, c]);
            new_faces.push([f[1], b, a]);
            new_faces.push([f[2], c, b]);
            new_faces.push([a, b, c]);
        }
        faces = new_faces;
    }

    (verts, faces)
}

/// Clip a spherical polygon (list of unit directions) against one cap, keeping the
/// part with `dir·axis ≤ cos_theta`. Standard Sutherland–Hodgman; crossings are
/// found by bisection along the great-circle arc so they land exactly on the cut
/// circle (`dir·axis == cos_theta`). Tiny icosphere triangles make the
/// at-most-one-crossing-per-edge assumption safe (a vanishing low-`subdiv`
/// approximation when a cut circle is smaller than a triangle).
fn clip_against_cap<Scalar>(poly: &[Vector3<Scalar>], cap: &Cap<Scalar>) -> Vec<Vector3<Scalar>>
where
    Scalar: RealField + Copy,
{
    let n = poly.len();
    if n == 0 {
        return Vec::new();
    }
    let mut out: Vec<Vector3<Scalar>> = Vec::with_capacity(n + 2);
    for k in 0..n {
        let cur = poly[k];
        let nxt = poly[(k + 1) % n];
        let cur_in = cur.dot(&cap.axis) <= cap.cos_theta;
        let nxt_in = nxt.dot(&cap.axis) <= cap.cos_theta;
        if nxt_in {
            if !cur_in {
                out.push(arc_cross(&cur, &nxt, cap));
            }
            out.push(nxt);
        } else if cur_in {
            out.push(arc_cross(&cur, &nxt, cap));
        }
    }
    out
}

/// Bisection for the point on the arc `a → b` where `dir·axis == cos_theta`. The
/// arc is traced by the normalized linear interpolant (which lies on the great
/// circle through a, b). `a` and `b` must straddle the cut; the bracket is kept
/// regardless of which endpoint is the inside one, so either orientation works.
fn arc_cross<Scalar>(a: &Vector3<Scalar>, b: &Vector3<Scalar>, cap: &Cap<Scalar>) -> Vector3<Scalar>
where
    Scalar: RealField + Copy,
{
    let half = Scalar::from_f64(0.5).unwrap();
    let a_inside = a.dot(&cap.axis) - cap.cos_theta <= Scalar::zero();
    let mut lo = Scalar::zero();
    let mut hi = Scalar::one();
    let mut mid = half;
    for _ in 0..24 {
        mid = (lo + hi) * half;
        let dir = (a * (Scalar::one() - mid) + b * mid).normalize();
        let inside = dir.dot(&cap.axis) - cap.cos_theta <= Scalar::zero();
        // Move the endpoint that is on the same side as `a` toward the root.
        if inside == a_inside {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (a * (Scalar::one() - mid) + b * mid).normalize()
}

/// Tessellate one atom's exposed SAS patch into `mesh`. The unit icosphere
/// (`ico_v`/`ico_f`) is clipped triangle-by-triangle against `caps`; survivors are
/// fan-triangulated. `world_center`/`rad` place the patch in world space.
pub(crate) fn tessellate_atom<Scalar>(
    ico_v: &[Vector3<Scalar>],
    ico_f: &[[usize; 3]],
    caps: &[Cap<Scalar>],
    world_center: Vector3<Scalar>,
    rad: Scalar,
    atom_id: u32,
    mesh: &mut SurfaceMesh<Scalar>,
) where
    Scalar: RealField + Copy,
{
    for f in ico_f {
        let mut poly: Vec<Vector3<Scalar>> = vec![ico_v[f[0]], ico_v[f[1]], ico_v[f[2]]];
        for cap in caps {
            if poly.len() < 3 {
                break;
            }
            poly = clip_against_cap(&poly, cap);
        }
        if poly.len() < 3 {
            continue;
        }
        // Fan-triangulate the surviving (convex) spherical polygon.
        let base = mesh.vertices.len() as u32;
        for dir in &poly {
            mesh.vertices.push(world_center + dir * rad);
            mesh.normals.push(*dir);
            mesh.atom_ids.push(atom_id);
        }
        for j in 1..(poly.len() as u32 - 1) {
            mesh.indices.push(base);
            mesh.indices.push(base + j);
            mesh.indices.push(base + j + 1);
        }
    }
}
