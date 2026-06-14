//! Analytic solvent-excluded surface (SES / Connolly / "rolling-probe" surface).
//!
//! Unlike the SAS (a creased union of spheres, see [`crate::surface`]), the SES is
//! smooth: it is what the *surface* of a probe sphere of radius `p` traces as the
//! probe rolls over the atoms. It is made of three patch types that join smoothly:
//!
//! * **Convex (contact) patches** — the part of each atom's vdW sphere the probe
//!   can touch. This is exactly the atom's SAS-exposed angular region, emitted at
//!   the vdW radius `R_i - p` (the SAS uses `R_i`). Reuses the icosphere clip.
//! * **Toroidal (saddle) patches** — where the probe rolls touching two atoms, its
//!   surface sweeps a torus; the SES is the inner tube arc between the two atoms'
//!   contact circles.
//! * **Concave (reentrant) patches** — where the probe rests against three atoms it
//!   sits still; the SES is a spherical triangle on that probe sphere.
//!
//! The contour (which neighbor arcs bound an atom, and where three atoms meet) is
//! derived self-contained as a circle arrangement on the sphere — each neighbor is
//! a circle `dir·e = cosθ`, and the exposed arcs of a circle are the angles not
//! inside any other neighbor's cap.

use nalgebra::{RealField, Vector3};

use crate::surface::{tessellate_atom, Cap, SurfaceMesh};

/// One neighbor of the current atom that actually cuts its sphere: the intersection
/// circle (`e`, `costheta`, `sintheta`) plus the neighbor's center and SAS radius.
pub(crate) struct NeighborCap<S>
where
    S: RealField + Copy,
{
    /// Unit direction from the current atom center to the neighbor.
    pub e: Vector3<S>,
    pub costheta: S,
    pub sintheta: S,
    /// Neighbor center (same — centered — frame as the current atom).
    pub center: Vector3<S>,
    /// Neighbor global atom index (for triple-point de-duplication / coloring).
    pub id: u32,
}

/// Orthonormal basis (u, v) spanning the plane perpendicular to unit vector `e`.
fn perp_basis<S>(e: &Vector3<S>) -> (Vector3<S>, Vector3<S>)
where
    S: RealField + Copy,
{
    let a = if e.x.abs() < S::from_f64(0.9).unwrap() {
        Vector3::new(S::one(), S::zero(), S::zero())
    } else {
        Vector3::new(S::zero(), S::one(), S::zero())
    };
    let u = (a - e * a.dot(e)).normalize();
    let v = e.cross(&u);
    (u, v)
}

fn clamp_unit<S: RealField + Copy>(x: S) -> S {
    if x > S::one() {
        S::one()
    } else if x < -S::one() {
        -S::one()
    } else {
        x
    }
}

/// Great-circle interpolation between two unit vectors.
fn slerp<S>(a: &Vector3<S>, b: &Vector3<S>, t: S) -> Vector3<S>
where
    S: RealField + Copy,
{
    let d = clamp_unit(a.dot(b));
    let omega = d.acos();
    let so = omega.sin();
    if so.abs() < S::from_f64(1e-7).unwrap() {
        return *a;
    }
    *a * (((S::one() - t) * omega).sin() / so) + *b * ((t * omega).sin() / so)
}

/// A point on intersection circle `k` (with basis u,v) at polar angle `phi`,
/// expressed as a unit direction from the current atom center.
fn circle_point<S>(cap: &NeighborCap<S>, u: &Vector3<S>, v: &Vector3<S>, phi: S) -> Vector3<S>
where
    S: RealField + Copy,
{
    cap.e * cap.costheta + (*u * phi.cos() + *v * phi.sin()) * cap.sintheta
}

const TWO_PI_F64: f64 = std::f64::consts::TAU;

/// Exposed arcs of circle `k`: the angular intervals `(phi_start, phi_end)` whose
/// points lie outside every other neighbor cap. A wrapping arc is returned with a
/// negative `phi_start` (it is not split at the seam). Used to sweep toroidal
/// patches; triple points (concave patches) are found separately by
/// [`circle_circle_dirs`], so endpoint provenance is not tracked here.
fn exposed_arcs<S>(
    k: usize,
    caps: &[NeighborCap<S>],
    u: &Vector3<S>,
    v: &Vector3<S>,
) -> Vec<(S, S)>
where
    S: RealField + Copy,
{
    let two_pi = S::from_f64(TWO_PI_F64).unwrap();
    // Boundary events around the circle: (+1, id) = enters cap, (-1, id) = leaves.
    let mut events: Vec<(S, i32, u32)> = Vec::new();
    let mut base_coverage: i32 = 0; // caps covering phi = 0
    let capk = &caps[k];

    for (m, capm) in caps.iter().enumerate() {
        if m == k {
            continue;
        }
        // circle_point(phi)·e_m = A + R·cos(phi - delta); inside cap m when > costheta_m.
        let big_a = capk.costheta * capk.e.dot(&capm.e);
        let big_b = capk.sintheta * u.dot(&capm.e);
        let big_c = capk.sintheta * v.dot(&capm.e);
        let r = (big_b * big_b + big_c * big_c).sqrt();
        let thr = capm.costheta - big_a;
        if r <= S::from_f64(1e-9).unwrap() {
            // Circle k is parallel to cap m's plane: uniformly inside or outside.
            if big_a > capm.costheta {
                base_coverage += 1; // entire circle buried by m
            }
            continue;
        }
        let ct = thr / r;
        if ct >= S::one() {
            continue; // circle k never enters cap m
        }
        if ct <= -S::one() {
            base_coverage += 1; // circle k entirely inside cap m
            continue;
        }
        let delta = big_c.atan2(big_b); // phase so that cos(phi-delta) peaks
        let alpha = clamp_unit(ct).acos();
        // Inside when phi-delta in (-alpha, alpha): interval (delta-alpha, delta+alpha).
        let mut s = delta - alpha;
        let mut e = delta + alpha;
        // Normalize to [0, 2π) and split wraps; count coverage of phi=0.
        let norm = |mut x: S| {
            while x < S::zero() {
                x += two_pi;
            }
            while x >= two_pi {
                x -= two_pi;
            }
            x
        };
        s = norm(s);
        e = norm(e);
        if s <= e {
            events.push((s, 1, capm.id));
            events.push((e, -1, capm.id));
        } else {
            // Inside-interval wraps past 0: covers [s,2π) ∪ [0,e]. The [0,e] part is
            // accounted by base_coverage (coverage at φ=0); the [s,2π) part enters at
            // s and never closes before the sweep ends.
            events.push((s, 1, capm.id));
            events.push((e, -1, capm.id));
            base_coverage += 1;
        }
    }

    // Sort by angle; at coincident angles process leaves (−1) before enters (+1) so
    // a point shared by two caps stays covered (closed-interval union), with a final
    // deterministic tie-break on id. Without this total order the sweep is
    // order-dependent and produces asymmetric results on symmetric inputs.
    events.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap()
            .then(a.1.cmp(&b.1))
            .then(a.2.cmp(&b.2))
    });

    let mut arcs: Vec<(S, S)> = Vec::new();
    if events.is_empty() {
        if base_coverage == 0 {
            arcs.push((S::zero(), two_pi)); // whole circle exposed
        }
        return arcs;
    }

    let mut coverage = base_coverage;
    let mut open_start: Option<S> = if coverage == 0 { Some(S::zero()) } else { None };
    for &(phi, delta, _id) in &events {
        if coverage == 0 {
            if let Some(start) = open_start.take() {
                if phi > start {
                    arcs.push((start, phi));
                }
            }
        }
        coverage += delta;
        if coverage == 0 {
            open_start = Some(phi);
        }
    }
    // Close a trailing exposed arc that runs to 2π, merging with a seam arc at 0.
    if coverage == 0 {
        if let Some(start) = open_start.take() {
            if !arcs.is_empty() && arcs[0].0 == S::zero() {
                arcs[0].0 = start - two_pi; // wrap: extend the leading arc backwards
            } else if two_pi > start {
                arcs.push((start, two_pi));
            }
        }
    }
    arcs
}

/// Append a quad (two triangles) given four vertex indices in order.
fn push_quad(mesh: &mut SurfaceMesh<impl RealField + Copy>, a: u32, b: u32, c: u32, d: u32) {
    mesh.indices.extend_from_slice(&[a, b, c, a, c, d]);
}

/// Running counts of emitted SES geometry (diagnostics).
#[derive(Default)]
pub struct SesStats {
    pub spindle_skips: usize,
    pub convex_tris: usize,
    pub toroidal_tris: usize,
    pub concave_tris: usize,
}

/// Build the SES of one atom's neighborhood into `mesh` (centered frame).
///
/// `caps` are the neighbors that cut atom `i`; `center_i`/`sas_r_i` are its center
/// and SAS radius; `vdw_i = sas_r_i - probe`. `ico_*` is a unit icosphere for the
/// convex patch. `stats` accumulates patch/triangle counts and singularity skips.
#[allow(clippy::too_many_arguments)]
pub(crate) fn atom_ses<S>(
    ico_v: &[Vector3<S>],
    ico_f: &[[usize; 3]],
    caps: &[NeighborCap<S>],
    center_i: Vector3<S>,
    sas_r_i: S,
    probe: S,
    atom_i: u32,
    subdiv: usize,
    stats: &mut SesStats,
    mesh: &mut SurfaceMesh<S>,
) where
    S: RealField + Copy,
{
    let vdw_i = sas_r_i - probe;

    // --- Convex (contact) patch: the SAS-exposed region at the vdW radius. ---
    let simple: Vec<Cap<S>> = caps
        .iter()
        .map(|c| Cap {
            axis: c.e,
            cos_theta: c.costheta,
        })
        .collect();
    let before_convex = mesh.indices.len();
    tessellate_atom(ico_v, ico_f, &simple, center_i, vdw_i, atom_i, mesh);
    stats.convex_tris += (mesh.indices.len() - before_convex) / 3;

    // --- Toroidal + concave patches along each neighbor's exposed arcs. ---
    let tube_steps = 2 + 2 * subdiv as i32;
    let seg_ang = S::from_f64(0.5).unwrap() / S::from_f64((1 + subdiv) as f64).unwrap();

    // Patches meet exactly at shared curves (no overlap): the convex↔toroidal and
    // toroidal↔concave seams are left as boundary loops and closed by a hole-filling
    // pass on the welded mesh, which is robust to whatever produced the gap.
    for k in 0..caps.len() {
        let capk = &caps[k];
        // Each atom pair shares one toroidal patch; emit it once, from the
        // lower-index atom, to avoid a doubled (z-fighting) tube.
        if atom_i > capk.id {
            continue;
        }
        // Spindle singularity: tube self-intersects when the contact-circle radius
        // (distance from the i–k axis to the circle) is below the probe radius.
        let rho = sas_r_i * capk.sintheta;
        if rho < probe {
            stats.spindle_skips += 1;
            continue;
        }
        let (u, v) = perp_basis(&capk.e);
        let arcs = exposed_arcs(k, caps, &u, &v);

        for (phi_a, phi_b) in arcs {
            let span = phi_b - phi_a;
            if span <= S::from_f64(1e-5).unwrap() {
                continue;
            }
            // ceil(span / seg_ang), without a generic Scalar -> primitive cast.
            let mut n_phi = 1i32;
            while seg_ang * S::from_f64(n_phi as f64).unwrap() < span && n_phi < 512 {
                n_phi += 1;
            }

            // Sweep the tube as a (phi, psi) grid; emit a quad strip over the exact
            // exposed arc (no overlap — seams are sealed later by hole-filling).
            let before_tor = mesh.indices.len();
            let mut prev_row: Option<Vec<u32>> = None;
            for ip in 0..=n_phi {
                let t = S::from_f64(ip as f64).unwrap() / S::from_f64(n_phi as f64).unwrap();
                let phi = phi_a + span * t;
                let dir = circle_point(capk, &u, &v, phi); // unit, from center_i
                let pc = center_i + dir * sas_r_i; // probe center (on both SAS spheres)
                // Tube end directions (probe-center -> contact point), both unit.
                let n_i = -dir; // toward contact on atom i
                let n_k = (capk.center - pc).normalize(); // toward contact on atom k
                let mut row: Vec<u32> = Vec::with_capacity(tube_steps as usize + 1);
                for it in 0..=tube_steps {
                    let psi = S::from_f64(it as f64).unwrap()
                        / S::from_f64(tube_steps as f64).unwrap();
                    let nrm = slerp(&n_i, &n_k, psi); // unit, probe-center -> surface
                    let pos = pc + nrm * probe;
                    let idx = mesh.vertices.len() as u32;
                    mesh.vertices.push(pos);
                    mesh.normals.push(-nrm); // concave: face toward the probe center
                    // Color by whichever atom this side of the tube belongs to.
                    mesh.atom_ids.push(if psi < S::from_f64(0.5).unwrap() {
                        atom_i
                    } else {
                        capk.id
                    });
                    row.push(idx);
                }
                if let Some(prev) = &prev_row {
                    for it in 0..tube_steps as usize {
                        push_quad(mesh, prev[it], prev[it + 1], row[it + 1], row[it]);
                    }
                }
                prev_row = Some(row);
            }
            stats.toroidal_tris += (mesh.indices.len() - before_tor) / 3;
        }
    }

    // --- Concave (reentrant) patches at triple points. ---
    // Enumerate triple points directly as circle–circle intersections (robust at
    // quadruple junctions, where arc-endpoint detection drops zero-length arcs):
    // for each neighbor pair (ka, kb), the two circles meet at up to 2 directions;
    // a meeting that lies outside every other cap is an exposed triple (i, ka, kb).
    // Emitted once per triple by the lowest-index atom (which always sees the pair).
    let before_conc = mesh.indices.len();
    for a in 0..caps.len() {
        for b in (a + 1)..caps.len() {
            if atom_i > caps[a].id || atom_i > caps[b].id {
                continue; // de-dup: only the lowest-index atom of the triple emits
            }
            for dir in circle_circle_dirs(&caps[a], &caps[b]) {
                // Exposed iff outside every other neighbor cap.
                let exposed = caps.iter().enumerate().all(|(j, c)| {
                    j == a || j == b || dir.dot(&c.e) <= c.costheta + S::from_f64(1e-6).unwrap()
                });
                if !exposed {
                    continue;
                }
                let pc = center_i + dir * sas_r_i;
                let c_i = -dir; // probe-center -> contact on atom i
                let c_a = (caps[a].center - pc).normalize();
                let c_b = (caps[b].center - pc).normalize();
                // Skip degenerate triangles (near-tangent circles → coincident
                // contacts → a zero-area patch whose normal is garbage and renders
                // as a stray colored spot).
                let area2 = (c_a - c_i).cross(&(c_b - c_i)).norm();
                if area2 < S::from_f64(1e-5).unwrap() {
                    continue;
                }
                emit_concave_triangle(pc, probe, c_i, c_a, c_b, atom_i, subdiv, mesh);
            }
        }
    }
    stats.concave_tris += (mesh.indices.len() - before_conc) / 3;
}

/// The (up to two) unit directions where intersection circles `a` and `b` meet on
/// the sphere — solving `dir·e_a = cosθ_a`, `dir·e_b = cosθ_b`, `|dir| = 1`.
fn circle_circle_dirs<S>(a: &NeighborCap<S>, b: &NeighborCap<S>) -> Vec<Vector3<S>>
where
    S: RealField + Copy,
{
    let d = a.e.dot(&b.e);
    let denom = S::one() - d * d;
    if denom <= S::from_f64(1e-9).unwrap() {
        return Vec::new(); // axes (near-)parallel: circles don't cross cleanly
    }
    // dir = ca·e_a + cb·e_b + gamma·(e_a×e_b normalized); n·e_a = n·e_b = 0.
    let ca = (a.costheta - b.costheta * d) / denom;
    let cb = (b.costheta - a.costheta * d) / denom;
    let g2 = S::one() - (ca * ca + cb * cb + (S::one() + S::one()) * ca * cb * d);
    if g2 < S::zero() {
        return Vec::new(); // circles do not intersect
    }
    let n = a.e.cross(&b.e).normalize();
    let base = a.e * ca + b.e * cb;
    let g = g2.sqrt();
    if g <= S::from_f64(1e-7).unwrap() {
        return vec![base]; // tangent: single point
    }
    vec![base + n * g, base - n * g]
}

/// A concave spherical triangle on the probe sphere (center `pc`, radius `probe`)
/// through the three contact directions, recursively subdivided `subdiv` times.
/// Normals point toward the probe center (the surface is reentrant/concave).
fn emit_concave_triangle<S>(
    pc: Vector3<S>,
    probe: S,
    a: Vector3<S>,
    b: Vector3<S>,
    c: Vector3<S>,
    atom_i: u32,
    subdiv: usize,
    mesh: &mut SurfaceMesh<S>,
) where
    S: RealField + Copy,
{
    // Recursive 4-way subdivision on the sphere of unit directions.
    fn rec<S: RealField + Copy>(
        a: Vector3<S>,
        b: Vector3<S>,
        c: Vector3<S>,
        depth: usize,
        pc: Vector3<S>,
        probe: S,
        atom_i: u32,
        mesh: &mut SurfaceMesh<S>,
    ) {
        if depth == 0 {
            let base = mesh.vertices.len() as u32;
            for d in [a, b, c] {
                mesh.vertices.push(pc + d * probe);
                mesh.normals.push(-d); // toward probe center
                mesh.atom_ids.push(atom_i);
            }
            mesh.indices.extend_from_slice(&[base, base + 1, base + 2]);
            return;
        }
        let ab = (a + b).normalize();
        let bc = (b + c).normalize();
        let ca = (c + a).normalize();
        rec(a, ab, ca, depth - 1, pc, probe, atom_i, mesh);
        rec(ab, b, bc, depth - 1, pc, probe, atom_i, mesh);
        rec(ca, bc, c, depth - 1, pc, probe, atom_i, mesh);
        rec(ab, bc, ca, depth - 1, pc, probe, atom_i, mesh);
    }
    rec(a, b, c, subdiv.min(3), pc, probe, atom_i, mesh);
}
