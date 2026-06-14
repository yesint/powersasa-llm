//! Tessellated-surface tests for `PowerSasa::surface_mesh`.

use nalgebra::Vector3;
use powersasa::PowerSasa;

/// Surface area of an indexed triangle mesh.
fn mesh_area(verts: &[Vector3<f64>], indices: &[u32]) -> f64 {
    let mut area = 0.0;
    for tri in indices.chunks_exact(3) {
        let a = verts[tri[0] as usize];
        let b = verts[tri[1] as usize];
        let c = verts[tri[2] as usize];
        area += 0.5 * (b - a).cross(&(c - a)).norm();
    }
    area
}

#[test]
fn lone_sphere_is_closed_and_correct_area() {
    let coords = vec![Vector3::new(0.0_f64, 0.0, 0.0)];
    let radii = vec![1.5_f64];
    let mut ps = PowerSasa::<f64>::new(
        coords.into_iter(),
        radii.into_iter(),
        true,
        false,
        false,
        false,
    );
    ps.calc_all().unwrap();
    let m = ps.surface_mesh(3);

    assert!(!m.vertices.is_empty(), "lone sphere produced no geometry");
    assert_eq!(m.vertices.len(), m.normals.len());
    assert_eq!(m.vertices.len(), m.atom_ids.len());
    assert_eq!(m.indices.len() % 3, 0);

    // Every vertex sits on the sphere of radius R, with an outward unit normal.
    let r = 1.5;
    for (v, n) in m.vertices.iter().zip(&m.normals) {
        assert!((v.norm() - r).abs() < 1e-9, "vertex off the sphere: {}", v.norm());
        assert!((n.norm() - 1.0).abs() < 1e-9, "normal not unit");
    }

    // Triangle area approximates the analytic sphere area (4πR²) from below.
    let exact = 4.0 * std::f64::consts::PI * r * r;
    let got = mesh_area(&m.vertices, &m.indices);
    assert!(
        got > 0.95 * exact && got <= exact * 1.001,
        "sphere mesh area {got:.4} vs analytic {exact:.4}"
    );
}

#[test]
fn two_overlapping_spheres_have_no_buried_vertices() {
    let r = 1.0_f64;
    let d = 1.5_f64; // < 2r → the spheres overlap
    let c0 = Vector3::new(0.0, 0.0, 0.0);
    let c1 = Vector3::new(d, 0.0, 0.0);
    let coords = vec![c0, c1];
    let radii = vec![r, r];
    let mut ps = PowerSasa::<f64>::new(
        coords.into_iter(),
        radii.into_iter(),
        true,
        false,
        false,
        false,
    );
    ps.calc_all().unwrap();
    let m = ps.surface_mesh(3);
    assert!(!m.vertices.is_empty());

    // Both atoms contribute exposed patches.
    assert!(m.atom_ids.iter().any(|&a| a == 0), "atom 0 missing");
    assert!(m.atom_ids.iter().any(|&a| a == 1), "atom 1 missing");

    // No kept vertex lies strictly inside the *other* ball (the buried region
    // was clipped away; boundary verts are snapped exactly onto the cut circle).
    let eps = 1e-3;
    for (v, &a) in m.vertices.iter().zip(&m.atom_ids) {
        let other = if a == 0 { c1 } else { c0 };
        assert!(
            (v - other).norm() >= r - eps,
            "vertex of atom {a} is buried inside the other sphere: dist {}",
            (v - other).norm()
        );
    }
}

#[test]
fn ses_mesh_is_finite_nonempty_and_in_range() {
    // A small tetrahedral cluster (nm), SAS radii = vdW(0.16) + probe(0.14).
    let coords = vec![
        Vector3::new(0.00_f64, 0.0, 0.0),
        Vector3::new(0.24, 0.0, 0.0),
        Vector3::new(0.12, 0.208, 0.0),
        Vector3::new(0.12, 0.069, 0.196),
    ];
    let radii = vec![0.30_f64; 4];
    let mut ps = PowerSasa::<f64>::new(
        coords.into_iter(),
        radii.into_iter(),
        true,
        false,
        false,
        false,
    );
    ps.calc_all().unwrap();
    let m = ps.ses_mesh(0.14, 2);

    assert!(!m.vertices.is_empty(), "SES mesh is empty");
    assert_eq!(m.vertices.len(), m.normals.len());
    assert_eq!(m.vertices.len(), m.atom_ids.len());
    assert_eq!(m.indices.len() % 3, 0);
    for v in &m.vertices {
        assert!(
            v.x.is_finite() && v.y.is_finite() && v.z.is_finite(),
            "non-finite SES vertex {v:?}"
        );
    }
    for nrm in &m.normals {
        assert!((nrm.norm() - 1.0).abs() < 1e-3, "non-unit SES normal");
    }
    assert!(*m.indices.iter().max().unwrap() < m.vertices.len() as u32);
    assert!(m.atom_ids.iter().all(|&a| (a as usize) < 4));
}

#[test]
fn fully_buried_atom_emits_nothing() {
    // A tiny atom sitting well inside a big one contributes no surface.
    let coords = vec![Vector3::new(0.0_f64, 0.0, 0.0), Vector3::new(0.5, 0.0, 0.0)];
    let radii = vec![2.0_f64, 0.1_f64];
    let mut ps = PowerSasa::<f64>::new(
        coords.into_iter(),
        radii.into_iter(),
        true,
        false,
        false,
        false,
    );
    ps.calc_all().unwrap();
    let m = ps.surface_mesh(2);
    assert!(!m.vertices.is_empty(), "big atom should still have a surface");
    assert!(
        m.atom_ids.iter().all(|&a| a == 0),
        "the engulfed atom 1 must emit no geometry"
    );
}
