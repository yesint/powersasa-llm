use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use nalgebra::Vector3;
use powersasa::PowerSasa;

#[derive(Debug, Clone)]
struct GoldenCase {
    file_path: String,
    expected_sasa: f32,
    expected_vol: f32,
}

fn load_case_file(path: &str) -> Result<(Vec<Vector3<f32>>, Vec<f32>), String> {
    let file = fs::File::open(path).map_err(|e| format!("Could not open {path}: {e}"))?;
    let reader = BufReader::new(file);

    let mut coords = Vec::new();
    let mut weights = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading {path}: {e}"))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let Some(x) = parts.next() else { continue };
        let Some(y) = parts.next() else { continue };
        let Some(z) = parts.next() else { continue };
        let Some(radius) = parts.next() else {
            continue;
        };

        let x: f32 = x.parse().map_err(|e| format!("Bad x in {path}: {e}"))?;
        let y: f32 = y.parse().map_err(|e| format!("Bad y in {path}: {e}"))?;
        let z: f32 = z.parse().map_err(|e| format!("Bad z in {path}: {e}"))?;
        let radius: f32 = radius
            .parse()
            .map_err(|e| format!("Bad radius in {path}: {e}"))?;

        coords.push(Vector3::new(x, y, z));
        weights.push(radius + 0.14);
    }

    if coords.is_empty() {
        return Err(format!("No atoms loaded from {path}"));
    }

    Ok((coords, weights))
}

fn load_golden_manifest(manifest_path: &Path) -> Result<Vec<GoldenCase>, String> {
    let file = fs::File::open(manifest_path)
        .map_err(|e| format!("Could not open manifest {}: {e}", manifest_path.display()))?;
    let reader = BufReader::new(file);

    let manifest_dir = manifest_path
        .parent()
        .ok_or_else(|| format!("Manifest has no parent: {}", manifest_path.display()))?;

    let mut cases = Vec::new();

    for line in reader.lines() {
        let line =
            line.map_err(|e| format!("Error reading manifest {}: {e}", manifest_path.display()))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let Some(file_path) = parts.next() else {
            continue;
        };
        let Some(expected_sasa) = parts.next() else {
            continue;
        };
        let Some(expected_vol) = parts.next() else {
            continue;
        };

        let expected_sasa: f32 = expected_sasa
            .parse()
            .map_err(|e| format!("Bad expected_sasa '{expected_sasa}': {e}"))?;
        let expected_vol: f32 = expected_vol
            .parse()
            .map_err(|e| format!("Bad expected_vol '{expected_vol}': {e}"))?;

        let listed = PathBuf::from(file_path);
        let resolved = if listed.is_absolute() {
            listed
        } else {
            manifest_dir.join(&listed)
        };

        cases.push(GoldenCase {
            file_path: resolved.to_string_lossy().into_owned(),
            expected_sasa,
            expected_vol,
        });
    }

    if cases.is_empty() {
        return Err(format!("No cases in manifest {}", manifest_path.display()));
    }

    Ok(cases)
}

fn find_manifest() -> PathBuf {
    PathBuf::from("tests/testdata/sasa_cases/golden_912340c.txt")
}

#[test]
fn sasatest_matches_golden_912340c() -> Result<(), String> {
    const EPSILON: f32 = 1e-3;
    let manifest = find_manifest();

    let cases = load_golden_manifest(&manifest)?;

    for c in &cases {
        let (coords, weights) = load_case_file(&c.file_path)?;

        let mut ps = PowerSasa::<f32>::new(
            coords.clone().into_iter(),
            weights.clone().into_iter(),
            true,
            false,
            true,
            false,
        );
        
        ps.calc_all()
            .map_err(|e| format!("calc_sasa_all failed for {}: {e}", c.file_path))?;

        let total_sasa: f32 = ps.per_atom_sasa().iter().copied().sum();
        let total_vol: f32 = ps.per_atom_vol().iter().copied().sum();

        if (total_sasa - c.expected_sasa).abs() > EPSILON {
            return Err(format!(
                "SASA mismatch for {}: got {:.6}, expected {:.6}",
                c.file_path, total_sasa, c.expected_sasa
            ));
        }
        if (total_vol - c.expected_vol).abs() > EPSILON {
            return Err(format!(
                "Volume mismatch for {}: got {:.6}, expected {:.6}",
                c.file_path, total_vol, c.expected_vol
            ));
        }
    }

    Ok(())
}
