use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use nalgebra::Vector3;
use powersasa_rust::{PowerSasa, PowerSasaSettings};

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
        let Some(radius) = parts.next() else { continue };

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
    let repo_root = manifest_dir
        .parent()
        .and_then(|p| p.parent())
        .ok_or_else(|| format!("Could not derive repo root from {}", manifest_path.display()))?;

    let mut cases = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading manifest {}: {e}", manifest_path.display()))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let Some(file_path) = parts.next() else { continue };
        let Some(expected_sasa) = parts.next() else { continue };
        let Some(expected_vol) = parts.next() else { continue };

        let expected_sasa: f32 = expected_sasa
            .parse()
            .map_err(|e| format!("Bad expected_sasa '{expected_sasa}': {e}"))?;
        let expected_vol: f32 = expected_vol
            .parse()
            .map_err(|e| format!("Bad expected_vol '{expected_vol}': {e}"))?;

        let listed = PathBuf::from(file_path);
        let resolved = if listed.is_absolute() {
            listed
        } else if listed.exists() {
            listed
        } else if (repo_root.join(&listed)).exists() {
            repo_root.join(&listed)
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
    let rel = PathBuf::from("testdata/sasa_cases/golden_912340c.txt");
    for candidate in [rel.clone(), PathBuf::from("../").join(&rel), PathBuf::from("../../").join(&rel)] {
        if candidate.exists() {
            return candidate;
        }
    }
    rel
}

fn main() -> Result<(), String> {
    const EPSILON: f32 = 1e-3;
    let manifest = find_manifest();
    println!("Precision: float");
    println!("Manifest: {}", manifest.display());

    let cases = load_golden_manifest(&manifest)?;
    let mut all_passed = true;

    for c in &cases {
        let (coords, weights) = load_case_file(&c.file_path)?;

        let mut ps = PowerSasa::<f32>::new(
            coords.clone().into_iter(),
            weights.clone().into_iter(),
            PowerSasaSettings {
                with_sasa: true,
                with_dsasa: false,
                with_vol: true,
                with_dvol: false,
            },
        );
        ps.calc_sasa_all().map_err(|e| format!("calc_sasa_all failed: {e}"))?;

        let total_sasa: f32 = ps.get_sasa().iter().copied().sum();
        let total_vol: f32 = ps.get_vol().iter().copied().sum();

        let sasa_ok = (total_sasa - c.expected_sasa).abs() <= EPSILON;
        let vol_ok = (total_vol - c.expected_vol).abs() <= EPSILON;
        all_passed = all_passed && sasa_ok && vol_ok;

        println!("Case: {}", c.file_path);
        println!(
            "  SASA:   {:>12.6}  expected {:>12.6}  [{}]",
            total_sasa,
            c.expected_sasa,
            if sasa_ok { "OK" } else { "FAIL" }
        );
        println!(
            "  Volume: {:>12.6}  expected {:>12.6}  [{}]",
            total_vol,
            c.expected_vol,
            if vol_ok { "OK" } else { "FAIL" }
        );
    }

    if !all_passed {
        return Err("REGRESSION FAILURE: one or more cases differ from 912340c golden values.".to_string());
    }

    println!("REGRESSION PASS: all cases match 912340c golden values.");
    Ok(())
}
