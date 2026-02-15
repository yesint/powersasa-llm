use std::fmt;

use nalgebra::Vector3;

use crate::power_diagram::PowerDiagram;

#[derive(Debug, Clone)]
pub struct PowerSasaError;

impl fmt::Display for PowerSasaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PowerSasaError")
    }
}

impl std::error::Error for PowerSasaError {}

pub struct PowerSasa<Scalar>
where
    Scalar: nalgebra::RealField + Copy,
{
    power_diagram: PowerDiagram<Scalar>,
    with_sasa: bool,
    with_dsasa: bool,
    with_vol: bool,
    with_dvol: bool,
    sasa: Vec<Scalar>,
    vol: Vec<Scalar>,
    dsasa: Vec<Vector3<Scalar>>,
    dvol: Vec<Vector3<Scalar>>,
}

impl<Scalar> PowerSasa<Scalar>
where
    Scalar: nalgebra::RealField + Copy,
{
    pub fn new(
        coords: impl Iterator<Item = Vector3<Scalar>> + Clone,
        radii: impl Iterator<Item = Scalar> + Clone,
        with_sasa: bool,
        with_dsasa: bool,
        with_vol: bool,
        with_dvol: bool,
    ) -> Self {
        let coords_vec: Vec<Vector3<Scalar>> = coords.collect();
        let radii_vec: Vec<Scalar> = radii.collect();
        let mut bond_to = Vec::with_capacity(coords_vec.len());
        bond_to.push(0);
        for i in 1..coords_vec.len() {
            bond_to.push((i - 1) as i32);
        }

        let power_diagram = PowerDiagram::create(
            coords_vec.len(),
            coords_vec.clone().into_iter(),
            radii_vec.clone().into_iter(),
            bond_to.into_iter(),
        );

        let n = power_diagram.get_points().len();

        Self {
            power_diagram,
            with_sasa,
            with_dsasa,
            with_vol,
            with_dvol,
            sasa: vec![Scalar::zero(); n],
            vol: vec![Scalar::zero(); n],
            dsasa: vec![Vector3::zeros(); n],
            dvol: vec![Vector3::zeros(); n],
        }
    }

    pub fn calc_sasa_single(&mut self, iatom: usize) -> Result<(), PowerSasaError> {
        if iatom >= self.power_diagram.get_points().len() {
            return Err(PowerSasaError);
        }
        Ok(())
    }

    pub fn calc_sasa_all(&mut self) -> Result<(), PowerSasaError> {
        for i in 0..self.power_diagram.get_points().len() {
            self.calc_sasa_single(i)?;
        }
        Ok(())
    }

    pub fn update_coords(
        &mut self,
        coords: impl Iterator<Item = Vector3<Scalar>>,
        radii: impl Iterator<Item = Scalar>,
        size: usize,
    ) {
        self.power_diagram.recalculate(coords, radii, size);
        let n = self.power_diagram.get_points().len();
        self.sasa.resize(n, Scalar::zero());
        self.vol.resize(n, Scalar::zero());
        self.dsasa.resize(n, Vector3::zeros());
        self.dvol.resize(n, Vector3::zeros());
    }

    pub fn get_sasa(&self) -> &[Scalar] {
        &self.sasa
    }

    pub fn get_vol(&self) -> &[Scalar] {
        &self.vol
    }

    pub fn get_dsasa(&self) -> &[Vector3<Scalar>] {
        &self.dsasa
    }

    pub fn get_dvol(&self) -> &[Vector3<Scalar>] {
        &self.dvol
    }

    pub fn with_sasa(&self) -> bool {
        self.with_sasa
    }

    pub fn with_dsasa(&self) -> bool {
        self.with_dsasa
    }

    pub fn with_vol(&self) -> bool {
        self.with_vol
    }

    pub fn with_dvol(&self) -> bool {
        self.with_dvol
    }
}
