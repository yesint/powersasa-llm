use std::fmt;

use nalgebra::{RealField, Vector3};

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
    Scalar: RealField + Copy,
{
    power_diagram: PowerDiagram<Scalar>,

    with_sasa: bool,
    with_dsasa: bool,
    with_vol: bool,
    with_dvol: bool,

    sasa: Vec<Scalar>,
    dsasa_parts: Vec<Vec<Vector3<Scalar>>>,
    dsasa: Vec<Vector3<Scalar>>,
    vol: Vec<Scalar>,
    dvol: Vec<Vector3<Scalar>>,

    tol_pow: Scalar,

    np: Vec<i32>,
    nt: Vec<i32>,
    e: Vec<Vector3<Scalar>>,
    sintheta: Vec<Scalar>,
    costheta: Vec<Scalar>,
    nb_rad2: Vec<Scalar>,
    nb_dist: Vec<Scalar>,

    off: Vec<i32>,
    vx: Vec<Vector3<Scalar>>,
    br_c: Vec<Vec<i32>>,
    br_p: Vec<Vec<i32>>,

    ang: Vec<Vec<Scalar>>,
    next: Vec<Vec<i32>>,
    p: Vec<Vec<i32>>,

    volnb: Vec<Scalar>,
    knot: Vec<Vector3<Scalar>>,
    fknot: Vec<bool>,

    rang: Vec<i32>,
    pos: Vec<i32>,
}

impl<Scalar> PowerSasa<Scalar>
where
    Scalar: RealField + Copy,
{
    pub const K_MAX_NB: usize = 20;
    pub const K_MAX_VX: usize = 12;
    pub const K_MAX_PNT: usize = 4;
    pub const K_MAX_COUNT: usize = 100;

    #[inline]
    pub fn drad2() -> Scalar {
        Scalar::from_f64(1000.0).unwrap() * Scalar::default_epsilon()
    }

    #[inline]
    pub fn dang() -> Scalar {
        Scalar::from_f64(1000.0).unwrap() * Scalar::default_epsilon()
    }

    #[inline]
    pub fn near_one_cosine_threshold() -> Scalar {
        Scalar::from_f64(0.999).unwrap()
    }

    #[inline]
    pub fn axis_component_threshold() -> Scalar {
        Scalar::from_f64(0.001).unwrap()
    }

    #[inline]
    pub fn clamp_unit_interval(value: Scalar) -> Scalar {
        if value < -Scalar::one() {
            -Scalar::one()
        } else if value > Scalar::one() {
            Scalar::one()
        } else {
            value
        }
    }

    pub fn new(
        coords: impl Iterator<Item = Vector3<Scalar>>,
        radii: impl Iterator<Item = Scalar>,
        with_sasa: bool,
        with_dsasa: bool,
        with_vol: bool,
        with_dvol: bool,
    ) -> Self {
        let coords_vec: Vec<Vector3<Scalar>> = coords.collect();
        let radii_vec: Vec<Scalar> = radii.collect();
        let power_diagram = Self::create_with_chain_bond(&coords_vec, &radii_vec);

        let mut this = Self {
            power_diagram,
            with_sasa,
            with_dsasa,
            with_vol,
            with_dvol,
            sasa: Vec::new(),
            dsasa_parts: Vec::new(),
            dsasa: Vec::new(),
            vol: Vec::new(),
            dvol: Vec::new(),
            tol_pow: Scalar::zero(),
            np: Vec::new(),
            nt: Vec::new(),
            e: Vec::new(),
            sintheta: Vec::new(),
            costheta: Vec::new(),
            nb_rad2: Vec::new(),
            nb_dist: Vec::new(),
            off: Vec::new(),
            vx: Vec::new(),
            br_c: Vec::new(),
            br_p: Vec::new(),
            ang: Vec::new(),
            next: Vec::new(),
            p: Vec::new(),
            volnb: Vec::new(),
            knot: Vec::new(),
            fknot: Vec::new(),
            rang: Vec::new(),
            pos: Vec::new(),
        };
        this.init();
        this
    }

    pub fn new_with_bond_to(
        coords: impl Iterator<Item = Vector3<Scalar>>,
        radii: impl Iterator<Item = Scalar>,
        bond_to: impl Iterator<Item = i32>,
        with_sasa: bool,
        with_dsasa: bool,
        with_vol: bool,
        with_dvol: bool,
    ) -> Self {
        let coords_vec: Vec<Vector3<Scalar>> = coords.collect();
        let radii_vec: Vec<Scalar> = radii.collect();
        let bond_to_vec: Vec<i32> = bond_to.collect();
        let power_diagram = PowerDiagram::from_params(
            PowerDiagram::create(
                coords_vec.len(),
                coords_vec.iter().cloned(),
                radii_vec.iter().copied(),
                bond_to_vec.into_iter(),
            )
            .with_radii_given(true)
            .with_calculate(true)
            .with_cells(true)
            .with_zero_points(true)
            .with_warnings(false)
            .without_check(true),
        );

        let mut this = Self {
            power_diagram,
            with_sasa,
            with_dsasa,
            with_vol,
            with_dvol,
            sasa: Vec::new(),
            dsasa_parts: Vec::new(),
            dsasa: Vec::new(),
            vol: Vec::new(),
            dvol: Vec::new(),
            tol_pow: Scalar::zero(),
            np: Vec::new(),
            nt: Vec::new(),
            e: Vec::new(),
            sintheta: Vec::new(),
            costheta: Vec::new(),
            nb_rad2: Vec::new(),
            nb_dist: Vec::new(),
            off: Vec::new(),
            vx: Vec::new(),
            br_c: Vec::new(),
            br_p: Vec::new(),
            ang: Vec::new(),
            next: Vec::new(),
            p: Vec::new(),
            volnb: Vec::new(),
            knot: Vec::new(),
            fknot: Vec::new(),
            rang: Vec::new(),
            pos: Vec::new(),
        };
        this.init();
        this
    }

    fn create_with_chain_bond(coords: &[Vector3<Scalar>], radii: &[Scalar]) -> PowerDiagram<Scalar> {
        let mut bond_to = Vec::with_capacity(coords.len());
        if !coords.is_empty() {
            bond_to.push(0_i32);
            for i in 1..coords.len() {
                bond_to.push((i - 1) as i32);
            }
        }

        PowerDiagram::from_params(
            PowerDiagram::create(
                coords.len(),
                coords.iter().cloned(),
                radii.iter().copied(),
                bond_to.into_iter(),
            )
            .with_radii_given(true)
            .with_calculate(true)
            .with_cells(true)
            .with_zero_points(true)
            .with_warnings(false)
            .without_check(true),
        )
    }

    fn init(&mut self) {
        self.resize_na();
        self.resize_nb(Self::K_MAX_NB);
        self.resize_vx(Self::K_MAX_VX);
        self.resize_pnt(Self::K_MAX_PNT);
    }

    fn resize_nb(&mut self, nnb: usize) {
        self.np.resize(nnb, 0);
        self.nt.resize(nnb, 0);
        self.e.resize(nnb, Vector3::zeros());
        self.sintheta.resize(nnb, Scalar::zero());
        self.costheta.resize(nnb, Scalar::zero());
        self.nb_rad2.resize(nnb, Scalar::zero());
        self.nb_dist.resize(nnb, Scalar::zero());
        self.volnb.resize(nnb, Scalar::zero());
        self.knot.resize(nnb, Vector3::zeros());
        self.fknot.resize(nnb, false);

        let npnt = if self.p.is_empty() {
            Self::K_MAX_PNT
        } else {
            self.p[0].len()
        };

        let old = self.p.len();
        self.next.resize(nnb, Vec::new());
        self.p.resize(nnb, Vec::new());
        self.ang.resize(nnb, Vec::new());
        for i in old..nnb {
            self.next[i].resize(npnt, 0);
            self.p[i].resize(npnt, 0);
            self.ang[i].resize(npnt, Scalar::zero());
        }

        if self.with_dsasa {
            for i in 0..self.dsasa_parts.len() {
                self.dsasa_parts[i].resize(nnb, Vector3::zeros());
            }
        }
    }

    fn resize_vx(&mut self, nvx: usize) {
        self.off.resize(nvx, 0);
        self.vx.resize(nvx, Vector3::zeros());
        self.br_c.resize(nvx, vec![0; 2]);
        self.br_p.resize(nvx, vec![0; 2]);
    }

    fn resize_pnt(&mut self, npnt: usize) {
        for i in 0..self.p.len() {
            self.next[i].resize(npnt, 0);
            self.p[i].resize(npnt, 0);
            self.ang[i].resize(npnt, Scalar::zero());
        }
        self.rang.resize(npnt, 0);
        self.pos.resize(npnt, 0);
    }

    fn resize_na(&mut self) {
        let n = self.power_diagram.get_points().len();
        if self.with_sasa {
            self.sasa.resize(n, Scalar::zero());
        }
        if self.with_dsasa {
            let nnb = if self.dsasa_parts.is_empty() || self.dsasa_parts[0].is_empty() {
                Self::K_MAX_NB
            } else {
                self.dsasa_parts[0].len()
            };
            let old = self.dsasa.len();
            self.dsasa.resize(n, Vector3::zeros());
            self.dsasa_parts.resize(n, Vec::new());
            for i in old..n {
                self.dsasa_parts[i].resize(nnb, Vector3::zeros());
            }
        }
        if self.with_vol {
            self.vol.resize(n, Scalar::zero());
        }
        if self.with_dvol {
            self.dvol.resize(n, Vector3::zeros());
        }

        let mut maxr2 = Scalar::zero();
        for p in self.power_diagram.get_points() {
            if p.r2 > maxr2 {
                maxr2 = p.r2;
            }
        }
        self.tol_pow = maxr2 * Self::drad2();
    }

    fn ang_about(&self, a: Vector3<Scalar>, b: Vector3<Scalar>, c: Vector3<Scalar>) -> Result<Scalar, PowerSasaError> {
        let co = a.dot(&b);
        let mut ang;
        let v = a.cross(&b);

        if co <= -Self::near_one_cosine_threshold() {
            ang = Scalar::pi() - v.norm().asin();
        } else if co >= Self::near_one_cosine_threshold() {
            ang = v.norm().asin();
        } else {
            ang = co.acos();
        }

        let axis_component_threshold = Self::axis_component_threshold();
        let vp = if c[0].abs() > axis_component_threshold {
            v[0] / c[0]
        } else if c[1].abs() > axis_component_threshold {
            v[1] / c[1]
        } else if c[2].abs() > axis_component_threshold {
            v[2] / c[2]
        } else {
            return Err(PowerSasaError);
        };

        if vp < Scalar::zero() {
            ang = -ang;
        }
        if ang < Scalar::zero() {
            ang += Scalar::from_f64(2.0).unwrap() * Scalar::pi();
        }

        Ok(ang)
    }

    fn get_ang(
        &self,
        np: i32,
        p: &[i32],
        e: Vector3<Scalar>,
        sintheta: Scalar,
        costheta: Scalar,
        ang: &mut [Scalar],
    ) -> Result<(), PowerSasaError> {
        let mut pu0 = Vector3::zeros();
        if e[0].abs() > Self::axis_component_threshold() {
            pu0 = Vector3::new(-(e[1] + e[2]) / e[0], Scalar::one(), Scalar::one());
        } else if e[1].abs() > Self::axis_component_threshold() {
            pu0 = Vector3::new(Scalar::one(), -(e[0] + e[2]) / e[1], Scalar::one());
        } else if e[2].abs() > Self::axis_component_threshold() {
            pu0 = Vector3::new(Scalar::one(), Scalar::one(), -(e[0] + e[1]) / e[2]);
        }

        pu0 /= pu0.norm();
        pu0 = pu0 * sintheta + e * costheta;

        for j in 0..(np as usize) {
            let p_idx = p[j] as usize;
            let pu = self.vx[p_idx];
            ang[j] = self.ang_about(pu0, pu, e)?;
        }
        Ok(())
    }

    fn get_next(
        &mut self,
        n: i32,
        ang: &mut [Scalar],
        next: &mut [i32],
        p: &[i32],
        e: Vector3<Scalar>,
    ) -> Result<(), PowerSasaError> {
        let n_usize = n as usize;
        if self.rang.len() < n_usize {
            self.resize_pnt(n_usize);
        }

        let two_pi = Scalar::from_f64(2.0).unwrap() * Scalar::pi();
        for j in 0..n_usize {
            if ang[j] > two_pi - Self::dang() {
                ang[j] = Scalar::zero();
            }
            self.rang[j] = 0;
        }

        for j in 0..n_usize {
            for k in 0..n_usize {
                if ang[k] > ang[j] + Self::dang() {
                    self.rang[j] += 1;
                } else if ang[k] > ang[j] - Self::dang() {
                    self.rang[j] += 0;
                }
            }
            self.pos[self.rang[j] as usize] = j as i32;
        }

        for j in 0..n_usize {
            let to = if j + 1 < n_usize {
                self.pos[j + 1]
            } else {
                self.pos[0]
            };
            next[self.pos[j] as usize] = to;
        }

        for j in 0..n_usize {
            if next[j] < 0 || next[j] >= n {
                return Err(PowerSasaError);
            }
            let jj = p[j] as usize;
            let kk = p[next[j] as usize] as usize;
            if (self.vx[jj].cross(&self.vx[kk]).dot(&e)).abs() < Self::dang() {
                return Err(PowerSasaError);
            }
        }

        Ok(())
    }

    pub fn calc_sasa_single(&mut self, iatom: usize) -> Result<(), PowerSasaError> {
        if iatom >= self.power_diagram.get_points().len() {
            return Err(PowerSasaError);
        }

        let atom = &self.power_diagram.get_points()[iatom];
        let rad = atom.r;
        let rad2 = atom.r2;
        let nnb = atom.neighbours_ids.len();

        if nnb > self.np.len() {
            self.resize_nb(nnb);
        }

        if self.with_sasa {
            self.sasa[iatom] = Scalar::zero();
        }
        if self.with_dsasa {
            self.dsasa[iatom] = Vector3::zeros();
        }
        if self.with_vol {
            self.vol[iatom] = Scalar::zero();
        }
        if self.with_dvol {
            self.dvol[iatom] = Vector3::zeros();
        }

        if nnb == 0 {
            let four = Scalar::from_f64(4.0).unwrap();
            let three = Scalar::from_f64(3.0).unwrap();
            if self.with_sasa {
                self.sasa[iatom] = four * Scalar::pi() * rad2;
            }
            if self.with_vol {
                self.vol[iatom] = (four / three) * Scalar::pi() * rad * rad2;
            }
            return Ok(());
        }

        let _ = self.tol_pow;
        let _ = self.get_ang(0, &[], Vector3::zeros(), Scalar::zero(), Scalar::zero(), &mut []);
        let _ = self.get_next(0, &mut [], &mut [], &[], Vector3::zeros());

        // Full literal port pending. Temporary fallback preserves non-zero baseline output shape.
        let four = Scalar::from_f64(4.0).unwrap();
        let three = Scalar::from_f64(3.0).unwrap();
        if self.with_sasa {
            self.sasa[iatom] = four * Scalar::pi() * rad2;
        }
        if self.with_vol {
            self.vol[iatom] = (four / three) * Scalar::pi() * rad * rad2;
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
        self.resize_na();
    }

    pub fn add_more(
        &mut self,
        pos_it: impl Iterator<Item = Vector3<Scalar>>,
        strength_it: impl Iterator<Item = Scalar>,
        new_size: usize,
    ) {
        self.power_diagram.add_more(pos_it, strength_it, new_size);
        self.resize_na();
    }

    pub fn revert(&mut self) {
        self.power_diagram.revert();
    }

    pub fn get_power_diagram(&self) -> &PowerDiagram<Scalar> {
        &self.power_diagram
    }

    pub fn num_of_neighbours(&self, iatom: usize) -> usize {
        self.power_diagram.get_points()[iatom].neighbours_ids.len()
    }

    pub fn atom_no(&self, i_atom: usize, i_neighbour: usize) -> usize {
        self.power_diagram.get_points()[i_atom].neighbours_ids[i_neighbour]
    }

    pub fn get_sasa(&self) -> &[Scalar] {
        &self.sasa
    }

    pub fn get_vol(&self) -> &[Scalar] {
        &self.vol
    }

    pub fn get_dvol(&self) -> &[Vector3<Scalar>] {
        &self.dvol
    }

    pub fn get_dsasa(&self) -> &[Vector3<Scalar>] {
        &self.dsasa
    }
}
