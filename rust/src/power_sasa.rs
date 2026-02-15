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

        let (rad, rad2, atom_pos, neighbours_ids, my_vertices_ids) = {
            let atom = &self.power_diagram.get_points()[iatom];
            (
                atom.r,
                atom.r2,
                atom.position,
                atom.neighbours_ids.clone(),
                atom.my_vertices_ids.clone(),
            )
        };
        let nnb = neighbours_ids.len();

        if nnb > self.np.len() {
            self.resize_nb(nnb);
        }

        if self.with_sasa {
            self.sasa[iatom] = Scalar::zero();
        }
        if self.with_dsasa {
            self.dsasa[iatom] = Vector3::zeros();
            for i in 0..nnb {
                self.dsasa_parts[iatom][i] = Vector3::zeros();
            }
        }
        if self.with_vol {
            self.vol[iatom] = Scalar::zero();
        }
        if self.with_dvol {
            self.dvol[iatom] = Vector3::zeros();
        }

        let do_sasa = self.with_sasa || self.with_vol;

        if nnb == 0 {
            let four = Scalar::from_f64(4.0).unwrap();
            let three = Scalar::from_f64(3.0).unwrap();
            let mut is_owner = false;
            if let Some(first_vertex) = self.power_diagram.get_vertices().first() {
                let gref = first_vertex.generator_refs[0];
                if gref.is_valid() && gref.kind == crate::power_diagram::GeneratorKind::Point && gref.index == iatom {
                    is_owner = true;
                }
            }
            if is_owner {
                if self.with_sasa {
                    self.sasa[iatom] = four * Scalar::pi() * rad2;
                }
                if self.with_vol {
                    self.vol[iatom] = (four / three) * Scalar::pi() * rad * rad2;
                }
            }
            return Ok(());
        }

        let mut covered = true;
        for vid in &my_vertices_ids {
            if *vid >= self.power_diagram.get_vertices().len() {
                continue;
            }
            let atom_vertex = &self.power_diagram.get_vertices()[*vid];
            if atom_vertex.power_value.abs() < self.tol_pow {
                return Err(PowerSasaError);
            }
            if atom_vertex.power_value > Scalar::zero() {
                covered = false;
            }
        }
        if covered && !self.with_vol {
            return Ok(());
        }

        for (i, nb_id) in neighbours_ids.iter().copied().enumerate() {
            if let Some(nb) = self.power_diagram.get_points_mut().get_mut(nb_id) {
                nb.visited_as = i as i32;
            } else {
                return Err(PowerSasaError);
            }
        }

        let mut n_apart = 0_usize;
        let mut contributing = 0_usize;
        for (i, nb_id) in neighbours_ids.iter().copied().enumerate() {
            let neighbour = &self.power_diagram.get_points()[nb_id];
            let nb_rad = neighbour.r;
            let nb_rad2 = neighbour.r2;
            self.nb_rad2[i] = nb_rad2;
            if self.with_vol {
                self.volnb[i] = Scalar::zero();
                self.fknot[i] = false;
            }
            self.np[i] = 0;
            self.nt[i] = 0;

            let rel_pos = neighbour.position - atom_pos;
            let dist = rel_pos.norm();
            self.nb_dist[i] = dist;

            if dist <= nb_rad - rad {
                // Completely covered by one larger neighbor.
                if self.with_sasa {
                    self.sasa[iatom] = Scalar::zero();
                }
                if self.with_vol {
                    self.vol[iatom] = Scalar::zero();
                }
                for nb_id in neighbours_ids {
                    if let Some(nb) = self.power_diagram.get_points_mut().get_mut(nb_id) {
                        nb.visited_as = 0;
                    }
                }
                return Ok(());
            }

            if dist >= rad + nb_rad || dist <= rad - nb_rad {
                n_apart += 1;
                self.np[i] = -1;
                continue;
            }

            self.costheta[i] = (Scalar::from_f64(0.5).unwrap() * (dist + (rad2 - nb_rad2) / dist)) / rad;
            self.sintheta[i] = (Scalar::one() - self.costheta[i] * self.costheta[i]).sqrt();
            self.e[i] = rel_pos / dist;
            contributing += 1;
        }

        if n_apart == nnb || contributing == 0 {
            let four = Scalar::from_f64(4.0).unwrap();
            let three = Scalar::from_f64(3.0).unwrap();
            if self.with_sasa {
                self.sasa[iatom] = four * Scalar::pi() * rad2;
            }
            if self.with_vol {
                self.vol[iatom] = (four / three) * Scalar::pi() * rad * rad2;
            }
            for nb_id in neighbours_ids {
                if let Some(nb) = self.power_diagram.get_points_mut().get_mut(nb_id) {
                    nb.visited_as = 0;
                }
            }
            return Ok(());
        }

        let mut nvx = 0_usize;
        let mut partner = [0_i32; 2];
        for zp_i in self.power_diagram.get_points()[iatom].my_zero_points.clone() {
            if zp_i < 0 {
                continue;
            }
            let zp_idx = zp_i as usize;
            if zp_idx >= self.power_diagram.get_zero_points().len() {
                continue;
            }
            let zp = &self.power_diagram.get_zero_points()[zp_idx];
            if !self.power_diagram.zero_point_valid(zp) {
                continue;
            }
            let zp_pos = self.power_diagram.zero_point_pos(zp);
            let mut ptn = 0_usize;
            for kg in 0..3 {
                let zp_generator_ref = zp.generator_refs[kg];
                if !zp_generator_ref.is_valid() {
                    continue;
                }
                if zp_generator_ref.kind != crate::power_diagram::GeneratorKind::Point || zp_generator_ref.index != iatom {
                    partner[ptn] = self.power_diagram.get_generator(zp_generator_ref).visited_as;
                    ptn += 1;
                }
            }
            if ptn != 2 {
                return Err(PowerSasaError);
            }
            let ptn0 = partner[0] as usize;
            let ptn1 = partner[1] as usize;
            if ptn0 >= nnb || ptn1 >= nnb {
                return Err(PowerSasaError);
            }

            if zp.pos < Scalar::zero() || zp.pos > Scalar::one() {
                self.nt[ptn0] += 1;
                self.nt[ptn1] += 1;
                continue;
            }

            if self.np[ptn0] < 0 || self.np[ptn1] < 0 {
                return Err(PowerSasaError);
            }
            if (self.np[ptn0] as usize) >= self.rang.len() || (self.np[ptn1] as usize) >= self.rang.len() {
                let new_pnt = if self.np[ptn0] > self.np[ptn1] { self.np[ptn0] + 1 } else { self.np[ptn1] + 1 };
                self.resize_pnt(new_pnt as usize);
            }
            if nvx >= self.vx.len() {
                self.resize_vx(nvx + 1);
            }

            self.vx[nvx] = (zp_pos - atom_pos) / rad;
            self.p[ptn0][self.np[ptn0] as usize] = nvx as i32;
            self.p[ptn1][self.np[ptn1] as usize] = nvx as i32;

            self.br_c[nvx][0] = ptn0 as i32;
            self.br_c[nvx][1] = ptn1 as i32;
            self.br_p[nvx][0] = self.np[ptn0];
            self.br_p[nvx][1] = self.np[ptn1];
            nvx += 1;
            self.np[ptn0] += 1;
            self.np[ptn1] += 1;
        }

        for i in 0..nnb {
            if self.np[i] <= 0 {
                continue;
            }
            if self.np[i] % 2 != 0 {
                return Err(PowerSasaError);
            }
            let np_i = self.np[i] as usize;
            let mut ang_i = vec![Scalar::zero(); np_i];
            let p_i = self.p[i][..np_i].to_vec();
            self.get_ang(
                self.np[i],
                &p_i,
                self.e[i],
                self.sintheta[i],
                self.costheta[i],
                &mut ang_i,
            )?;
            let mut next_i = vec![0_i32; np_i];
            self.get_next(self.np[i], &mut ang_i, &mut next_i, &p_i, self.e[i])?;
            self.ang[i][..np_i].copy_from_slice(&ang_i);
            self.next[i][..np_i].copy_from_slice(&next_i);
        }

        let two = Scalar::from_f64(2.0).unwrap();
        let three = Scalar::from_f64(3.0).unwrap();
        let four = Scalar::from_f64(4.0).unwrap();
        let six = Scalar::from_f64(6.0).unwrap();
        let half = Scalar::from_f64(0.5).unwrap();
        let two_pi = two * Scalar::pi();

        let mut vol2 = Scalar::zero();
        let mut vol3 = Scalar::zero();
        let mut sasa_ia = Scalar::zero();

        for iv in 0..nvx {
            self.off[iv] = 0;
        }

        for iv in 0..nvx {
            if self.off[iv] != 0 {
                continue;
            }

            let p_ini_idx = iv;
            let ic_0 = self.br_c[iv][0] as usize;
            let ic_1 = self.br_c[iv][1] as usize;

            let dirdet = self.e[ic_1].cross(&self.e[ic_0]).dot(&self.vx[p_ini_idx]);
            if dirdet == Scalar::zero() {
                return Err(PowerSasaError);
            }

            let (mut ic1, mut ic2, mut ip2): (usize, usize, i32) = if dirdet > Scalar::zero() {
                (ic_0, ic_1, self.br_p[iv][1])
            } else {
                (ic_1, ic_0, self.br_p[iv][0])
            };

            let mut pt_idx = p_ini_idx;
            if do_sasa {
                sasa_ia += two_pi;
            }
            let mut count = 0_usize;

            loop {
                count += 1;
                if count > Self::K_MAX_COUNT {
                    return Err(PowerSasaError);
                }

                let ip_next = self.next[ic2][ip2 as usize];
                let mut phi = self.ang[ic2][ip_next as usize] - self.ang[ic2][ip2 as usize];
                if phi < Scalar::zero() {
                    phi += two_pi;
                }

                if do_sasa {
                    let mut co = (self.e[ic1].dot(&self.e[ic2]) - self.costheta[ic1] * self.costheta[ic2])
                        / (self.sintheta[ic1] * self.sintheta[ic2]);
                    co = Self::clamp_unit_interval(co);
                    sasa_ia += phi * self.costheta[ic2] - co.acos();
                }

                self.off[self.p[ic2][ip2 as usize] as usize] = 1;
                ic1 = ic2;

                let ivx = self.p[ic1][ip_next as usize] as usize;
                let pt0_idx = pt_idx;
                pt_idx = ivx;

                if self.br_c[ivx][0] == ic1 as i32 {
                    ic2 = self.br_c[ivx][1] as usize;
                    ip2 = self.br_p[ivx][1];
                } else {
                    ic2 = self.br_c[ivx][0] as usize;
                    ip2 = self.br_p[ivx][0];
                }

                if self.with_dsasa {
                    let ds1 = half
                        * rad
                        * phi
                        * (Scalar::one() + (self.nb_rad2[ic1] - rad2) / (self.nb_dist[ic1] * self.nb_dist[ic1]));
                    let ds2 = rad2 / self.nb_dist[ic1];
                    self.dsasa_parts[iatom][ic1] += self.e[ic1] * ds1
                        - (self.vx[pt_idx] - self.vx[pt0_idx]).cross(&self.e[ic1]) * ds2;
                }

                let mut scone = Scalar::zero();
                let mut vv = Vector3::zeros();
                if self.with_vol || self.with_dvol {
                    vv = (atom_pos + self.vx[pt0_idx] * rad).cross(&(atom_pos + self.vx[pt_idx] * rad));
                    scone = self.sintheta[ic1] * self.sintheta[ic1] * (phi - phi.sin());
                }

                if self.with_vol {
                    vol2 += self.costheta[ic1] * scone;
                    if !self.fknot[ic1] {
                        self.fknot[ic1] = true;
                        self.knot[ic1] = atom_pos + self.vx[pt0_idx] * rad;
                    } else {
                        let a = atom_pos + self.vx[pt0_idx] * rad - self.knot[ic1];
                        let b = atom_pos + self.vx[pt_idx] * rad - self.knot[ic1];
                        self.volnb[ic1] += a.cross(&b).dot(&self.e[ic1]).abs();
                    }
                }

                if self.with_dvol {
                    self.dvol[iatom] -= (vv + self.e[ic1] * (rad2 * scone)) * half;
                }

                if pt_idx == p_ini_idx {
                    break;
                }
            }

            if do_sasa && sasa_ia > four * Scalar::pi() {
                sasa_ia -= four * Scalar::pi();
            }
        }

        for i in 0..nnb {
            if self.np[i] != 0 || self.nt[i] != 0 {
                continue;
            }

            let mut ok = true;
            let cc = atom_pos + self.e[i] * (rad * self.costheta[i]);
            let pw_i = -self.sintheta[i] * self.sintheta[i] * rad2;

            for j in 0..nnb {
                if j == i {
                    continue;
                }
                let nb_j = self.power_diagram.get_points()[neighbours_ids[j]].position;
                let pw_j = (nb_j - cc).norm_squared() - self.nb_rad2[j];
                if pw_j <= pw_i {
                    ok = false;
                    break;
                }
            }

            if ok {
                if do_sasa {
                    sasa_ia += two_pi * (Scalar::one() + self.costheta[i]);
                    if sasa_ia > four * Scalar::pi() {
                        sasa_ia -= four * Scalar::pi();
                    }
                }

                if self.with_dsasa {
                    self.dsasa_parts[iatom][i] = self.e[i]
                        * (rad
                            * Scalar::pi()
                            * (Scalar::one() + (self.nb_rad2[i] - rad2) / (self.nb_dist[i] * self.nb_dist[i])));
                }

                let mut scone = Scalar::zero();
                if self.with_vol || self.with_dvol {
                    scone = self.sintheta[i] * self.sintheta[i] * two_pi;
                }
                if self.with_vol {
                    vol2 += self.costheta[i] * scone;
                }
                if self.with_dvol {
                    self.dvol[iatom] -= self.e[i] * (half * rad2 * scone);
                }
            }
        }

        if self.with_sasa {
            self.sasa[iatom] = rad2 * sasa_ia;
        }
        if self.with_dsasa {
            for i in 0..nnb {
                self.dsasa[iatom] -= self.dsasa_parts[iatom][i];
            }
        }
        if self.with_vol {
            for i in 0..nnb {
                if self.fknot[i] {
                    vol3 += rad * self.volnb[i] * self.costheta[i];
                }
            }
            self.vol[iatom] = rad * rad2 * sasa_ia / three + rad * rad2 * vol2 / six + vol3 / six;
        }

        for nb_id in neighbours_ids {
            if let Some(nb) = self.power_diagram.get_points_mut().get_mut(nb_id) {
                nb.visited_as = 0;
            }
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
