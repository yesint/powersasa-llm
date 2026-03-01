use nalgebra::{RealField, Vector3};

use crate::error::SasaError;
use crate::power_diagram::PowerDiagram;

/// PowerSasa orchestrates per-atom SASA/volume evaluation over a power diagram, including optional derivatives.
pub struct PowerSasa<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Weighted power diagram over solvent-expanded atomic spheres.
    power_diagram: PowerDiagram<Scalar>,

    /// Enable per-atom SASA computation.
    with_sasa: bool,
    /// Enable per-atom SASA gradient (`dsasa`) computation.
    with_dsasa: bool,
    /// Enable per-atom volume computation.
    with_vol: bool,
    /// Enable per-atom volume gradient (`dvol`) computation.
    with_dvol: bool,

    /// Per-atom solvent-accessible surface area.
    sasa: Vec<Scalar>,
    /// Per-neighbor partial SASA-gradient contributions used to assemble `dsasa` (scratch, reset per atom).
    dsasa_parts: Vec<Vector3<Scalar>>,
    /// Final per-atom SASA gradient after summing neighbor contributions.
    dsasa: Vec<Vector3<Scalar>>,
    /// Per-atom volume contribution.
    vol: Vec<Scalar>,
    /// Per-atom volume gradient.
    dvol: Vec<Vector3<Scalar>>,

    /// Tolerance for near-zero power values when classifying local topology.
    tol_pow: Scalar,

    /// Number of contour points registered on each neighbor circle.
    np: Vec<i32>,
    /// Marker/counter for "single circle" handling on each neighbor.
    nt: Vec<i32>,
    /// Unit direction from current atom center to each neighbor center.
    e: Vec<Vector3<Scalar>>,
    /// `sin(theta)` for each atom-neighbor intersection circle.
    sintheta: Vec<Scalar>,
    /// `cos(theta)` for each atom-neighbor intersection circle.
    costheta: Vec<Scalar>,
    /// Squared solvent-expanded neighbor radii.
    nb_rad2: Vec<Scalar>,
    /// Center-to-center distance from current atom to each neighbor.
    nb_dist: Vec<Scalar>,

    /// Per-contour-vertex visitation/ownership flag during loop traversal.
    off: Vec<i32>,
    /// Unit vectors of contour vertices on the current atom sphere.
    vx: Vec<Vector3<Scalar>>,
    /// For each contour vertex, the two incident neighbor-circle ids.
    br_c: Vec<[i32; 2]>,
    /// For each contour vertex, local point indices on the two incident circles.
    br_p: Vec<[i32; 2]>,

    /// Polar angle of each registered point on each neighbor circle.
    ang: Vec<Vec<Scalar>>,
    /// Successor index for ordered traversal around each neighbor circle.
    next: Vec<Vec<i32>>,
    /// Mapping from per-circle point slot to global contour-vertex index.
    p: Vec<Vec<i32>>,

    /// Per-neighbor volume accumulator from triangulated contour sectors.
    volnb: Vec<Scalar>,
    /// First anchor point used to build fan triangles for each neighbor volume term.
    knot: Vec<Vector3<Scalar>>,
    /// Indicates whether `knot` has been initialized for a neighbor.
    fknot: Vec<bool>,

    /// Temporary rank buffer used to sort circle points by angle.
    rang: Vec<i32>,
    /// Temporary position/index buffer paired with `rang` for angle ordering.
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

    #[inline(always)]
    /// Returns the power-value tolerance scale used to detect near-degenerate radius/power comparisons.
    pub(crate) fn drad2() -> Scalar {
        Scalar::from_f64(1000.0).unwrap() * Scalar::default_epsilon()
    }

    #[inline(always)]
    /// Returns the angular tolerance used when sorting contour points on a neighbor circle.
    pub(crate) fn dang() -> Scalar {
        Scalar::from_f64(1000.0).unwrap() * Scalar::default_epsilon()
    }

    #[inline(always)]
    /// Threshold for stable angle evaluation when cosine is very close to +-1.
    pub(crate) fn near_one_cosine_threshold() -> Scalar {
        Scalar::from_f64(0.999).unwrap()
    }

    #[inline(always)]
    /// Minimum axis component magnitude used to pick a numerically stable orientation divisor.
    pub(crate) fn axis_component_threshold() -> Scalar {
        Scalar::from_f64(0.001).unwrap()
    }

    /// Builds a PowerSasa with default chain bond hints and initializes all scratch/output buffers.
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

    /// Builds a PowerSasa with explicit bond-parent hints for power-diagram insertion order.
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

    /// Creates a PowerDiagram using a simple i->i-1 bond chain when no bond map is supplied.
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

    /// Initializes atom-, neighbor-, vertex-, and point-ordering work arrays to default capacities.
    fn init(&mut self) {
        self.resize_nb(Self::K_MAX_NB);
        self.resize_vx(Self::K_MAX_VX);
        self.resize_pnt(Self::K_MAX_PNT);
        self.update_tol_pow();
    }

    /// Resizes all neighbor-scoped scratch buffers used by per-atom contour integration.
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
            self.dsasa_parts.resize(nnb, Vector3::zeros());
        }
    }

    /// Resizes temporary contour-vertex storage and bridge tables linking circles through vertices.
    fn resize_vx(&mut self, nvx: usize) {
        self.off.resize(nvx, 0);
        self.vx.resize(nvx, Vector3::zeros());
        self.br_c.resize(nvx, [0; 2]);
        self.br_p.resize(nvx, [0; 2]);
    }

    /// Resizes per-circle point-order arrays used to sort angles and build successor links.
    fn resize_pnt(&mut self, npnt: usize) {
        for i in 0..self.p.len() {
            self.next[i].resize(npnt, 0);
            self.p[i].resize(npnt, 0);
            self.ang[i].resize(npnt, Scalar::zero());
        }
        self.rang.resize(npnt, 0);
        self.pos.resize(npnt, 0);
    }

    /// Refreshes the power ambiguity tolerance from the current maximum atom radius.
    fn update_tol_pow(&mut self) {
        let maxr2 = self.power_diagram.points.iter()
            .map(|p| p.r2)
            .fold(Scalar::zero(), |a, b| if b > a { b } else { a });
        self.tol_pow = maxr2 * Self::drad2();
    }

    /// Computes the oriented angle from vector a to b around axis c, normalized to [0, 2*pi).
    fn ang_about(a: Vector3<Scalar>, b: Vector3<Scalar>, c: Vector3<Scalar>) -> Result<Scalar, SasaError> {
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
            return Err(SasaError::AxisTooShort);
        };

        if vp < Scalar::zero() {
            ang = -ang;
        }
        if ang < Scalar::zero() {
            ang += Scalar::from_f64(2.0).unwrap() * Scalar::pi();
        }

        Ok(ang)
    }

    /// Computes polar angles of contour points projected onto one neighbor intersection circle.
    fn get_ang(
        np: i32,
        p: &[i32],
        e: Vector3<Scalar>,
        sintheta: Scalar,
        costheta: Scalar,
        ang: &mut [Scalar],
        vx: &[Vector3<Scalar>],
    ) -> Result<(), SasaError> {
        if np <= 0 {
            return Ok(());
        }
        let n = np as usize;
        if n > p.len() || n > ang.len() {
            return Err(SasaError::InvalidGetAngInput);
        }
        ang[0] = Scalar::zero();
        let p0 = p[0] as usize;
        if p0 >= vx.len() {
            return Err(SasaError::SurfaceVertexIndexOutOfBounds);
        }
        let pu0 = (vx[p0] - e * costheta) / sintheta;
        for j in 1..n {
            let pj = p[j] as usize;
            if pj >= vx.len() {
                return Err(SasaError::SurfaceVertexIndexOutOfBounds);
            }
            let pu = (vx[pj] - e * costheta) / sintheta;
            ang[j] = Self::ang_about(pu0, pu, e)?;
        }
        Ok(())
    }

    /// Sorts circle points by angle (with tolerances) and builds cyclic next-point traversal links.
    fn get_next(
        n: i32,
        ang: &mut [Scalar],
        next: &mut [i32],
        p: &[i32],
        e: Vector3<Scalar>,
        vx: &[Vector3<Scalar>],
        rang: &mut [i32],
        pos: &mut [i32],
    ) -> Result<(), SasaError> {
        if n <= 0 {
            return Ok(());
        }
        let n_usize = n as usize;
        if n_usize > ang.len() || n_usize > next.len() || n_usize > p.len() {
            return Err(SasaError::InvalidGetNextInput);
        }
        let two_pi = Scalar::from_f64(2.0).unwrap() * Scalar::pi();

        for j in 1..n_usize {
            if ang[j] <= two_pi - Self::dang() {
                continue;
            }
            let pj = p[j] as usize;
            let p0 = p[0] as usize;
            if pj >= vx.len() || p0 >= vx.len() {
                return Err(SasaError::SurfaceVertexIndexOutOfBounds);
            }
            let dp = vx[pj].cross(&vx[p0]).dot(&e);
            if dp < Scalar::zero() {
                ang[j] = Scalar::zero();
            } else if dp == Scalar::zero() {
                return Err(SasaError::DegenerateCircleOrdering);
            }
        }

        rang[0] = 0;
        if n_usize > 1 {
            rang[1] = 1;
        }
        for j in 2..n_usize {
            let mut m = j as i32;
            for k in 1..j {
                if ang[k] > ang[j] + Self::dang() {
                    rang[k] += 1;
                    m -= 1;
                } else if ang[k] > ang[j] - Self::dang() {
                    let pj = p[j] as usize;
                    let pk = p[k] as usize;
                    if pj >= vx.len() || pk >= vx.len() {
                        return Err(SasaError::SurfaceVertexIndexOutOfBounds);
                    }
                    let dp = vx[pj].cross(&vx[pk]).dot(&e);
                    if dp > Scalar::zero() {
                        rang[k] += 1;
                        m -= 1;
                    } else if dp == Scalar::zero() {
                        return Err(SasaError::DegenerateCircleOrdering);
                    }
                }
            }
            rang[j] = m;
        }

        for j in 0..n_usize {
            pos[j] = -1;
        }
        for j in 0..n_usize {
            let rj = rang[j];
            if rj < 0 || (rj as usize) >= n_usize {
                return Err(SasaError::RankOutOfRange);
            }
            pos[rj as usize] = j as i32;
        }
        for j in 0..n_usize {
            let rj = rang[j];
            if rj == (n_usize as i32 - 1) {
                next[j] = pos[0];
            } else {
                next[j] = pos[(rj + 1) as usize];
            }
        }
        Ok(())
    }

    /// Evaluates SASA, volume, and optional derivatives for one atom from local power-diagram topology.
    fn calc_single(&mut self, iatom: usize) -> Result<(), SasaError> {
        if iatom >= self.power_diagram.points.len() {
            return Err(SasaError::AtomIndexOutOfBounds);
        }

        let (rad, rad2, atom_pos) = {
            let atom = &self.power_diagram.points[iatom];
            (atom.r, atom.r2, atom.position)
        };
        let nnb  = self.power_diagram.points[iatom].neighbours_ids.len();
        let nvid = self.power_diagram.points[iatom].my_vertices_ids.len();

        if nnb > self.np.len() {
            self.resize_nb(nnb);
        }

        let mut local_sasa = Scalar::zero();
        let mut local_vol = Scalar::zero();
        let mut local_dsasa: Vector3<Scalar> = Vector3::zeros();
        let mut local_dvol: Vector3<Scalar> = Vector3::zeros();

        if self.with_dsasa {
            for i in 0..nnb {
                self.dsasa_parts[i] = Vector3::zeros();
            }
        }

        let do_sasa = self.with_sasa || self.with_vol;

        if nnb == 0 {
            let four = Scalar::from_f64(4.0).unwrap();
            let three = Scalar::from_f64(3.0).unwrap();
            let mut is_owner = false;
            if let Some(first_vertex) = self.power_diagram.vertices.first() {
                let gref = first_vertex.generator_refs[0];
                if gref.is_valid() && gref.kind == crate::power_diagram::GeneratorKind::Point && gref.index == iatom {
                    is_owner = true;
                }
            }
            if is_owner {
                if self.with_sasa { local_sasa = four * Scalar::pi() * rad2; }
                if self.with_vol  { local_vol  = (four / three) * Scalar::pi() * rad * rad2; }
            }
            if self.with_sasa  { self.sasa.push(local_sasa); }
            if self.with_vol   { self.vol.push(local_vol);   }
            if self.with_dsasa { self.dsasa.push(local_dsasa); }
            if self.with_dvol  { self.dvol.push(local_dvol);  }
            return Ok(());
        }

        let mut covered = true;
        for vi in 0..nvid {
            let vid = self.power_diagram.points[iatom].my_vertices_ids[vi];
            if vid >= self.power_diagram.vertices.len() {
                continue;
            }
            let atom_vertex = &self.power_diagram.vertices[vid];
            if atom_vertex.power_value.abs() < self.tol_pow {
                return Err(SasaError::AmbiguousVertexPower);
            }
            if atom_vertex.power_value > Scalar::zero() {
                covered = false;
            }
        }
        if covered && !self.with_vol {
            if self.with_sasa  { self.sasa.push(local_sasa); }
            if self.with_vol   { self.vol.push(local_vol);   }
            if self.with_dsasa { self.dsasa.push(local_dsasa); }
            if self.with_dvol  { self.dvol.push(local_dvol);  }
            return Ok(());
        }

        for i in 0..nnb {
            let nb_id = self.power_diagram.points[iatom].neighbours_ids[i];
            if nb_id < self.power_diagram.points.len() {
                self.power_diagram.get_cell_mut(nb_id).visited_as = i as i32;
            } else {
                return Err(SasaError::NeighbourIndexOutOfBounds);
            }
        }

        let mut n_apart = 0_usize;
        let mut contributing = 0_usize;
        for i in 0..nnb {
            let nb_id = self.power_diagram.points[iatom].neighbours_ids[i];
            let neighbour = &self.power_diagram.points[nb_id];
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
                for k in 0..nnb {
                    let nb_k = self.power_diagram.points[iatom].neighbours_ids[k];
                    if nb_k < self.power_diagram.points.len() {
                        self.power_diagram.get_cell_mut(nb_k).visited_as = 0;
                    }
                }
                if self.with_sasa  { self.sasa.push(local_sasa); }
                if self.with_vol   { self.vol.push(local_vol);   }
                if self.with_dsasa { self.dsasa.push(local_dsasa); }
                if self.with_dvol  { self.dvol.push(local_dvol);  }
                return Ok(());
            }

            if dist >= rad + nb_rad || dist <= rad - nb_rad {
                n_apart += 1;
                self.np[i] = -1;
                continue;
            }

            self.costheta[i] = (Scalar::from_f64(0.5).unwrap() * (dist + (rad2 - nb_rad2) / dist)) / rad;
            if self.costheta[i] <= -Scalar::one() {
                for k in 0..nnb {
                    let nb_k = self.power_diagram.points[iatom].neighbours_ids[k];
                    if nb_k < self.power_diagram.points.len() {
                        self.power_diagram.get_cell_mut(nb_k).visited_as = 0;
                    }
                }
                if self.with_sasa  { self.sasa.push(local_sasa); }
                if self.with_vol   { self.vol.push(local_vol);   }
                if self.with_dsasa { self.dsasa.push(local_dsasa); }
                if self.with_dvol  { self.dvol.push(local_dvol);  }
                return Ok(());
            }
            if self.costheta[i] >= Scalar::one() {
                n_apart += 1;
                self.np[i] = -1;
                continue;
            }
            self.sintheta[i] = (Scalar::one() - self.costheta[i] * self.costheta[i]).sqrt();
            self.e[i] = rel_pos / dist;
            contributing += 1;
        }

        if n_apart == nnb || contributing == 0 {
            let four = Scalar::from_f64(4.0).unwrap();
            let three = Scalar::from_f64(3.0).unwrap();
            if self.with_sasa { local_sasa = four * Scalar::pi() * rad2; }
            if self.with_vol  { local_vol  = (four / three) * Scalar::pi() * rad * rad2; }
            for k in 0..nnb {
                let nb_k = self.power_diagram.points[iatom].neighbours_ids[k];
                if nb_k < self.power_diagram.points.len() {
                    self.power_diagram.get_cell_mut(nb_k).visited_as = 0;
                }
            }
            if self.with_sasa  { self.sasa.push(local_sasa); }
            if self.with_vol   { self.vol.push(local_vol);   }
            if self.with_dsasa { self.dsasa.push(local_dsasa); }
            if self.with_dvol  { self.dvol.push(local_dvol);  }
            return Ok(());
        }

        let mut nvx = 0_usize;
        let mut partner = [0_i32; 2];
        let nzp = self.power_diagram.points[iatom].my_zero_points.len();
        for zi in 0..nzp {
            let zp_i = self.power_diagram.points[iatom].my_zero_points[zi];
            if zp_i < 0 {
                continue;
            }
            let zp_idx = zp_i as usize;
            if zp_idx >= self.power_diagram.zeros.len() {
                continue;
            }
            if !self.power_diagram.zero_point_valid(&self.power_diagram.zeros[zp_idx]) {
                continue;
            }
            let zp_pos = self.power_diagram.zero_point_pos(&self.power_diagram.zeros[zp_idx]);
            let zp_gen_refs   = self.power_diagram.zeros[zp_idx].generator_refs;
            let zp_pos_scalar = self.power_diagram.zeros[zp_idx].pos;
            let zp_from_id    = self.power_diagram.zeros[zp_idx].from_id;
            let zp_branch     = self.power_diagram.zeros[zp_idx].branch;
            let mut ptn = 0_usize;
            for zp_generator_ref in zp_gen_refs {
                if !zp_generator_ref.is_valid() {
                    continue;
                }
                if zp_generator_ref.kind != crate::power_diagram::GeneratorKind::Point || zp_generator_ref.index != iatom {
                    partner[ptn] = self.power_diagram.get_generator(zp_generator_ref).visited_as;
                    ptn += 1;
                }
            }
            if ptn != 2 {
                return Err(SasaError::InvalidZeroPointPartnerCount);
            }
            let ptn0 = partner[0] as usize;
            let ptn1 = partner[1] as usize;
            if ptn0 >= nnb || ptn1 >= nnb {
                return Err(SasaError::PartnerOutOfRange);
            }

            if zp_pos_scalar < Scalar::zero() || zp_pos_scalar > Scalar::one() {
                self.nt[ptn0] += 1;
                self.nt[ptn1] += 1;
                continue;
            }

            if self.np[ptn0] < 0 || self.np[ptn1] < 0 {
                return Err(SasaError::InvalidPartnerPointCounts);
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

            if self.with_vol {
                let pd_vertices = &self.power_diagram.vertices;
                let node1_id = zp_from_id;
                if node1_id == crate::power_diagram::GeneratorRef::INVALID_ID || node1_id >= pd_vertices.len() {
                    return Err(SasaError::VolumeNodeOutOfBounds);
                }
                let node1 = &pd_vertices[node1_id];
                let node2_id = node1.end_point_ids[zp_branch as usize];
                if node2_id == crate::power_diagram::GeneratorRef::INVALID_ID || node2_id >= pd_vertices.len() {
                    return Err(SasaError::VolumeNodeOutOfBounds);
                }
                let node2 = &pd_vertices[node2_id];

                if node1.power_value < Scalar::zero() && node2.power_value > Scalar::zero() {
                    if !self.fknot[ptn0] {
                        self.fknot[ptn0] = true;
                        self.knot[ptn0] = node1.position;
                    } else {
                        self.volnb[ptn0] +=
                            (node1.position - self.knot[ptn0]).cross(&(zp_pos - self.knot[ptn0])).dot(&self.e[ptn0]).abs();
                    }
                    if !self.fknot[ptn1] {
                        self.fknot[ptn1] = true;
                        self.knot[ptn1] = node1.position;
                    } else {
                        self.volnb[ptn1] +=
                            (node1.position - self.knot[ptn1]).cross(&(zp_pos - self.knot[ptn1])).dot(&self.e[ptn1]).abs();
                    }
                } else if node1.power_value > Scalar::zero() && node2.power_value < Scalar::zero() {
                    if !self.fknot[ptn0] {
                        self.fknot[ptn0] = true;
                        self.knot[ptn0] = node2.position;
                    } else {
                        self.volnb[ptn0] +=
                            (node2.position - self.knot[ptn0]).cross(&(zp_pos - self.knot[ptn0])).dot(&self.e[ptn0]).abs();
                    }
                    if !self.fknot[ptn1] {
                        self.fknot[ptn1] = true;
                        self.knot[ptn1] = node2.position;
                    } else {
                        self.volnb[ptn1] +=
                            (node2.position - self.knot[ptn1]).cross(&(zp_pos - self.knot[ptn1])).dot(&self.e[ptn1]).abs();
                    }
                } else if node1.power_value > Scalar::zero() && node2.power_value > Scalar::zero() {
                    let denom =
                        node2.power_value * zp_pos_scalar + node1.power_value * (Scalar::one() - zp_pos_scalar);
                    if denom == Scalar::zero() {
                        return Err(SasaError::ZeroVolumeDenominator);
                    }
                    let dpos = node1.power_value * (Scalar::one() - zp_pos_scalar) / denom - zp_pos_scalar;
                    let half = Scalar::from_f64(0.5).unwrap();
                    if !self.fknot[ptn0] {
                        self.fknot[ptn0] = true;
                        self.knot[ptn0] = zp_pos;
                    } else {
                        self.volnb[ptn0] +=
                            (half * dpos * (zp_pos - self.knot[ptn0]).cross(&(node2.position - node1.position)).dot(&self.e[ptn0]))
                                .abs();
                    }
                    if !self.fknot[ptn1] {
                        self.fknot[ptn1] = true;
                        self.knot[ptn1] = zp_pos;
                    } else {
                        self.volnb[ptn1] +=
                            (half * dpos * (zp_pos - self.knot[ptn1]).cross(&(node2.position - node1.position)).dot(&self.e[ptn1]))
                                .abs();
                    }
                } else {
                    return Err(SasaError::InvalidVolumeSignConfiguration);
                }
            }
        }

        let pd_vertices = &self.power_diagram.vertices;
        for vi in 0..nvid {
            let node1_id = self.power_diagram.points[iatom].my_vertices_ids[vi];
            if node1_id >= pd_vertices.len() {
                continue;
            }
            let node1 = &pd_vertices[node1_id];
            for kn in 0..4 {
                let node2_id = node1.end_point_ids[kn];
                if node2_id == crate::power_diagram::GeneratorRef::INVALID_ID || node2_id == node1_id {
                    continue;
                }
                if node2_id >= pd_vertices.len() {
                    continue;
                }
                if node2_id > node1_id {
                    continue;
                }
                let node2 = &pd_vertices[node2_id];
                if node1.power_value > Scalar::zero() || node2.power_value > Scalar::zero() {
                    continue;
                }

                let mut node2_contains_atom = false;
                for kg in 0..4 {
                    let g2ref = node2.generator_refs[kg];
                    if g2ref.is_valid()
                        && g2ref.kind == crate::power_diagram::GeneratorKind::Point
                        && g2ref.index == iatom
                    {
                        node2_contains_atom = true;
                        break;
                    }
                }
                if !node2_contains_atom {
                    continue;
                }

                let mut ptn = 0_usize;
                for kg in 0..4 {
                    let g1ref = node1.generator_refs[kg];
                    if !g1ref.is_valid() {
                        continue;
                    }
                    if g1ref.kind == crate::power_diagram::GeneratorKind::Point && g1ref.index == iatom {
                        continue;
                    }
                    let mut shared_generator = false;
                    for kh in 0..4 {
                        let g2ref = node2.generator_refs[kh];
                        if g1ref.is_valid() && g2ref.is_valid() && g1ref.kind == g2ref.kind && g1ref.index == g2ref.index {
                            shared_generator = true;
                            break;
                        }
                    }
                    if shared_generator {
                        if ptn >= 2 {
                            return Err(SasaError::SharedPartnerOverflow);
                        }
                        partner[ptn] = self.power_diagram.get_generator(g1ref).visited_as;
                        ptn += 1;
                    }
                }
                if ptn != 2 {
                    return Err(SasaError::InvalidZeroPointPartnerCount);
                }

                let ptn0 = partner[0] as usize;
                let ptn1 = partner[1] as usize;
                if ptn0 >= nnb || ptn1 >= nnb {
                    return Err(SasaError::PartnerOutOfRange);
                }
                self.nt[ptn0] += 1;
                self.nt[ptn1] += 1;

                if !self.fknot[ptn0] {
                    self.fknot[ptn0] = true;
                    self.knot[ptn0] = node1.position;
                } else {
                    self.volnb[ptn0] +=
                        (node1.position - self.knot[ptn0]).cross(&(node2.position - self.knot[ptn0])).dot(&self.e[ptn0]).abs();
                }
                if !self.fknot[ptn1] {
                    self.fknot[ptn1] = true;
                    self.knot[ptn1] = node1.position;
                } else {
                    self.volnb[ptn1] +=
                        (node1.position - self.knot[ptn1]).cross(&(node2.position - self.knot[ptn1])).dot(&self.e[ptn1]).abs();
                }
            }
        }

        for i in 0..nnb {
            if self.np[i] <= 0 {
                continue;
            }
            if self.np[i] % 2 != 0 {
                self.np[i] = -1;
                continue;
            }
            let np_i = self.np[i] as usize;
            if self.rang.len() < np_i || self.pos.len() < np_i {
                self.resize_pnt(np_i);
            }
            let e_i        = self.e[i];
            let sintheta_i = self.sintheta[i];
            let costheta_i = self.costheta[i];
            if Self::get_ang(
                self.np[i], &self.p[i][..np_i], e_i,
                sintheta_i, costheta_i,
                &mut self.ang[i][..np_i], &self.vx,
            ).is_err() {
                self.np[i] = -1;
                continue;
            }
            if Self::get_next(
                self.np[i],
                &mut self.ang[i][..np_i], &mut self.next[i][..np_i],
                &self.p[i][..np_i], e_i, &self.vx,
                &mut self.rang[..np_i], &mut self.pos[..np_i],
            ).is_err() {
                self.np[i] = -1;
                continue;
            }
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
                continue;
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
                    break;
                }

                let ip_next = self.next[ic2][ip2 as usize];
                let mut phi = self.ang[ic2][ip_next as usize] - self.ang[ic2][ip2 as usize];
                if phi < Scalar::zero() {
                    phi += two_pi;
                }

                if do_sasa {
                    let mut co = (self.e[ic1].dot(&self.e[ic2]) - self.costheta[ic1] * self.costheta[ic2])
                        / (self.sintheta[ic1] * self.sintheta[ic2]);
                    co = co.clamp(-Scalar::one(), Scalar::one());
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
                    self.dsasa_parts[ic1] += self.e[ic1] * ds1
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
                    local_dvol -= (vv + self.e[ic1] * (rad2 * scone)) * half;
                }

                if pt_idx == p_ini_idx {
                    break;
                }
            }

            if do_sasa && sasa_ia > four * Scalar::pi() {
                sasa_ia -= four * Scalar::pi();
            }
        }

        if nvx == 0 && nnb > 32 {
            for i in 0..nnb {
                if self.np[i] == 0 && self.costheta[i] > Scalar::from_f64(-0.75).unwrap() {
                    self.nt[i] = 1;
                }
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
                let nb_id_j = self.power_diagram.points[iatom].neighbours_ids[j];
                let nb_j = self.power_diagram.points[nb_id_j].position;
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
                    self.dsasa_parts[i] = self.e[i]
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
                    local_dvol -= self.e[i] * (half * rad2 * scone);
                }
            }
        }

        if self.with_sasa {
            local_sasa = rad2 * sasa_ia;
        }
        if self.with_dsasa {
            for i in 0..nnb {
                local_dsasa -= self.dsasa_parts[i];
            }
        }
        if self.with_vol {
            for i in 0..nnb {
                if self.fknot[i] {
                    vol3 += rad * self.volnb[i] * self.costheta[i];
                }
            }
            local_vol = rad * rad2 * sasa_ia / three + rad * rad2 * vol2 / six + vol3 / six;
        }

        for k in 0..nnb {
            let nb_id = self.power_diagram.points[iatom].neighbours_ids[k];
            if nb_id < self.power_diagram.points.len() {
                self.power_diagram.get_cell_mut(nb_id).visited_as = 0;
            }
        }
        if self.with_sasa  { self.sasa.push(local_sasa); }
        if self.with_vol   { self.vol.push(local_vol);   }
        if self.with_dsasa { self.dsasa.push(local_dsasa); }
        if self.with_dvol  { self.dvol.push(local_dvol);  }
        Ok(())
    }

    /// Runs per-atom SASA/volume evaluation for all atoms in the current power diagram.
    pub fn calc_all(&mut self) -> Result<(), SasaError> {
        let n = self.power_diagram.points.len();
        if self.with_sasa  { self.sasa.clear();  self.sasa.reserve(n);  }
        if self.with_vol   { self.vol.clear();   self.vol.reserve(n);   }
        if self.with_dsasa { self.dsasa.clear(); self.dsasa.reserve(n); }
        if self.with_dvol  { self.dvol.clear();  self.dvol.reserve(n);  }
        for i in 0..n {
            self.calc_single(i)?;
        }
        Ok(())
    }

    /// Recalculates the underlying power diagram for new coordinates/radii. Output arrays are
    /// repopulated by the next `calc_all()` call.
    pub fn update(
        &mut self,
        coords: impl Iterator<Item = Vector3<Scalar>>,
        radii: impl Iterator<Item = Scalar>,
        size: usize,
    ) {
        self.power_diagram.recalculate(coords, radii, size);
        self.update_tol_pow();
    }

    /// Returns per-atom SASA values from the most recent computation.
    pub fn per_atom_sasa(&self) -> &[Scalar] {
        assert!(self.with_sasa, "per_atom_sasa() called but with_sasa was not enabled");
        &self.sasa
    }

    /// Returns per-atom volume values from the most recent computation.
    pub fn per_atom_vol(&self) -> &[Scalar] {
        assert!(self.with_vol, "per_atom_vol() called but with_vol was not enabled");
        &self.vol
    }

    /// Returns per-atom volume gradients when derivative computation is enabled.
    pub fn per_atom_dvol(&self) -> &[Vector3<Scalar>] {
        assert!(self.with_dvol, "per_atom_dvol() called but with_dvol was not enabled");
        &self.dvol
    }

    /// Returns per-atom SASA gradients when derivative computation is enabled.
    pub fn per_atom_dsasa(&self) -> &[Vector3<Scalar>] {
        assert!(self.with_dsasa, "per_atom_dsasa() called but with_dsasa was not enabled");
        &self.dsasa
    }
}
