use nalgebra::{RealField, Vector3};

#[derive(Debug, Clone)]
pub struct MyException;

#[derive(Debug, Clone)]
pub struct IdenticalPointException;

#[derive(Debug, Clone)]
pub struct VerticesFullException;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneratorKind {
    Point,
    Side,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneratorRef {
    pub kind: GeneratorKind,
    pub index: usize,
}

impl GeneratorRef {
    pub const INVALID_ID: usize = usize::MAX;

    pub const fn new(kind: GeneratorKind, index: usize) -> Self {
        Self { kind, index }
    }

    pub const fn invalid() -> Self {
        Self {
            kind: GeneratorKind::Point,
            index: Self::INVALID_ID,
        }
    }

    pub const fn is_valid(self) -> bool {
        self.index != Self::INVALID_ID
    }
}

#[derive(Debug, Clone)]
pub struct PowerDiagramRuntimeParams {
    pub radii_given: bool,
    pub fill_my_vertices: bool,
    pub fill_neighbours: bool,
    pub fill_zero_points: bool,
    pub with_warnings: bool,
    pub without_check: bool,
}

impl Default for PowerDiagramRuntimeParams {
    fn default() -> Self {
        Self {
            radii_given: true,
            fill_my_vertices: true,
            fill_neighbours: true,
            fill_zero_points: true,
            with_warnings: false,
            without_check: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PowerDiagramParams<Scalar>
where
    Scalar: RealField + Copy,
{
    pub create_vertices: bool,
    pub runpar: PowerDiagramRuntimeParams,
    pub lowest_corner: Vector3<Scalar>,
    pub highest_corner: Vector3<Scalar>,
    pub size: usize,
    pub pos: Vec<Vector3<Scalar>>,
    pub strength: Vec<Scalar>,
    pub bond_to: Vec<i32>,
}

impl<Scalar> PowerDiagramParams<Scalar>
where
    Scalar: RealField + Copy,
{
    pub fn with_radii_given(mut self, yes: bool) -> Self {
        self.runpar.radii_given = yes;
        self
    }

    pub fn with_calculate(mut self, yes: bool) -> Self {
        self.create_vertices = yes;
        self
    }

    pub fn with_my_vertices(mut self, yes: bool) -> Self {
        self.runpar.fill_my_vertices = yes;
        self
    }

    pub fn with_cells(mut self, yes: bool) -> Self {
        self.runpar.fill_neighbours = yes;
        self.runpar.fill_my_vertices = yes;
        self
    }

    pub fn with_zero_points(mut self, yes: bool) -> Self {
        self.runpar.fill_zero_points = yes;
        self.runpar.fill_my_vertices = yes;
        self
    }

    pub fn with_warnings(mut self, yes: bool) -> Self {
        self.runpar.with_warnings = yes;
        self
    }

    pub fn without_check(mut self, yes: bool) -> Self {
        self.runpar.without_check = yes;
        self
    }
}

#[derive(Debug, Clone)]
pub struct ZeroPoint<Scalar>
where
    Scalar: RealField + Copy,
{
    pub pos: Scalar,
    pub from_id: usize,
    pub branch: i32,
    pub generator_refs: [GeneratorRef; 3],
}

#[derive(Debug, Clone)]
pub struct Vertex<Scalar>
where
    Scalar: RealField + Copy,
{
    pub rrv: Scalar,
    pub invalid: bool,
    pub generator_refs: [GeneratorRef; 4],
    pub position: Vector3<Scalar>,
    pub power_value: Scalar,
    pub end_point_ids: [usize; 4],
}

impl<Scalar> Default for Vertex<Scalar>
where
    Scalar: RealField + Copy,
{
    fn default() -> Self {
        Self {
            rrv: Scalar::zero(),
            invalid: true,
            generator_refs: [GeneratorRef::invalid(); 4],
            position: Vector3::zeros(),
            power_value: Scalar::zero(),
            end_point_ids: [GeneratorRef::INVALID_ID; 4],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cell<Scalar>
where
    Scalar: RealField + Copy,
{
    pub visited_as: i32,
    pub position: Vector3<Scalar>,
    pub r: Scalar,
    pub r2: Scalar,
    pub bond_to_id: usize,
    pub neighbours_ids: Vec<usize>,
    pub my_vertices_ids: Vec<usize>,
    pub my_zero_points: Vec<i32>,
}

impl<Scalar> Cell<Scalar>
where
    Scalar: RealField + Copy,
{
    pub fn new(position: Vector3<Scalar>, root: Scalar) -> Self {
        Self {
            visited_as: 0,
            position,
            r: root,
            r2: root * root,
            bond_to_id: GeneratorRef::INVALID_ID,
            neighbours_ids: Vec::new(),
            my_vertices_ids: Vec::new(),
            my_zero_points: Vec::new(),
        }
    }

    pub fn with_power(position: Vector3<Scalar>, root: Scalar, power: Scalar) -> Self {
        Self {
            visited_as: 0,
            position,
            r: root,
            r2: power,
            bond_to_id: GeneratorRef::INVALID_ID,
            neighbours_ids: Vec::new(),
            my_vertices_ids: Vec::new(),
            my_zero_points: Vec::new(),
        }
    }

    pub fn power(&self, coord: Vector3<Scalar>) -> Scalar {
        (self.position - coord).norm_squared() - self.r2
    }
}

#[derive(Debug, Clone)]
pub struct EdgeEnds {
    pub a_id: usize,
    pub a_slot: i32,
}

impl Default for EdgeEnds {
    fn default() -> Self {
        Self {
            a_id: GeneratorRef::INVALID_ID,
            a_slot: -1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PowerDiagram<Scalar>
where
    Scalar: RealField + Copy,
{
    pub center: Vector3<Scalar>,
    pub maxr2: Scalar,

    params: PowerDiagramRuntimeParams,
    n_vertices: usize,
    n_unused: usize,

    unused: Vec<usize>,
    points: Vec<Cell<Scalar>>,
    vertices: Vec<Vertex<Scalar>>,
    zeros: Vec<ZeroPoint<Scalar>>,
    side_generators: Vec<Cell<Scalar>>,

    n_revert_vertices: usize,
    n_revert_zeros: usize,
    n_revert_points: usize,
    corner_owners: [usize; 8],

    power_err: Scalar,
    insertion_error_scale: Scalar,
    replaced_ids: Vec<usize>,
    invalids: Vec<usize>,
    involved_refs: Vec<GeneratorRef>,
    planes: Vec<EdgeEnds>,
}

impl<Scalar> PowerDiagram<Scalar>
where
    Scalar: RealField + Copy,
{
    pub const K_INVALID_ID: usize = GeneratorRef::INVALID_ID;

    pub fn create(
        size: usize,
        pos_begin: impl Iterator<Item = Vector3<Scalar>>,
        strength_begin: impl Iterator<Item = Scalar>,
        bond_to_begin: impl Iterator<Item = i32>,
    ) -> PowerDiagramParams<Scalar> {
        let pos: Vec<Vector3<Scalar>> = pos_begin.take(size).collect();
        let strength: Vec<Scalar> = strength_begin.take(size).collect();
        let bond_to: Vec<i32> = bond_to_begin.take(size).collect();

        let (lowest_corner, highest_corner) = get_bounding_box(&pos, &strength);

        PowerDiagramParams {
            create_vertices: true,
            runpar: PowerDiagramRuntimeParams::default(),
            lowest_corner,
            highest_corner,
            size,
            pos,
            strength,
            bond_to,
        }
    }

    pub fn from_params(params: PowerDiagramParams<Scalar>) -> Self {
        let center = (params.highest_corner + params.lowest_corner) * Scalar::from_f64(0.5).unwrap();

        let mut points = Vec::with_capacity(params.size);
        for i in 0..params.size {
            let pos = params.pos[i] - center;
            let s = params.strength[i];
            if params.runpar.radii_given {
                points.push(Cell::new(pos, s));
            } else {
                points.push(Cell::with_power(pos, s.sqrt(), s));
            }
            if i > 0 && i < params.bond_to.len() {
                let b = params.bond_to[i];
                points[i].bond_to_id = if b >= 0 { b as usize } else { GeneratorRef::INVALID_ID };
            }
        }

        let maxr2 = points
            .iter()
            .map(|p| p.r2)
            .fold(Scalar::zero(), |acc, v| if v > acc { v } else { acc });

        let mut this = Self {
            center,
            maxr2,
            params: params.runpar,
            n_vertices: 0,
            n_unused: 0,
            unused: Vec::new(),
            points,
            vertices: Vec::new(),
            zeros: Vec::new(),
            side_generators: Vec::new(),
            n_revert_vertices: 0,
            n_revert_zeros: 0,
            n_revert_points: 0,
            corner_owners: [GeneratorRef::INVALID_ID; 8],
            power_err: Scalar::default_epsilon(),
            insertion_error_scale: Scalar::zero(),
            replaced_ids: Vec::new(),
            invalids: Vec::new(),
            involved_refs: Vec::new(),
            planes: vec![EdgeEnds::default(); 64 * 64],
        };

        if params.create_vertices {
            this.build_vertices(this.points.len());
        }
        if this.params.fill_my_vertices {
            this.fill_all_my_vertices();
        }
        if this.params.fill_neighbours {
            this.fill_all_neighbours();
        }
        if this.params.fill_zero_points {
            this.fill_all_zero_points();
        }

        this
    }

    fn build_vertices(&mut self, _n: usize) {
        self.n_vertices = 0;
    }

    pub fn recalculate(
        &mut self,
        pos_it: impl Iterator<Item = Vector3<Scalar>>,
        strength_it: impl Iterator<Item = Scalar>,
        size: usize,
    ) {
        let new_pos: Vec<Vector3<Scalar>> = pos_it.take(size).collect();
        let new_str: Vec<Scalar> = strength_it.take(size).collect();

        self.points.clear();
        self.points.reserve(size);
        for i in 0..size {
            let p = new_pos[i] - self.center;
            let s = new_str[i];
            if self.params.radii_given {
                self.points.push(Cell::new(p, s));
            } else {
                self.points.push(Cell::with_power(p, s.sqrt(), s));
            }
        }

        self.vertices.clear();
        self.zeros.clear();
        self.replaced_ids.clear();
        self.invalids.clear();
        self.involved_refs.clear();

        if self.params.fill_my_vertices {
            self.fill_all_my_vertices();
        }
        if self.params.fill_neighbours {
            self.fill_all_neighbours();
        }
        if self.params.fill_zero_points {
            self.fill_all_zero_points();
        }
    }

    pub fn add_more(
        &mut self,
        pos_it: impl Iterator<Item = Vector3<Scalar>>,
        strength_it: impl Iterator<Item = Scalar>,
        new_size: usize,
    ) {
        let target_add = new_size.saturating_sub(self.points.len());
        let pos_vec: Vec<Vector3<Scalar>> = pos_it.take(target_add).collect();
        let str_vec: Vec<Scalar> = strength_it.take(target_add).collect();

        self.n_revert_vertices = self.n_vertices;
        self.n_revert_zeros = self.zeros.len();
        self.n_revert_points = self.points.len();

        for (p, s) in pos_vec.into_iter().zip(str_vec.into_iter()) {
            if self.params.radii_given {
                self.points.push(Cell::new(p - self.center, s));
            } else {
                self.points.push(Cell::with_power(p - self.center, s.sqrt(), s));
            }
        }

        if self.params.fill_my_vertices {
            self.fill_all_my_vertices();
        }
        if self.params.fill_neighbours {
            self.fill_all_neighbours();
        }
        if self.params.fill_zero_points {
            self.fill_all_zero_points();
        }
    }

    pub fn revert(&mut self) {
        if self.n_revert_points > 0 && self.n_revert_points <= self.points.len() {
            self.points.truncate(self.n_revert_points);
        }
        if self.n_revert_zeros <= self.zeros.len() {
            self.zeros.truncate(self.n_revert_zeros);
        }
        self.n_vertices = self.n_revert_vertices;

        self.n_revert_vertices = 0;
        self.n_revert_zeros = 0;
        self.n_revert_points = 0;
    }

    pub fn fill_all_my_vertices(&mut self) {}

    pub fn fill_all_neighbours(&mut self) {
        let n = self.points.len();
        for i in 0..self.points.len() {
            self.points[i].neighbours_ids.clear();
            self.points[i].neighbours_ids.reserve(n.saturating_sub(1));
            for j in 0..n {
                if i != j {
                    self.points[i].neighbours_ids.push(j);
                }
            }
        }
    }

    pub fn fill_all_zero_points(&mut self) {}

    pub fn get_points(&self) -> &[Cell<Scalar>] {
        &self.points
    }

    pub fn get_points_mut(&mut self) -> &mut [Cell<Scalar>] {
        &mut self.points
    }

    pub fn get_vertices(&self) -> &[Vertex<Scalar>] {
        &self.vertices
    }

    pub fn get_vertices_mut(&mut self) -> &mut [Vertex<Scalar>] {
        &mut self.vertices
    }

    pub fn get_zeros(&self) -> &[ZeroPoint<Scalar>] {
        &self.zeros
    }

    pub fn get_zero_points(&self) -> &[ZeroPoint<Scalar>] {
        &self.zeros
    }

    pub fn get_cell(&self, id: usize) -> &Cell<Scalar> {
        &self.points[id]
    }

    pub fn get_cell_mut(&mut self, id: usize) -> &mut Cell<Scalar> {
        &mut self.points[id]
    }

    pub fn get_generator(&self, aref: GeneratorRef) -> &Cell<Scalar> {
        match aref.kind {
            GeneratorKind::Point => &self.points[aref.index],
            GeneratorKind::Side => &self.side_generators[aref.index],
        }
    }

    pub fn get_points_compat(&self) -> &[Cell<Scalar>] {
        &self.points
    }

    pub fn get_points_cpp_name(&self) -> &[Cell<Scalar>] {
        &self.points
    }

    pub fn zero_point_valid(&self, zp: &ZeroPoint<Scalar>) -> bool {
        if zp.from_id == GeneratorRef::INVALID_ID || zp.from_id >= self.vertices.len() {
            return false;
        }
        if zp.branch < 0 || zp.branch > 3 {
            return false;
        }
        let from = &self.vertices[zp.from_id];
        let to_id = from.end_point_ids[zp.branch as usize];
        to_id != GeneratorRef::INVALID_ID && to_id < self.vertices.len()
    }

    pub fn zero_point_pos(&self, zp: &ZeroPoint<Scalar>) -> Vector3<Scalar> {
        let from = &self.vertices[zp.from_id];
        let to = &self.vertices[from.end_point_ids[zp.branch as usize]];
        to.position * zp.pos - from.position * (zp.pos - Scalar::one())
    }

    pub fn error(f: Scalar) -> Scalar {
        let min_over_eps = Scalar::min_value().unwrap() / Scalar::default_epsilon();
        if f > min_over_eps {
            f * Scalar::default_epsilon()
        } else if f < -min_over_eps {
            -f * Scalar::default_epsilon()
        } else {
            min_over_eps
        }
    }
}

pub fn nth(n: i32, without: i32) -> i32 {
    n + i32::from(without <= n)
}

pub fn get_bounding_box<Scalar>(
    pos: &[Vector3<Scalar>],
    strength: &[Scalar],
) -> (Vector3<Scalar>, Vector3<Scalar>)
where
    Scalar: RealField + Copy,
{
    if pos.is_empty() {
        return (Vector3::zeros(), Vector3::zeros());
    }

    let mut low = pos[0];
    let mut high = pos[0];

    for i in 0..pos.len().min(strength.len()) {
        for g in 0..3 {
            if pos[i][g] - strength[i] < low[g] {
                low[g] = pos[i][g] - strength[i];
            }
            if pos[i][g] + strength[i] > high[g] {
                high[g] = pos[i][g] + strength[i];
            }
        }
    }

    let center = (low + high) * Scalar::from_f64(0.5).unwrap();
    let additional_cube_size = Scalar::from_f64(2.0f64.powf(1.0 / 3.0) - 1.0).unwrap();
    low += (low - center) * additional_cube_size;
    high += (high - center) * additional_cube_size;

    (low, high)
}
