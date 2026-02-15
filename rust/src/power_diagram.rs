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
    revert_corner_owners: [usize; 8],
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
        let cube_lowest = params.lowest_corner - center;
        let cube_highest = params.highest_corner - center;

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
            revert_corner_owners: [GeneratorRef::INVALID_ID; 8],
            corner_owners: [GeneratorRef::INVALID_ID; 8],
            power_err: Scalar::default_epsilon(),
            insertion_error_scale: Scalar::zero(),
            replaced_ids: Vec::new(),
            invalids: Vec::new(),
            involved_refs: Vec::new(),
            planes: vec![EdgeEnds::default(); 64 * 64],
        };

        this.build_cube(cube_lowest, cube_highest);
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

    fn clear_interna(&mut self) {
        self.replaced_ids.clear();
        self.involved_refs.clear();
    }

    fn clear_involved(&mut self) {
        self.involved_refs.clear();
    }

    fn push_involved(&mut self, r#ref: GeneratorRef) {
        self.involved_refs.push(r#ref);
    }

    fn sort_involved_by_ref(&mut self) {
        self.involved_refs.sort_by_key(|r| {
            let kind_rank = match r.kind {
                GeneratorKind::Point => 0usize,
                GeneratorKind::Side => 1usize,
            };
            (kind_rank, r.index)
        });
    }

    fn insert_first(&mut self) {
        self.clear_interna();
        if self.points.is_empty() || self.vertices.len() < 8 {
            return;
        }
        for i in 0..8 {
            self.vertices[i].power_value = self.points[0].power(self.vertices[i].position);
            self.vertices[i].generator_refs[0] = GeneratorRef::new(GeneratorKind::Point, 0);
            self.points[0].my_vertices_ids.push(i);
            self.corner_owners[i] = 0;
        }
    }

    fn build_cube(&mut self, lowest: Vector3<Scalar>, highest: Vector3<Scalar>) {
        self.n_vertices = 1 << 3;
        self.side_generators.clear();
        for _ in 0..6 {
            self.side_generators.push(Cell::new(Vector3::zeros(), Scalar::zero()));
        }
        self.vertices.clear();
        self.vertices.resize(self.n_vertices, Vertex::default());

        let mut lhc = lowest;
        self.vertices[0].position = lowest;
        self.vertices[0].invalid = false;
        self.vertices[0].rrv = Scalar::zero();
        for j in (0..3).rev() {
            self.vertices[0].generator_refs[j + 1] = GeneratorRef::new(GeneratorKind::Side, j);
        }

        for i in 1..8 {
            let mut j = 0usize;
            while lhc[j] == highest[j] {
                lhc[j] = lowest[j];
                j += 1;
                if j >= 3 {
                    break;
                }
            }
            if j < 3 {
                lhc[j] = highest[j];
            }
            self.vertices[i].position = lhc;
            self.vertices[i].invalid = false;
            self.vertices[i].rrv = Scalar::zero();
            for j in (0..3).rev() {
                let side_idx = if lhc[j] == lowest[j] { j } else { j + 3 };
                self.vertices[i].generator_refs[j + 1] = GeneratorRef::new(GeneratorKind::Side, side_idx);
            }
        }

        for i in 0..8 {
            for d in 0..3 {
                let j = if (i >> d) % 2 == 1 { i - (1 << d) } else { i + (1 << d) };
                self.vertices[i].end_point_ids[d + 1] = j;
            }
            self.vertices[i].end_point_ids[0] = GeneratorRef::INVALID_ID;
        }
    }

    fn build_vertices(&mut self, _n: usize) {
        if self.points.is_empty() {
            self.n_vertices = 0;
            return;
        }
        self.maxr2 = self.points[0].r2;
        for p in &self.points {
            if p.r2 > self.maxr2 {
                self.maxr2 = p.r2;
            }
        }
        self.power_err = Scalar::from_f64(1000.0).unwrap() * Self::error(self.maxr2);
        self.insert_first();
        for i in 1..self.points.len() {
            let representative = self.prepare_insertion(i);
            if let Some(rep) = representative {
                let _ = self.do_insertion(i, rep);
            }
        }
        for c in 0..8.min(self.vertices.len()) {
            let corner_pos = self.vertices[c].position;
            let mut owner = 0usize;
            let mut best_power = self.points[0].power(corner_pos);
            for i in 1..self.points.len() {
                let p = self.points[i].power(corner_pos);
                if p < best_power {
                    best_power = p;
                    owner = i;
                }
            }
            self.corner_owners[c] = owner;
            self.vertices[c].generator_refs[0] = GeneratorRef::new(GeneratorKind::Point, owner);
            self.vertices[c].power_value = best_power;
        }
        self.build_edge_vertices_from_corners();
        // Incremental insertion for points[1..] pending full port.
    }

    fn build_edge_vertices_from_corners(&mut self) {
        let mut edge_vertices: Vec<Vertex<Scalar>> = Vec::new();
        for c in 0..8 {
            for d in 0..3 {
                let c2 = if (c >> d) % 2 == 1 { c - (1 << d) } else { c + (1 << d) };
                if c2 <= c || c2 >= 8 {
                    continue;
                }
                let owner_a = self.corner_owners[c];
                let owner_b = self.corner_owners[c2];
                if owner_a == GeneratorRef::INVALID_ID
                    || owner_b == GeneratorRef::INVALID_ID
                    || owner_a >= self.points.len()
                    || owner_b >= self.points.len()
                    || owner_a == owner_b
                {
                    continue;
                }

                let a = self.vertices[c].position;
                let b = self.vertices[c2].position;
                let dvec = b - a;
                let pa = self.points[owner_a].position;
                let pb = self.points[owner_b].position;
                let n = pb - pa;
                let rhs = Scalar::from_f64(0.5).unwrap()
                    * (pb.norm_squared() - pa.norm_squared() + self.points[owner_a].r2 - self.points[owner_b].r2);
                let denom = dvec.dot(&n);
                let mut t = if denom.abs() > self.power_err_scaled_epsilon() {
                    (rhs - a.dot(&n)) / denom
                } else {
                    Scalar::from_f64(0.5).unwrap()
                };
                if t < Scalar::zero() {
                    t = Scalar::zero();
                }
                if t > Scalar::one() {
                    t = Scalar::one();
                }
                let pos = a + dvec * t;

                let mut shared_sides: Vec<GeneratorRef> = Vec::new();
                for ga in 1..=3 {
                    let ra = self.vertices[c].generator_refs[ga];
                    for gb in 1..=3 {
                        let rb = self.vertices[c2].generator_refs[gb];
                        if ra.kind == GeneratorKind::Side && rb.kind == GeneratorKind::Side && ra.index == rb.index {
                            if !shared_sides.iter().any(|s| s.index == ra.index) {
                                shared_sides.push(ra);
                            }
                        }
                    }
                }
                if shared_sides.len() < 2 {
                    continue;
                }

                let mut v = Vertex::default();
                v.invalid = false;
                v.position = pos;
                v.generator_refs[0] = GeneratorRef::new(GeneratorKind::Point, owner_a);
                v.generator_refs[1] = GeneratorRef::new(GeneratorKind::Point, owner_b);
                v.generator_refs[2] = shared_sides[0];
                v.generator_refs[3] = shared_sides[1];
                v.end_point_ids[0] = c;
                v.end_point_ids[1] = c2;
                v.end_point_ids[2] = GeneratorRef::INVALID_ID;
                v.end_point_ids[3] = c2;
                v.power_value = self.points[owner_a].power(pos);
                edge_vertices.push(v);
            }
        }

        for v in edge_vertices {
            self.vertices.push(v);
        }
        self.n_vertices = self.vertices.len();
    }

    fn prepare_insertion(&self, point_id: usize) -> Option<usize> {
        if point_id >= self.points.len() {
            return None;
        }
        let hint = self.points[point_id].bond_to_id;
        if hint != GeneratorRef::INVALID_ID && hint < point_id && hint < self.points.len() {
            return Some(hint);
        }
        if let Some(id) = self.find_cell_inside_cube(self.points[point_id].position, None) {
            return Some(id);
        }
        self.find_replaced_vertex_corner(point_id)
            .and_then(|corner_id| self.corner_owners.get(corner_id).copied())
            .filter(|&owner| owner < self.points.len())
    }

    fn do_insertion(&mut self, point_id: usize, representative: usize) -> bool {
        if point_id >= self.points.len() || representative >= self.points.len() || point_id == representative {
            return false;
        }
        self.clear_interna();
        self.push_involved(GeneratorRef::new(GeneratorKind::Point, point_id));
        self.push_involved(GeneratorRef::new(GeneratorKind::Point, representative));
        self.points[point_id].bond_to_id = representative;

        if !self.points[point_id].neighbours_ids.contains(&representative) {
            self.points[point_id].neighbours_ids.push(representative);
        }
        if !self.points[representative].neighbours_ids.contains(&point_id) {
            self.points[representative].neighbours_ids.push(point_id);
        }
        true
    }

    fn find_replaced_vertex_corner(&self, point_id: usize) -> Option<usize> {
        if point_id >= self.points.len() || self.vertices.is_empty() {
            return None;
        }
        let insertion = &self.points[point_id];
        let mut best_corner = None;
        let mut best_val = Scalar::zero();
        for c in 0..8.min(self.vertices.len()) {
            let v = &self.vertices[c];
            if v.invalid {
                continue;
            }
            let val = insertion.power(v.position);
            if best_corner.is_none() || val < best_val {
                best_corner = Some(c);
                best_val = val;
            }
        }
        best_corner
    }

    fn has_virtual_generators(&self, vtx: &Vertex<Scalar>) -> bool {
        let refg = vtx.generator_refs[3];
        !(refg.kind == GeneratorKind::Point && refg.index < self.points.len())
    }

    fn n_virtual_generators(&self, vtx: &Vertex<Scalar>) -> usize {
        if !self.has_virtual_generators(vtx) {
            return 0;
        }
        if vtx.end_point_ids[0] == GeneratorRef::INVALID_ID {
            3
        } else {
            2
        }
    }

    fn above_power_err(&self, value: Scalar) -> bool {
        value > self.power_err
    }

    fn within_power_err(&self, value: Scalar) -> bool {
        value.abs() <= self.power_err
    }

    fn below_power_err(&self, value: Scalar) -> bool {
        value < self.power_err
    }

    fn below_neg_power_err(&self, value: Scalar) -> bool {
        value < -self.power_err
    }

    fn power_err_scaled_epsilon(&self) -> Scalar {
        self.power_err * Scalar::default_epsilon()
    }

    pub fn n_vertices_count(&self) -> usize {
        self.n_vertices
    }

    fn push_zero_from_edge(&mut self, source_id: usize, branch: i32, sol: Scalar) {
        let source = &self.vertices[source_id];
        let g0 = source.generator_refs[nth(0, branch) as usize];
        let g1 = source.generator_refs[nth(1, branch) as usize];
        let g2 = source.generator_refs[nth(2, branch) as usize];
        self.zeros.push(ZeroPoint {
            pos: sol,
            from_id: source_id,
            branch,
            generator_refs: [g0, g1, g2],
        });
    }

    pub fn find_cell_inside_cube(&self, pos: Vector3<Scalar>, hint_id: Option<usize>) -> Option<usize> {
        if self.points.is_empty() {
            return None;
        }
        let mut current = hint_id.unwrap_or(self.points.len() / 2);
        if current >= self.points.len() {
            current = self.points.len() / 2;
        }

        loop {
            let mut improved = false;
            let current_power = self.points[current].power(pos);
            for &neighbour_id in &self.points[current].neighbours_ids {
                if neighbour_id >= self.points.len() {
                    continue;
                }
                let neighbour_power = self.points[neighbour_id].power(pos);
                if neighbour_power < current_power {
                    current = neighbour_id;
                    improved = true;
                    break;
                }
            }
            if !improved {
                return Some(current);
            }
        }
    }

    pub fn recalculate(
        &mut self,
        pos_it: impl Iterator<Item = Vector3<Scalar>>,
        strength_it: impl Iterator<Item = Scalar>,
        size: usize,
    ) {
        let new_pos: Vec<Vector3<Scalar>> = pos_it.take(size).collect();
        let new_str: Vec<Scalar> = strength_it.take(size).collect();

        let (lowest, highest) = get_bounding_box(&new_pos, &new_str);
        self.center = (lowest + highest) * Scalar::from_f64(0.5).unwrap();

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
        self.n_vertices = 0;

        self.build_cube(lowest - self.center, highest - self.center);
        self.build_vertices(self.points.len());

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
        self.revert_corner_owners = self.corner_owners;

        for (p, s) in pos_vec.into_iter().zip(str_vec.into_iter()) {
            if self.params.radii_given {
                self.points.push(Cell::new(p - self.center, s));
            } else {
                self.points.push(Cell::with_power(p - self.center, s.sqrt(), s));
            }
        }

        self.build_vertices(self.points.len());

        if self.params.fill_my_vertices {
            self.fill_all_my_vertices_from(self.n_revert_points, self.n_revert_vertices.max(8));
        }
        if self.params.fill_neighbours {
            self.fill_all_neighbours_of_involved();
        }
        if self.params.fill_zero_points {
            self.fill_all_zero_points_from(self.n_revert_vertices.max(8), self.n_revert_zeros);
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
        self.corner_owners = self.revert_corner_owners;
        for c in 0..8.min(self.vertices.len()) {
            let owner = self.corner_owners[c];
            if owner != GeneratorRef::INVALID_ID && owner < self.points.len() {
                self.vertices[c].generator_refs[0] = GeneratorRef::new(GeneratorKind::Point, owner);
                self.vertices[c].power_value = self.points[owner].power(self.vertices[c].position);
            }
        }
        if self.params.fill_my_vertices {
            self.fill_all_my_vertices();
        }
        if self.params.fill_neighbours {
            self.fill_all_neighbours();
        }
        if self.params.fill_zero_points {
            self.fill_all_zero_points_from(8, 0);
        }

        self.n_revert_vertices = 0;
        self.n_revert_zeros = 0;
        self.n_revert_points = 0;
    }

    pub fn fill_all_my_vertices(&mut self) {
        self.fill_all_my_vertices_from(0, 8);
    }

    pub fn fill_all_my_vertices_from(&mut self, from_point: usize, from_vertex: usize) {
        for p in &mut self.points {
            if from_point == 0 {
                p.my_vertices_ids.clear();
                p.my_zero_points.clear();
            }
        }
        if from_point > 0 {
            self.clear_involved();
            for p in &mut self.points[from_point..] {
                p.my_vertices_ids.clear();
                p.my_zero_points.clear();
            }
        }
        for vi in from_vertex..self.n_vertices.min(self.vertices.len()) {
            let vtx = &self.vertices[vi];
            if vtx.invalid {
                continue;
            }
            let refs = vtx.generator_refs;
            let is_corner = vtx.end_point_ids[0] == GeneratorRef::INVALID_ID;
            if !self.has_virtual_generators(vtx) {
                for &refg in &refs {
                    if refg.kind == GeneratorKind::Point && refg.index < self.points.len() {
                        self.points[refg.index].my_vertices_ids.push(vi);
                        if from_point > 0 && self.points[refg.index].visited_as == 0 {
                            self.points[refg.index].visited_as = -1;
                            self.push_involved(GeneratorRef::new(GeneratorKind::Point, refg.index));
                        }
                    }
                }
            } else if !is_corner {
                for &refg in &refs {
                    if refg.kind == GeneratorKind::Side {
                        break;
                    }
                    if refg.kind == GeneratorKind::Point && refg.index < self.points.len() {
                        self.points[refg.index].my_vertices_ids.push(vi);
                        if from_point > 0 && self.points[refg.index].visited_as == 0 {
                            self.points[refg.index].visited_as = -1;
                            self.push_involved(GeneratorRef::new(GeneratorKind::Point, refg.index));
                        }
                    }
                }
            } else if !self.params.without_check {
                // Mirrors C++ assert path: corners should not be propagated to SASA ownership sets.
                return;
            }
        }
    }

    pub fn fill_all_neighbours(&mut self) {
        let has_any_my_vertices = self.points.iter().any(|p| !p.my_vertices_ids.is_empty());
        for p in &mut self.points {
            p.neighbours_ids.clear();
            p.visited_as = -1;
        }

        if has_any_my_vertices {
            for i in 0..self.points.len() {
                for vid in self.points[i].my_vertices_ids.clone() {
                    if vid >= self.vertices.len() {
                        continue;
                    }
                    let v = &self.vertices[vid];
                    for g in 0..=3 {
                        let gref = v.generator_refs[g];
                        if !gref.is_valid() || gref.kind != GeneratorKind::Point {
                            continue;
                        }
                        if gref.index == i || gref.index >= self.points.len() {
                            continue;
                        }
                        if self.points[gref.index].visited_as == i as i32 {
                            continue;
                        }
                        self.points[gref.index].visited_as = i as i32;
                        self.points[i].neighbours_ids.push(gref.index);
                    }
                }
            }
        }
        for i in 0..self.points.len() {
            if !self.points[i].neighbours_ids.is_empty() {
                continue;
            }
            for j in 0..self.points.len() {
                if i == j {
                    continue;
                }
                let dist = (self.points[j].position - self.points[i].position).norm();
                if dist <= self.points[i].r + self.points[j].r {
                    self.points[i].neighbours_ids.push(j);
                }
            }
        }

        for p in &mut self.points {
            p.visited_as = 0;
        }
    }

    pub fn fill_all_neighbours_of_involved(&mut self) {
        self.sort_involved_by_ref();
        for i in 0..self.involved_refs.len() {
            let involved = self.involved_refs[i];
            if involved.kind != GeneratorKind::Point || involved.index >= self.points.len() {
                continue;
            }
            let current_id = involved.index;
            let mut k = 0usize;
            while k < self.points[current_id].neighbours_ids.len() {
                let nb = self.points[current_id].neighbours_ids[k];
                if nb >= self.points.len() || self.points[nb].visited_as == -1 {
                    self.points[current_id].neighbours_ids.remove(k);
                } else {
                    k += 1;
                }
            }
        }

        for i in 0..self.involved_refs.len() {
            let involved = self.involved_refs[i];
            if involved.kind != GeneratorKind::Point || involved.index >= self.points.len() {
                continue;
            }
            let current_id = involved.index;
            let vids = self.points[current_id].my_vertices_ids.clone();
            for vid in vids {
                if vid >= self.vertices.len() {
                    continue;
                }
                let v = &self.vertices[vid];
                if v.end_point_ids[0] == GeneratorRef::INVALID_ID {
                    continue;
                }
                for g in (0..=3).rev() {
                    let refg = v.generator_refs[g];
                    if refg.kind == GeneratorKind::Point
                        && refg.index < self.points.len()
                        && refg.index != current_id
                        && self.points[refg.index].visited_as <= current_id as i32
                    {
                        self.points[current_id].neighbours_ids.push(refg.index);
                        self.points[refg.index].visited_as = current_id as i32 + 1;
                    }
                }
            }
        }

        for i in 0..self.involved_refs.len() {
            let involved = self.involved_refs[i];
            if involved.kind == GeneratorKind::Point && involved.index < self.points.len() {
                self.points[involved.index].visited_as = 0;
            }
        }
    }

    pub fn clear_all_my_vertices(&mut self) {
        for p in &mut self.points {
            p.my_vertices_ids.clear();
            p.my_zero_points.clear();
        }
    }

    pub fn fill_all_zero_points(&mut self) {
        self.fill_all_zero_points_from(8, 0);
    }

    pub fn fill_all_zero_points_from(&mut self, from_vertex: usize, from_zero: usize) {
        if from_zero <= self.zeros.len() {
            self.zeros.truncate(from_zero);
        } else {
            self.zeros.clear();
        }
        for p in &mut self.points {
            p.my_zero_points.clear();
        }
        for vertex_index in from_vertex..self.n_vertices.min(self.vertices.len()) {
            let current = self.vertices[vertex_index].clone();
            if current.invalid {
                continue;
            }
            let boundary_ref = current.generator_refs[2];
            if !(boundary_ref.kind == GeneratorKind::Point && boundary_ref.index < self.points.len()) {
                continue;
            }
            let endpoint_start = if self.has_virtual_generators(&current) { 3 } else { 0 };
            for endpoint_idx in endpoint_start..=3 {
                let endpoint_id = current.end_point_ids[endpoint_idx];
                if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id <= vertex_index || endpoint_id >= self.n_vertices {
                    continue;
                }
                let endpoint = self.vertices[endpoint_id].clone();
                if current.power_value > Scalar::zero() {
                    let branch = endpoint_idx as i32;
                    let v3 = endpoint.power_value;
                    let v2 = current.power_value;
                    let refg = current.generator_refs[if branch == 0 { 1 } else { 0 }];
                    let v1 = self.get_generator(refg).power(current.position * Scalar::from_f64(2.0).unwrap() - endpoint.position);
                    let quot = Scalar::from_f64(2.0).unwrap() * (v1 + v3 - Scalar::from_f64(2.0).unwrap() * v2);
                    let rootsq = (v1 - v3) * (v1 - v3) - Scalar::from_f64(4.0).unwrap() * quot * v2;
                    if rootsq <= Scalar::zero() {
                        continue;
                    }
                    if self.below_power_err(quot) && v1 >= Scalar::zero() && v2 >= Scalar::zero() && v3 >= Scalar::zero() {
                        continue;
                    }
                    let rootquot = rootsq.sqrt() / quot;
                    let min = (v1 - v3) / quot;
                    let sol1 = min + rootquot;
                    let sol2 = min - rootquot;
                    if sol1 > Scalar::zero() && sol1 < Scalar::one() {
                        if endpoint.power_value > Scalar::zero() {
                            self.push_zero_from_edge(vertex_index, branch, sol1);
                            self.push_zero_from_edge(vertex_index, branch, sol2);
                        } else {
                            self.push_zero_from_edge(vertex_index, branch, sol1);
                        }
                    } else if sol2 > Scalar::zero() && sol2 < Scalar::one() {
                        self.push_zero_from_edge(vertex_index, branch, sol2);
                    } else {
                        self.push_zero_from_edge(vertex_index, branch, sol1);
                        self.push_zero_from_edge(vertex_index, branch, sol2);
                    }
                } else if endpoint.power_value > Scalar::zero() {
                    let branch = endpoint_idx as i32;
                    let v3 = endpoint.power_value;
                    let v2 = current.power_value;
                    let refg = current.generator_refs[if branch == 0 { 1 } else { 0 }];
                    let v1 = self.get_generator(refg).power(current.position * Scalar::from_f64(2.0).unwrap() - endpoint.position);
                    let quot = Scalar::from_f64(2.0).unwrap() * (v1 + v3 - Scalar::from_f64(2.0).unwrap() * v2);
                    let rootsq = (v1 - v3) * (v1 - v3) - Scalar::from_f64(4.0).unwrap() * quot * v2;
                    if rootsq <= Scalar::zero() {
                        continue;
                    }
                    if self.below_power_err(quot) && v1 >= Scalar::zero() && v2 >= Scalar::zero() && v3 >= Scalar::zero() {
                        continue;
                    }
                    let rootquot = rootsq.sqrt() / quot;
                    let min = (v1 - v3) / quot;
                    let sol1 = min + rootquot;
                    let sol2 = min - rootquot;
                    if sol1 > Scalar::zero() && sol1 < Scalar::one() {
                        self.push_zero_from_edge(vertex_index, branch, sol1);
                    } else {
                        self.push_zero_from_edge(vertex_index, branch, sol2);
                    }
                }
            }
        }

        for i in from_zero..self.zeros.len() {
            for g in 0..3 {
                let refg = self.zeros[i].generator_refs[g];
                if refg.is_valid() && refg.kind == GeneratorKind::Point && refg.index < self.points.len() {
                    self.points[refg.index].my_zero_points.push(i as i32);
                }
            }
        }
    }

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
        if zp.from_id == GeneratorRef::INVALID_ID || zp.from_id >= self.n_vertices || zp.from_id >= self.vertices.len() {
            return false;
        }
        if zp.branch < 0 || zp.branch > 3 {
            return false;
        }
        let from = &self.vertices[zp.from_id];
        if from.invalid {
            return false;
        }
        let to_id = from.end_point_ids[zp.branch as usize];
        if to_id == GeneratorRef::INVALID_ID || to_id >= self.n_vertices || to_id >= self.vertices.len() {
            return false;
        }
        !self.vertices[to_id].invalid
    }

    pub fn zero_point_pos(&self, zp: &ZeroPoint<Scalar>) -> Vector3<Scalar> {
        if !self.zero_point_valid(zp) {
            return Vector3::zeros();
        }
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
