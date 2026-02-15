use nalgebra::Vector3;

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
pub struct ZeroPoint<Scalar> {
    pub pos: Scalar,
    pub from_id: usize,
    pub branch: i32,
    pub generator_refs: [GeneratorRef; 3],
}

#[derive(Debug, Clone)]
pub struct Vertex<Scalar> {
    pub rrv: Scalar,
    pub invalid: bool,
    pub generator_refs: [GeneratorRef; 4],
    pub position: Vector3<Scalar>,
    pub power_value: Scalar,
    pub end_point_ids: [usize; 4],
}

#[derive(Debug, Clone)]
pub struct Cell<Scalar> {
    pub visited_as: i32,
    pub position: Vector3<Scalar>,
    pub r: Scalar,
    pub r2: Scalar,
    pub bond_to_id: usize,
    pub neighbours_ids: Vec<usize>,
    pub my_vertices_ids: Vec<usize>,
    pub my_zero_points: Vec<i32>,
}

#[derive(Debug, Clone)]
pub struct PowerDiagram<Scalar> {
    pub center: Vector3<Scalar>,
    pub maxr2: Scalar,
    points: Vec<Cell<Scalar>>,
    vertices: Vec<Vertex<Scalar>>,
    zeros: Vec<ZeroPoint<Scalar>>,
}

impl<Scalar> PowerDiagram<Scalar>
where
    Scalar: nalgebra::RealField + Copy,
{
    pub fn create(
        size: usize,
        pos_begin: impl Iterator<Item = Vector3<Scalar>>,
        strength_begin: impl Iterator<Item = Scalar>,
        bond_to_begin: impl Iterator<Item = i32>,
    ) -> Self {
        let _ = bond_to_begin.count();
        let mut points = Vec::with_capacity(size);
        for (position, r) in pos_begin.zip(strength_begin).take(size) {
            let r2 = r * r;
            points.push(Cell {
                visited_as: 0,
                position,
                r,
                r2,
                bond_to_id: GeneratorRef::INVALID_ID,
                neighbours_ids: Vec::new(),
                my_vertices_ids: Vec::new(),
                my_zero_points: Vec::new(),
            });
        }

        let maxr2 = points
            .iter()
            .map(|p| p.r2)
            .fold(Scalar::zero(), |acc, v| if v > acc { v } else { acc });

        Self {
            center: Vector3::zeros(),
            maxr2,
            points,
            vertices: Vec::new(),
            zeros: Vec::new(),
        }
    }

    pub fn recalculate(
        &mut self,
        pos_it: impl Iterator<Item = Vector3<Scalar>>,
        strength_it: impl Iterator<Item = Scalar>,
        size: usize,
    ) {
        self.points.clear();
        for (position, r) in pos_it.zip(strength_it).take(size) {
            let r2 = r * r;
            self.points.push(Cell {
                visited_as: 0,
                position,
                r,
                r2,
                bond_to_id: GeneratorRef::INVALID_ID,
                neighbours_ids: Vec::new(),
                my_vertices_ids: Vec::new(),
                my_zero_points: Vec::new(),
            });
        }
        self.vertices.clear();
        self.zeros.clear();
    }

    pub fn add_more(
        &mut self,
        pos_it: impl Iterator<Item = Vector3<Scalar>>,
        strength_it: impl Iterator<Item = Scalar>,
        new_size: usize,
    ) {
        let start = self.points.len();
        for (position, r) in pos_it.zip(strength_it).take(new_size.saturating_sub(start)) {
            let r2 = r * r;
            self.points.push(Cell {
                visited_as: 0,
                position,
                r,
                r2,
                bond_to_id: GeneratorRef::INVALID_ID,
                neighbours_ids: Vec::new(),
                my_vertices_ids: Vec::new(),
                my_zero_points: Vec::new(),
            });
        }
    }

    pub fn revert(&mut self) {}

    pub fn fill_all_my_vertices(&mut self) {}

    pub fn fill_all_neighbours(&mut self) {}

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

    pub fn get_zeros(&self) -> &[ZeroPoint<Scalar>] {
        &self.zeros
    }
}
