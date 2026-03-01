use nalgebra::{RealField, Vector3};

/// Reference kind used to distinguish real point generators from boundary side generators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GeneratorKind {
    Point,
    Side,
}

/// Compact handle to a generator entry used by vertices and traversal state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct GeneratorRef {
    /// Source table kind (`Point` for `points`, `Side` for `side_generators`).
    pub kind: GeneratorKind,
    /// Index within the source table selected by `kind`.
    pub index: usize,
}

impl GeneratorRef {
    pub const INVALID_ID: usize = usize::MAX;

    #[inline(always)]
    /// Constructs a valid generator reference with the given kind and table index.
    pub const fn new(kind: GeneratorKind, index: usize) -> Self {
        Self {
            kind,
            index,
        }
    }

    #[inline(always)]
    /// Returns a sentinel reference with `index == INVALID_ID`; used where no generator is assigned yet.
    pub const fn invalid() -> Self {
        Self {
            kind: GeneratorKind::Point,
            index: Self::INVALID_ID,
        }
    }

    #[inline(always)]
    /// Returns `true` if this reference points to an actual generator (index ≠ INVALID_ID).
    pub const fn is_valid(self) -> bool {
        self.index != Self::INVALID_ID
    }
}

/// Runtime feature switches controlling which derived topology caches are maintained.
#[derive(Debug, Clone)]
pub struct PowerDiagramRuntimeParams {
    /// Interpret input strengths as radii (`true`) or raw power weights (`false`).
    pub radii_given: bool,
    /// Populate `Cell::my_vertices_ids` ownership lists.
    pub fill_my_vertices: bool,
    /// Populate `Cell::neighbours_ids` adjacency lists.
    pub fill_neighbours: bool,
    /// Populate zero-crossing caches (`zeros` and `Cell::my_zero_points`).
    pub fill_zero_points: bool,
    /// Enable warning output for corrective insertion fallbacks.
    pub with_warnings: bool,
    /// Skip strict consistency checks in edge/corner cases.
    pub without_check: bool,
}

impl Default for PowerDiagramRuntimeParams {
    /// Executes the `default` step of the power-diagram algorithm/state machine.
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

/// Input bundle for creating/rebuilding a power diagram from weighted points.
#[derive(Debug, Clone)]
pub struct PowerDiagramParams<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Build vertices/topology immediately when `true`.
    pub create_vertices: bool,
    /// Runtime cache/behavior toggles.
    pub runpar: PowerDiagramRuntimeParams,
    /// Lower corner of initial clipping cube.
    pub lowest_corner: Vector3<Scalar>,
    /// Upper corner of initial clipping cube.
    pub highest_corner: Vector3<Scalar>,
    /// Number of input points to consume.
    pub size: usize,
    /// Input point coordinates.
    pub pos: Vec<Vector3<Scalar>>,
    /// Input strengths (radii or powers, depending on `runpar.radii_given`).
    pub strength: Vec<Scalar>,
    /// Parent-hint map for incremental insertion seed selection.
    pub bond_to: Vec<i32>,
}

impl<Scalar> PowerDiagramParams<Scalar>
where
    Scalar: RealField + Copy,
{
    #[inline(always)]
    /// Sets the `radii_given` flag and returns `self` for chaining.
    pub(crate) fn with_radii_given(mut self, yes: bool) -> Self {
        self.runpar.radii_given = yes;
        self
    }

    #[inline(always)]
    /// Sets the `create_vertices` flag and returns `self` for chaining.
    pub(crate) fn with_calculate(mut self, yes: bool) -> Self {
        self.create_vertices = yes;
        self
    }

    #[inline(always)]
    /// Sets the `fill_neighbours` flag and returns `self` for chaining.
    pub(crate) fn with_cells(mut self, yes: bool) -> Self {
        self.runpar.fill_neighbours = yes;
        self.runpar.fill_my_vertices = yes;
        self
    }

    #[inline(always)]
    /// Sets the `fill_zero_points` flag and returns `self` for chaining.
    pub(crate) fn with_zero_points(mut self, yes: bool) -> Self {
        self.runpar.fill_zero_points = yes;
        self.runpar.fill_my_vertices = yes;
        self
    }

    #[inline(always)]
    /// Sets the `with_warnings` flag and returns `self` for chaining.
    pub(crate) fn with_warnings(mut self, yes: bool) -> Self {
        self.runpar.with_warnings = yes;
        self
    }

    #[inline(always)]
    /// Sets the `without_check` flag and returns `self` for chaining.
    pub(crate) fn without_check(mut self, yes: bool) -> Self {
        self.runpar.without_check = yes;
        self
    }
}

/// Power-zero crossing sampled along a vertex edge, used for SASA contour extraction.
#[derive(Debug, Clone)]
pub(crate) struct ZeroPoint<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Interpolation parameter in `[0, 1]` along the source edge.
    pub pos: Scalar,
    /// Source vertex id owning the edge branch.
    pub from_id: usize,
    /// Branch index at `from_id` that defines the edge.
    pub branch: i32,
    /// Generator triple defining the local zero-crossing context.
    pub generator_refs: [GeneratorRef; 3],
}

/// One power-diagram vertex with generator tuple, adjacency links, and replacement flags.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Vertex<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Cached signed power difference used during insertion to classify this vertex as replaced or persisting.
    pub rrv: Scalar,
    /// Whether this vertex slot is currently inactive/disconnected.
    pub invalid: bool,
    /// Generator tuple defining the vertex.
    pub generator_refs: [GeneratorRef; 4],
    /// Vertex position in centered coordinates.
    pub position: Vector3<Scalar>,
    /// Cached power value against the current insertion front.
    pub power_value: Scalar,
    /// Adjacent vertex ids by branch slot (slot 0 is corner sentinel semantics).
    pub end_point_ids: [usize; 4],
}

impl<Scalar> Default for Vertex<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Executes the `default` step of the power-diagram algorithm/state machine.
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

impl<Scalar> Vertex<Scalar>
where
    Scalar: RealField + Copy,
{
    #[inline(always)]
    /// Executes the `is_corner` step of the power-diagram algorithm/state machine.
    pub(crate) fn is_corner(&self) -> bool {
        self.end_point_ids[0] == GeneratorRef::INVALID_ID
    }

    #[inline(always)]
    /// Executes the `is_connected` step of the power-diagram algorithm/state machine.
    pub(crate) fn is_connected(&self) -> bool {
        !self.invalid
    }

    #[inline(always)]
    /// Executes the `resolved_endpoint_id` step of the power-diagram algorithm/state machine.
    pub(crate) fn resolved_endpoint_id(&self, g: usize) -> usize {
        self.end_point_ids[g]
    }

    #[inline(always)]
    /// Executes the `resolved_generator_ref` step of the power-diagram algorithm/state machine.
    pub(crate) fn resolved_generator_ref(&self, g: usize) -> GeneratorRef {
        self.generator_refs[g]
    }

    /// Computes the signed power difference at this vertex position between generators `a` and `b`: positive means `a` has lower power here (new generator wins), negative means `b` wins.
    pub(crate) fn powerdiff3d(&self, a_cell: &Cell<Scalar>, b_cell: &Cell<Scalar>) -> Scalar {
        // power_a(v) - power_b(v) = (|v − a.center|² − a.r²) − (|v − b.center|² − b.r²)
        // Positive ⇒ generator a has lower power at v (a's cell claims this vertex).
        -b_cell.r2 + a_cell.r2 - (a_cell.position - b_cell.position).norm_squared()
            + Scalar::from_f64(2.0).unwrap() * (a_cell.position - b_cell.position).dot(&(self.position - b_cell.position))
    }

    /// Returns the endpoint slot index (0..4) whose value equals `comp_id`, searching from slot 3 downward. Returns 0 if not found.
    pub(crate) fn endpoint_slot_to(&self, comp_id: usize) -> usize {
        for g in (1..=3).rev() {
            if self.end_point_ids[g] == comp_id {
                return g;
            }
        }
        0
    }

    #[inline(always)]
    /// Interpolates the position on the edge from `self` to `persist` where the power difference is zero, using the stored RRV values.
    pub(crate) fn get_power_point_on_line2(&self, persist: &Self) -> Vector3<Scalar> {
        // RRV is linear along the edge by construction.
        // Zero of r(t) = (1-t)*self.rrv + t*persist.rrv is at t = self.rrv / (self.rrv − persist.rrv).
        (persist.position - self.position) * (self.rrv / (self.rrv - persist.rrv)) + self.position
    }

    /// Executes the `end_points_and_position_overwrite` step of the power-diagram algorithm/state machine.
    pub(crate) fn end_points_and_position_overwrite(&mut self, endpoint_id: usize, pos: Vector3<Scalar>) {
        self.end_point_ids[0] = endpoint_id;
        self.rrv = Scalar::zero();
        self.invalid = false;
        self.position = pos;
    }
}

/// One weighted generator cell and its local topology caches.
#[derive(Debug, Clone)]
pub struct Cell<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Temporary traversal mark for incremental update passes.
    pub visited_as: i32,
    /// Generator position in centered coordinates.
    pub position: Vector3<Scalar>,
    /// Radius (or sqrt(weight) in power-input mode).
    pub r: Scalar,
    /// Squared radius/power used by distance comparisons.
    pub r2: Scalar,
    /// Parent hint to a nearby already-inserted cell.
    pub bond_to_id: usize,
    /// Adjacent point-generator ids.
    pub neighbours_ids: Vec<usize>,
    /// Owned vertex ids of this cell.
    pub my_vertices_ids: Vec<usize>,
    /// Zero-crossing ids associated with this cell.
    pub my_zero_points: Vec<i32>,
}

impl<Scalar> Cell<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Executes the `new` step of the power-diagram algorithm/state machine.
    pub(crate) fn new(position: Vector3<Scalar>, root: Scalar) -> Self {
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

    /// Executes the `with_power` step of the power-diagram algorithm/state machine.
    pub(crate) fn with_power(position: Vector3<Scalar>, root: Scalar, power: Scalar) -> Self {
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

    #[inline(always)]
    /// Evaluates the weighted power of this sphere at `coord`: `|center − coord|² − r²`. Negative inside the sphere's power region, positive outside.
    pub(crate) fn power(&self, coord: Vector3<Scalar>) -> Scalar {
        (self.position - coord).norm_squared() - self.r2
    }
}

/// Temporary edge-pair holder used to connect newly created finite vertices.
#[derive(Debug, Clone)]
pub struct EdgeEnds {
    /// First endpoint id waiting for its pair.
    pub a_id: usize,
    /// Slot on `a_id` that must be linked when the pair arrives.
    pub a_slot: i32,
}

impl Default for EdgeEnds {
    /// Executes the `default` step of the power-diagram algorithm/state machine.
    fn default() -> Self {
        Self {
            a_id: GeneratorRef::INVALID_ID,
            a_slot: -1,
        }
    }
}

impl EdgeEnds {
    /// Pairs two half-edge endpoints and writes reciprocal connectivity links.
    fn store_or_connect<Scalar>(&mut self, owner: &mut PowerDiagram<Scalar>, pvertex_id: usize, slot: usize)
    where
        Scalar: RealField + Copy,
    {
        if self.a_id == GeneratorRef::INVALID_ID {
            self.a_id = pvertex_id;
            self.a_slot = slot as i32;
        } else {
            let other_id = self.a_id;
            let other_slot = self.a_slot;
            owner.set_vertex_endpoint_deferred(pvertex_id, slot, other_id);
            owner.set_vertex_endpoint_deferred(other_id, other_slot as usize, pvertex_id);
            self.a_id = GeneratorRef::INVALID_ID;
            self.a_slot = -1;
        }
    }
}

/// Core weighted power-diagram state and incremental insertion/revert engine.
#[derive(Debug, Clone)]
pub struct PowerDiagram<Scalar>
where
    Scalar: RealField + Copy,
{
    /// Global center shift applied to all coordinates.
    pub center: Vector3<Scalar>,
    /// Maximum squared radius/weight among current point generators.
    pub maxr2: Scalar,

    /// Runtime toggles for optional caches and validation behavior.
    params: PowerDiagramRuntimeParams,
    /// Number of active vertices in `vertices`.
    n_vertices: usize,
    /// Number of reusable ids currently reserved in `unused`.
    n_unused: usize,

    /// Free-list of detached vertex ids.
    unused: Vec<usize>,
    /// Real point generators.
    pub(crate) points: Vec<Cell<Scalar>>,
    /// Vertex pool for the current diagram.
    pub(crate) vertices: Vec<Vertex<Scalar>>,
    /// Cached zero-power crossings.
    pub(crate) zeros: Vec<ZeroPoint<Scalar>>,
    /// Virtual side generators that encode clipping cube planes.
    side_generators: Vec<Cell<Scalar>>,

    /// Revert snapshot: vertex count before `add_more`.
    n_revert_vertices: usize,
    /// Current point owner for each corner vertex.
    corner_owners: [usize; 8],

    /// Tolerance used for replacement/power-sign decisions.
    power_err: Scalar,
    /// Scale of the last insertion-failure correction.
    insertion_error_scale: Scalar,
    /// Vertex ids marked as replaced in current insertion.
    replaced_ids: Vec<usize>,
    /// Previously existing vertices invalidated by insertion.
    invalids: Vec<usize>,
    /// Generators participating in current local insertion neighborhood.
    involved_refs: Vec<GeneratorRef>,
    /// Pairing table for connecting newly created finite vertices.
    planes: Vec<EdgeEnds>,
}

/// Classification of a diagram vertex with respect to a new generator being inserted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReplaceState {
    /// The existing generators still dominate (negative power difference); the vertex is kept.
    Persisting,
    /// The new generator dominates at this vertex (positive power difference); the vertex will be removed.
    Replaced,
    /// The power difference is within numerical tolerance; treated as persisting.
    Ambiguous,
}

impl<Scalar> PowerDiagram<Scalar>
where
    Scalar: RealField + Copy,
{
    pub const K_INVALID_ID: usize = GeneratorRef::INVALID_ID;

    /// Collects atom positions and radii from iterators, computes a centered bounding box, and returns a params bundle ready for `from_params`.
    pub(crate) fn create(
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

    /// Builds the initial power diagram: centers coordinates, optionally converts radii to power values, constructs the bounding cube with 8 corner vertices, inserts all atom generators incrementally, then populates requested topology caches.
    pub(crate) fn from_params(params: PowerDiagramParams<Scalar>) -> Self {
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
            this.build_vertices(this.points.len(), 0);
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

    /// Clears the per-insertion scratch lists `replaced_ids` and `involved_refs` at the start and end of each insertion.
    fn clear_interna(&mut self) {
        self.replaced_ids.clear();
        self.involved_refs.clear();
    }

    /// Sets generator slot `slot` of vertex `id` to the given reference.
    fn set_vertex_generator(&mut self, vertex_id: usize, slot: usize, r#ref: GeneratorRef) {
        if vertex_id < self.vertices.len() && slot <= 3 {
            self.vertices[vertex_id].generator_refs[slot] = r#ref;
        }
    }

    /// Sets endpoint slot `slot` of vertex `id` to the given vertex index.
    fn set_vertex_endpoint_deferred(&mut self, vertex_id: usize, slot: usize, endpoint_id: usize) {
        if vertex_id < self.vertices.len() && slot <= 3 {
            self.vertices[vertex_id].end_point_ids[slot] = endpoint_id;
        }
    }

    /// Swaps both the generator reference and the endpoint id at two slots of a vertex, keeping both arrays consistent.
    fn swap_vertex_link_slots(&mut self, vertex_id: usize, a: usize, b: usize) {
        if vertex_id >= self.vertices.len() || a > 3 || b > 3 {
            return;
        }
        self.vertices[vertex_id].generator_refs.swap(a, b);
        self.vertices[vertex_id].end_point_ids.swap(a, b);
    }

    #[inline(always)]
    /// Clears the involved-generator list without deallocating its capacity.
    fn clear_involved(&mut self) {
        self.involved_refs.clear();
    }

    #[inline(always)]
    /// Appends a generator reference to the involved list.
    fn push_involved(&mut self, r#ref: GeneratorRef) {
        self.involved_refs.push(r#ref);
    }

    /// Handles the first atom specially: assigns it as the generator of all 8 cube corners and caches the initial power values at those vertices.
    fn insert_first(&mut self) {
        self.clear_interna();
        if self.points.is_empty() || self.vertices.len() < 8 {
            return;
        }
        self.points[0].my_vertices_ids.clear();
        for i in 0..8 {
            self.vertices[i].power_value = self.points[0].power(self.vertices[i].position);
            self.set_vertex_generator(i, 0, GeneratorRef::new(GeneratorKind::Point, 0));
            self.corner_owners[i] = 0;
            self.points[0].my_vertices_ids.push(i);
        }
    }

    /// Creates the 8 corner vertices of the bounding cube and the 6 planar side generators. Hardwires all adjacency slots so the cube forms a valid starting topology for incremental insertion.
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
                self.set_vertex_endpoint_deferred(i, d + 1, j);
            }
            self.set_vertex_endpoint_deferred(i, 0, GeneratorRef::INVALID_ID);
        }
        for i in 0..8 {
            for g in (1..=2).rev() {
                for j in (1..=g).rev() {
                    self.swap_vertex_link_slots(i, j, j + 1);
                }
            }
        }
    }

    /// Inserts atom generators one by one via `insert_first` (atom 0) and `insert_point` (atoms 1..n). On numerical failure retries with a scaled error tolerance; populates requested topology caches afterward.
    fn build_vertices(&mut self, n_points: usize, from: usize) {
        if self.points.is_empty() {
            self.n_vertices = 0;
            return;
        }
        self.maxr2 = self.points[0].r2;
        for i in from..n_points.min(self.points.len()) {
            if self.points[i].r2 > self.maxr2 {
                self.maxr2 = self.points[i].r2;
            }
        }
        self.power_err = Scalar::from_f64(1000.0).unwrap() * Self::error(self.maxr2);
        if from == 0 {
            self.insert_first();
        }
        let start = if from == 0 { 1 } else { from };
        for i in start..n_points.min(self.points.len()) {
            let mut done: u32 = 1;
            loop {
                let success = match self.prepare_insertion(i) {
                    Some(hint_id) => self.do_insertion(i, hint_id),
                    None => false,
                };
                if success {
                    break;
                }

                let error_scale = self.insertion_error_scale;

                let mut identical_point_id = GeneratorRef::INVALID_ID;
                done += 1;
                if done > 100 {
                    return;
                }

                {
                    let insertion_pos = self.points[i].position;
                    let mut closest_id = GeneratorRef::INVALID_ID;
                    let mut mindist = Scalar::zero();
                    for involved_idx in 1..self.involved_refs.len() {
                        let candidate_id = self.involved_id_at(involved_idx);
                        if candidate_id == GeneratorRef::INVALID_ID {
                            continue;
                        }
                        let dist = (self.points[candidate_id].position - insertion_pos).norm_squared();
                        if closest_id == GeneratorRef::INVALID_ID || dist < mindist {
                            mindist = dist;
                            closest_id = candidate_id;
                        }
                    }
                    if closest_id != GeneratorRef::INVALID_ID && Self::error(self.points[i].r) > mindist.sqrt() {
                        identical_point_id = closest_id;
                    }
                }

                if self.n_unused == 0 {
                    let my_vertices_len = self.points[i].my_vertices_ids.len();
                    let unused_len = self.unused.len();
                    let delta = (my_vertices_len as isize) - (unused_len as isize) + (self.n_unused as isize);
                    if delta >= 0 {
                        self.n_vertices = self.n_vertices.saturating_sub(delta as usize);
                    } else {
                        self.n_vertices = self.n_vertices.saturating_add((-delta) as usize);
                    }
                    for idx in 0..self.points[i].my_vertices_ids.len() {
                        let involved_vid = self.points[i].my_vertices_ids[idx];
                        if involved_vid == GeneratorRef::INVALID_ID || involved_vid >= self.vertices.len() {
                            continue;
                        }
                        if involved_vid < 8 {
                            self.n_vertices += 1;
                        } else {
                            self.vertices[involved_vid].invalid = true;
                        }
                    }
                }

                for unused_idx in self.n_unused..self.unused.len() {
                    let unused_id = self.unused[unused_idx];
                    if unused_id == GeneratorRef::INVALID_ID || unused_id >= self.vertices.len() {
                        continue;
                    }
                    self.vertices[unused_id].invalid = true;
                }
                self.n_unused = self.unused.len();

                for ridx in 0..self.replaced_ids.len() {
                    let replaced_id = self.replaced_ids[ridx];
                    if replaced_id == GeneratorRef::INVALID_ID || replaced_id >= self.vertices.len() {
                        continue;
                    }
                    let replaced_start = if self.vertices[replaced_id].is_corner() { 1 } else { 0 };
                    for endpoint_idx in replaced_start..=3 {
                        let endpoint_id = self.vertices[replaced_id].resolved_endpoint_id(endpoint_idx);
                        if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
                            continue;
                        }
                        if self.vertices[endpoint_id].rrv <= Scalar::zero() {
                            self.vertices[endpoint_id].rrv = Scalar::zero();
                            let endpoint_start = if self.vertices[endpoint_id].is_corner() { 1 } else { 0 };
                            for g1 in replaced_start..=3 {
                                for g2 in endpoint_start..=3 {
                                    let a0 = self.vertices[replaced_id].resolved_generator_ref(nth(0, g1 as i32) as usize);
                                    let a1 = self.vertices[replaced_id].resolved_generator_ref(nth(1, g1 as i32) as usize);
                                    let a2 = self.vertices[replaced_id].resolved_generator_ref(nth(2, g1 as i32) as usize);
                                    let b0 = self.vertices[endpoint_id].resolved_generator_ref(nth(0, g2 as i32) as usize);
                                    let b1 = self.vertices[endpoint_id].resolved_generator_ref(nth(1, g2 as i32) as usize);
                                    let b2 = self.vertices[endpoint_id].resolved_generator_ref(nth(2, g2 as i32) as usize);
                                    if a0 == b0 && a1 == b1 && a2 == b2 {
                                        self.set_vertex_endpoint_deferred(replaced_id, g1, endpoint_id);
                                        self.set_vertex_endpoint_deferred(endpoint_id, g2, replaced_id);
                                    }
                                }
                            }
                        }
                    }
                }

                let fallback_replaced_id = self.replaced_ids.first().copied().unwrap_or(GeneratorRef::INVALID_ID);
                for involved_idx in 1..self.involved_refs.len() {
                    let involved_id = self.involved_id_at(involved_idx);
                    if involved_id == GeneratorRef::INVALID_ID || self.points[involved_id].my_vertices_ids.is_empty() {
                        continue;
                    }
                    let representative_id = self.points[involved_id].my_vertices_ids[0];
                    if representative_id == GeneratorRef::INVALID_ID || representative_id >= self.vertices.len() {
                        continue;
                    }
                    if !self.vertices[representative_id].is_connected()
                        && fallback_replaced_id != GeneratorRef::INVALID_ID
                        && fallback_replaced_id < self.vertices.len()
                    {
                        self.points[involved_id].my_vertices_ids[0] = fallback_replaced_id;
                    }
                }

                self.set_involved_persisting_visited_to_zero();
                self.points[i].my_vertices_ids.clear();
                for ridx in 0..self.replaced_ids.len() {
                    let replaced_id = self.replaced_ids[ridx];
                    if replaced_id == GeneratorRef::INVALID_ID || replaced_id >= self.vertices.len() {
                        continue;
                    }
                    self.vertices[replaced_id].rrv = Scalar::zero();
                    self.vertices[replaced_id].invalid = false;
                }
                self.replaced_ids.clear();

                if identical_point_id != GeneratorRef::INVALID_ID {
                    if !self.points[identical_point_id].my_vertices_ids.is_empty() {
                        let identical_vid = self.points[identical_point_id].my_vertices_ids[0];
                        if identical_vid != GeneratorRef::INVALID_ID && identical_vid < self.vertices.len() {
                            self.points[i].my_vertices_ids.push(identical_vid);
                        }
                    }
                    break;
                } else {
                    self.points[i].r2 -= Scalar::from_f64(2.0f64.powi(done as i32)).unwrap() * error_scale;
                    if self.points[i].r2 >= Scalar::zero() {
                        self.points[i].r = self.points[i].r2.sqrt();
                    } else {
                        self.points[i].r = -(-self.points[i].r2).sqrt();
                    }
                }
            }
        }
        self.compact_unused_vertices();
    }

    /// Updates all endpoint back-pointers in adjacent vertices to reflect that vertex `from_id` has moved to `to_id`.
    fn move_vertex_network_update_only(&mut self, from_id: usize, to_id: usize) {
        if from_id >= self.vertices.len() || to_id >= self.vertices.len() {
            return;
        }
        let moved = self.vertices[from_id];
        let endpoint_dim = moved.resolved_endpoint_id(3);
        if endpoint_dim != GeneratorRef::INVALID_ID && endpoint_dim < self.vertices.len() {
            let slot = self.vertices[endpoint_dim].endpoint_slot_to(from_id);
            self.set_vertex_endpoint_deferred(endpoint_dim, slot, to_id);
        }
        for g in 0..3 {
            let endpoint_id = moved.resolved_endpoint_id(g);
            if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
                continue;
            }
            let slot = self.vertices[endpoint_id].endpoint_slot_to(from_id);
            self.set_vertex_endpoint_deferred(endpoint_id, slot, to_id);
        }
        self.vertices[to_id] = moved;
    }

    /// Compacts the vertex array by swapping each unused slot with the last live vertex and updating all adjacency pointers, then decrements `n_vertices`.
    fn compact_unused_vertices(&mut self) {
        self.unused.sort_unstable();
        let mut idx = self.unused.len();
        while idx > 0 {
            idx -= 1;
            let unused_id = self.unused[idx];
            if self.n_vertices == 0 {
                break;
            }
            if unused_id == GeneratorRef::INVALID_ID || unused_id >= self.vertices.len() {
                continue;
            }
            if unused_id != self.n_vertices - 1 {
                let source_id = self.n_vertices - 1;
                self.n_vertices -= 1;
                self.move_vertex_network_update_only(source_id, unused_id);
            } else {
                self.n_vertices -= 1;
            }
        }
        self.unused.clear();
        self.n_unused = 0;
    }

    /// Traverses the vertex graph from the hint cell to locate the replacement frontier for the new generator: identifies all vertices whose power becomes positive (they will be replaced) and collects all generators adjacent to those vertices (the involved set).
    fn prepare_insertion(&mut self, point_id: usize) -> Option<usize> {
        if point_id >= self.points.len() {
            return None;
        }
        let mut hint_id = self.get_representative_vertex(self.points[point_id].bond_to_id);
        if hint_id >= self.vertices.len() {
            hint_id = 0;
        }
        if hint_id >= self.vertices.len() {
            return None;
        }

        let mut value = {
            let base_ref = self.vertices[hint_id].resolved_generator_ref(0);
            self.vertices[hint_id].powerdiff3d(self.get_generator(base_ref), &self.points[point_id])
        };
        self.find_replaced_vertex(&mut hint_id, &mut value, point_id);

        let mut done: u32 = 1;
        loop {
            if done != 1 {
                let base_ref = self.vertices[hint_id].resolved_generator_ref(0);
                value = self.vertices[hint_id].powerdiff3d(self.get_generator(base_ref), &self.points[point_id]);
                self.find_replaced_vertex(&mut hint_id, &mut value, point_id);
            }

            if self.fill_replaced_persisting_and_involved(point_id, hint_id) {
                break;
            }

            self.set_involved_persisting_visited_to_zero();
            self.points[point_id].my_vertices_ids.clear();
            for ridx in 0..self.replaced_ids.len() {
                let replaced_id = self.replaced_ids[ridx];
                if replaced_id >= self.vertices.len() {
                    continue;
                }
                self.vertices[replaced_id].rrv = Scalar::zero();
                let start = if self.vertices[replaced_id].is_corner() { 1 } else { 0 };
                for g in start..=3 {
                    let endpoint_id = self.vertices[replaced_id].resolved_endpoint_id(g);
                    if endpoint_id != GeneratorRef::INVALID_ID && endpoint_id < self.vertices.len() {
                        self.vertices[endpoint_id].rrv = Scalar::zero();
                    }
                }
            }
            self.replaced_ids.clear();

            let delta = Scalar::from_f64(2.0f64.powi(done as i32)).unwrap() * self.power_err;
            self.points[point_id].r2 -= delta;
            if self.points[point_id].r2 > Scalar::zero() {
                self.points[point_id].r = self.points[point_id].r2.sqrt();
            } else {
                self.points[point_id].r = -(-self.points[point_id].r2).sqrt();
            }
            done += 1;
            if done > 100 {
                return None;
            }
        }
        Some(hint_id)
    }

    /// Follows the `bond_to` hint chain to find an already-inserted cell with a valid vertex, returning that vertex id as a traversal seed.
    fn get_representative_vertex(&self, start_id: usize) -> usize {
        let mut current_id = start_id;
        while current_id != GeneratorRef::INVALID_ID && current_id < self.points.len() {
            let current = &self.points[current_id];
            let mut first_connected = GeneratorRef::INVALID_ID;
            for &vid in &current.my_vertices_ids {
                if vid == GeneratorRef::INVALID_ID || vid >= self.vertices.len() {
                    continue;
                }
                let candidate = &self.vertices[vid];
                if !candidate.is_connected() {
                    continue;
                }
                if first_connected == GeneratorRef::INVALID_ID {
                    first_connected = vid;
                }
                for g in 0..=3 {
                    let refg = candidate.resolved_generator_ref(g);
                    if refg.kind == GeneratorKind::Point && refg.index == current_id {
                        return vid;
                    }
                }
            }
            if first_connected != GeneratorRef::INVALID_ID {
                return first_connected;
            }
            if current.bond_to_id == current_id {
                break;
            }
            current_id = current.bond_to_id;
        }
        0
    }

    /// Orchestrates vertex creation for the new generator: creates frontier vertices, connects them pairwise along shared generator planes, recycles replaced slots into the unused free-list, and assigns representative vertices to affected cells.
    fn do_insertion(&mut self, point_id: usize, hint_id: usize) -> bool {
        if point_id >= self.points.len() || hint_id >= self.vertices.len() {
            return false;
        }
        if !self.create_finite_vertices_from_replaced() {
            return false;
        }
        if !self.connect_new_finites_among_themselves_3d() {
            return false;
        }
        self.update_unused();
        self.assign_representative_vertices_to_cells(hint_id);
        self.set_involved_persisting_visited_to_zero();
        true
    }

    /// Appends a generator to the involved set and stamps its cell with its position in `involved_refs` (used to look up the index by reference during vertex connection).
    fn add_to_involved(&mut self, r#ref: GeneratorRef) {
        if !self.valid_generator_ref(r#ref) {
            return;
        }
        match r#ref.kind {
            GeneratorKind::Point => self.points[r#ref.index].visited_as = self.involved_refs.len() as i32,
            GeneratorKind::Side => self.side_generators[r#ref.index].visited_as = self.involved_refs.len() as i32,
        }
        self.push_involved(r#ref);
    }

    /// Returns the insertion-time mark stored on the generator's cell, which encodes its index in `involved_refs`.
    fn generator_visited_as(&self, r#ref: GeneratorRef) -> i32 {
        if !self.valid_generator_ref(r#ref) {
            return 0;
        }
        match r#ref.kind {
            GeneratorKind::Point => self.points[r#ref.index].visited_as,
            GeneratorKind::Side => self.side_generators[r#ref.index].visited_as,
        }
    }

    /// Returns the point index of the involved generator at the given position, or `INVALID_ID` if it is a Side generator.
    fn involved_id_at(&self, index: usize) -> usize {
        if index >= self.involved_refs.len() {
            return GeneratorRef::INVALID_ID;
        }
        let r = self.involved_refs[index];
        if r.kind == GeneratorKind::Point && r.index < self.points.len() {
            r.index
        } else {
            GeneratorRef::INVALID_ID
        }
    }

    /// Computes the signed power difference (rrv) of this vertex with respect to the new generator vs. the generator it shares the boundary with. Returns `Replaced` if rrv > power_err, `Persisting` if rrv < -power_err, or `Ambiguous` if within tolerance.
    fn finite_replaced(&mut self, vertex_id: usize, cell_id: usize) -> ReplaceState {
        if vertex_id >= self.vertices.len() || cell_id == GeneratorRef::INVALID_ID || cell_id >= self.points.len() {
            return ReplaceState::Ambiguous;
        }
        let base_ref = self.vertices[vertex_id].resolved_generator_ref(0);
        if !self.valid_generator_ref(base_ref) {
            self.vertices[vertex_id].rrv = Scalar::zero();
            return ReplaceState::Ambiguous;
        }
        let rrv = self.vertices[vertex_id].powerdiff3d(&self.points[cell_id], self.get_generator(base_ref));
        self.vertices[vertex_id].rrv = rrv;
        if self.vertices[vertex_id].rrv > self.power_err {
            return ReplaceState::Replaced;
        }
        if self.vertices[vertex_id].rrv < -self.power_err {
            return ReplaceState::Persisting;
        }
        self.vertices[vertex_id].rrv = Scalar::zero();
        ReplaceState::Ambiguous
    }

    /// Dispatches to corner or finite replacement logic depending on whether the vertex has virtual (Side) generators.
    fn replace_check(&mut self, self_id: usize) -> bool {
        if self_id >= self.vertices.len() {
            return false;
        }
        if self.vertices[self_id].is_corner() {
            self.corner_replace_check(self_id)
        } else {
            self.finite_replace_check(self_id)
        }
    }

    /// Checks whether a finite vertex is replaced by the new generator and, if so, recurses to all its endpoint neighbors.
    fn finite_replace_check(&mut self, self_id: usize) -> bool {
        let involved_front_id = self.involved_id_at(0);
        if involved_front_id == GeneratorRef::INVALID_ID {
            return false;
        }
        match self.finite_replaced(self_id, involved_front_id) {
            ReplaceState::Ambiguous => false,
            ReplaceState::Replaced => self.finite_to_replaced_and_go(self_id),
            ReplaceState::Persisting => true,
        }
    }

    /// Same as `finite_replace_check` but for corner vertices, which skip endpoint slot 0.
    fn corner_replace_check(&mut self, self_id: usize) -> bool {
        let involved_front_id = self.involved_id_at(0);
        if involved_front_id == GeneratorRef::INVALID_ID {
            return false;
        }
        match self.finite_replaced(self_id, involved_front_id) {
            ReplaceState::Ambiguous => false,
            ReplaceState::Replaced => self.corner_to_replaced_and_go(self_id),
            ReplaceState::Persisting => true,
        }
    }

    /// Same as `finite_to_replaced_and_go` for corner vertices, iterating only endpoint slots 1..3.
    fn corner_to_replaced_and_go(&mut self, self_id: usize) -> bool {
        if self_id >= self.vertices.len() {
            return false;
        }
        self.replaced_ids.push(self_id);
        for g in 0..=3 {
            let refg = self.vertices[self_id].resolved_generator_ref(g);
            if !self.valid_generator_ref(refg) {
                continue;
            }
            if self.generator_visited_as(refg) == 0 {
                self.add_to_involved(refg);
            }
        }
        let involved_front_id = self.involved_id_at(0);
        if involved_front_id != GeneratorRef::INVALID_ID {
            self.points[involved_front_id].my_vertices_ids.push(self_id);
        }
        for g in (1..=3).rev() {
            let endpoint_id = self.vertices[self_id].resolved_endpoint_id(g);
            if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
                continue;
            }
            if self.vertices[endpoint_id].rrv == Scalar::zero() && !self.replace_check(endpoint_id) {
                return false;
            }
        }
        true
    }

    /// Marks a finite vertex as replaced, adds all its generators to the involved set, and recursively checks its endpoint neighbors for replacement.
    fn finite_to_replaced_and_go(&mut self, self_id: usize) -> bool {
        if self_id >= self.vertices.len() {
            return false;
        }
        self.replaced_ids.push(self_id);
        for g in 0..=3 {
            let refg = self.vertices[self_id].resolved_generator_ref(g);
            if !self.valid_generator_ref(refg) {
                continue;
            }
            if self.generator_visited_as(refg) == 0 {
                self.add_to_involved(refg);
            }
        }
        for g in 0..=3 {
            let endpoint_id = self.vertices[self_id].resolved_endpoint_id(g);
            if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
                continue;
            }
            if self.vertices[endpoint_id].rrv == Scalar::zero() && !self.replace_check(endpoint_id) {
                return false;
            }
        }
        true
    }

    /// Seeds the replacement search: adds the new generator to the involved set, evaluates its power difference at the hint vertex, and dispatches to the appropriate replacement traversal (replaced, persisting, or ambiguous).
    fn fill_replaced_persisting_and_involved(&mut self, this_id: usize, start_id: usize) -> bool {
        if this_id == GeneratorRef::INVALID_ID || this_id >= self.points.len() {
            return false;
        }
        self.clear_interna();
        if start_id == GeneratorRef::INVALID_ID || start_id >= self.vertices.len() {
            return false;
        }
        self.push_involved(GeneratorRef::new(GeneratorKind::Point, this_id));
        match self.finite_replaced(start_id, this_id) {
            ReplaceState::Ambiguous => false,
            ReplaceState::Persisting => true,
            ReplaceState::Replaced => {
                if self.vertices[start_id].is_corner() {
                    self.corner_to_replaced_and_go(start_id)
                } else {
                    self.finite_to_replaced_and_go(start_id)
                }
            }
        }
    }

    /// Walks outward from `start_id` to find the minimum-power vertex reachable, used as the seed for the replacement-frontier search. Tries single-hop neighbors first, then two-hop, then a global scan as fallback.
    fn find_replaced_vertex(&mut self, this_id: &mut usize, value: &mut Scalar, insertion_id: usize) {
        if *value < Scalar::zero() || *this_id >= self.vertices.len() {
            return;
        }
        let mut small_val = Scalar::max_value().unwrap();

        let this_vertex_id = *this_id;
        let this_start = if self.vertices[this_vertex_id].is_corner() { 1 } else { 0 };
        let mut idx = this_start;
        while idx <= 3 {
            let endpoint_id = self.vertices[this_vertex_id].resolved_endpoint_id(idx);
            if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
                idx += 1;
                continue;
            }
            let base_ref = self.vertices[endpoint_id].resolved_generator_ref(0);
            if !self.valid_generator_ref(base_ref) {
                idx += 1;
                continue;
            }
            let new_value = self.vertices[endpoint_id].powerdiff3d(self.get_generator(base_ref), &self.points[insertion_id]);
            if new_value < *value {
                *value = new_value;
                *this_id = endpoint_id;
                if *value < Scalar::zero() {
                    return;
                }
                idx = this_start;
                continue;
            } else if new_value == *value {
                small_val = new_value;
            }
            idx += 1;
        }
        if small_val != *value {
            return;
        }

        small_val = Scalar::max_value().unwrap();
        for g in this_start..=3 {
            let ep1_id = self.vertices[this_vertex_id].resolved_endpoint_id(g);
            if ep1_id == GeneratorRef::INVALID_ID || ep1_id >= self.vertices.len() {
                continue;
            }
            let ep1_start = if self.vertices[ep1_id].is_corner() { 1 } else { 0 };
            for g2 in ep1_start..=3 {
                let candidate_id = self.vertices[ep1_id].resolved_endpoint_id(g2);
                if candidate_id == GeneratorRef::INVALID_ID || candidate_id >= self.vertices.len() || candidate_id == *this_id {
                    continue;
                }
                let base_ref = self.vertices[candidate_id].resolved_generator_ref(0);
                if !self.valid_generator_ref(base_ref) {
                    continue;
                }
                let new_value = self.vertices[candidate_id].powerdiff3d(self.get_generator(base_ref), &self.points[insertion_id]);
                if new_value < *value {
                    *value = new_value;
                    *this_id = candidate_id;
                    return self.find_replaced_vertex(this_id, value, insertion_id);
                } else if new_value == *value {
                    small_val = new_value;
                }
            }
        }
        if small_val != *value {
            return;
        }

        small_val = Scalar::max_value().unwrap();
        for g in this_start..=3 {
            let ep1_id = self.vertices[this_vertex_id].resolved_endpoint_id(g);
            if ep1_id == GeneratorRef::INVALID_ID || ep1_id >= self.vertices.len() {
                continue;
            }
            let ep1_start = if self.vertices[ep1_id].is_corner() { 1 } else { 0 };
            for g2 in ep1_start..=3 {
                let ep2_id = self.vertices[ep1_id].resolved_endpoint_id(g2);
                if ep2_id == GeneratorRef::INVALID_ID || ep2_id >= self.vertices.len() || ep2_id == *this_id {
                    continue;
                }
                let ep2_start = if self.vertices[ep2_id].is_corner() { 1 } else { 0 };
                for g3 in ep2_start..=3 {
                    let candidate_id = self.vertices[ep2_id].resolved_endpoint_id(g3);
                    if candidate_id == GeneratorRef::INVALID_ID
                        || candidate_id >= self.vertices.len()
                        || candidate_id == ep1_id
                        || candidate_id == *this_id
                    {
                        continue;
                    }
                    let base_ref = self.vertices[candidate_id].resolved_generator_ref(0);
                    if !self.valid_generator_ref(base_ref) {
                        continue;
                    }
                    let new_value = self.vertices[candidate_id].powerdiff3d(self.get_generator(base_ref), &self.points[insertion_id]);
                    if new_value < *value {
                        *value = new_value;
                        *this_id = candidate_id;
                        return self.find_replaced_vertex(this_id, value, insertion_id);
                    } else if new_value == *value {
                        small_val = new_value;
                    }
                }
            }
        }
        if small_val != *value {
            return;
        }
        for vi in 0..self.n_vertices.min(self.vertices.len()) {
            if !self.vertices[vi].is_connected() {
                continue;
            }
            let base_ref = self.vertices[vi].resolved_generator_ref(0);
            if !self.valid_generator_ref(base_ref) {
                continue;
            }
            let v = self.vertices[vi].powerdiff3d(self.get_generator(base_ref), &self.points[insertion_id]);
            if v < *value {
                *value = v;
                *this_id = vi;
            }
        }
    }

    /// Resets the `visited_as` mark and the `rrv` cache on all vertices of persisting generators, cleaning up insertion-time state.
    fn set_involved_persisting_visited_to_zero(&mut self) {
        for i in 1..self.involved_refs.len() {
            let refg = self.involved_refs[i];
            if !self.valid_generator_ref(refg) {
                continue;
            }
            match refg.kind {
                GeneratorKind::Point => self.points[refg.index].visited_as = 0,
                GeneratorKind::Side => self.side_generators[refg.index].visited_as = 0,
            }
        }
        let involved_front_id = self.involved_id_at(0);
        if involved_front_id == GeneratorRef::INVALID_ID {
            return;
        }
        for vidx in 0..self.points[involved_front_id].my_vertices_ids.len() {
            let vid = self.points[involved_front_id].my_vertices_ids[vidx];
            if vid == GeneratorRef::INVALID_ID || vid >= self.vertices.len() {
                continue;
            }
            self.vertices[vid].rrv = Scalar::zero();
            if !self.vertices[vid].is_corner() {
                let endpoint_id = self.vertices[vid].resolved_endpoint_id(0);
                if endpoint_id != GeneratorRef::INVALID_ID && endpoint_id < self.vertices.len() {
                    self.vertices[endpoint_id].rrv = Scalar::zero();
                }
            } else {
                for g in 1..=3 {
                    let endpoint_id = self.vertices[vid].resolved_endpoint_id(g);
                    if endpoint_id != GeneratorRef::INVALID_ID && endpoint_id < self.vertices.len() {
                        self.vertices[endpoint_id].rrv = Scalar::zero();
                    }
                }
            }
        }
    }

    /// Interpolates a new vertex on the edge from a replaced vertex to its persisting neighbor using the RRV ratio, allocates or reuses a vertex slot, and initializes its generator references from the replaced vertex.
    fn try_to_build_vertex_on_edge(&mut self, this_id: usize, here: usize) -> bool {
        if this_id >= self.vertices.len() {
            self.insertion_error_scale = self.power_err;
            return false;
        }
        let persisting_id = self.vertices[this_id].resolved_endpoint_id(here);
        if persisting_id == GeneratorRef::INVALID_ID || persisting_id >= self.vertices.len() {
            self.insertion_error_scale = self.power_err;
            return false;
        }
        let new_pos = {
            let this_v = &self.vertices[this_id];
            let persisting = &self.vertices[persisting_id];
            this_v.get_power_point_on_line2(persisting)
        };

        let built_vertex_id = if self.n_unused == 0 {
            if self.n_vertices == self.vertices.len() {
                self.vertices.push(Vertex::default());
            }
            let id = self.n_vertices;
            self.n_vertices += 1;
            id
        } else {
            let id = self.unused[self.n_unused - 1];
            self.n_unused -= 1;
            id
        };
        if built_vertex_id >= self.vertices.len() {
            self.insertion_error_scale = self.power_err;
            return false;
        }
        self.vertices[built_vertex_id].end_points_and_position_overwrite(persisting_id, new_pos);
        self.init_new_vertex_from_replaced(this_id, here, built_vertex_id)
    }

    /// Copies generator references from a replaced vertex to a new vertex, substituting the new generator at slot 0 and skipping the dropped generator slot.
    fn init_new_vertex_from_replaced(&mut self, this_id: usize, keep: usize, self_id: usize) -> bool {
        let involved_front_id = self.involved_id_at(0);
        if involved_front_id == GeneratorRef::INVALID_ID || this_id >= self.vertices.len() || self_id >= self.vertices.len() {
            self.insertion_error_scale = self.power_err;
            return false;
        }
        let pwr = self.points[involved_front_id].power(self.vertices[self_id].position);
        self.vertices[self_id].power_value = pwr;
        self.vertices[self_id].invalid = false;

        for g in (1..=3).rev() {
            let src = g - usize::from(g <= keep);
            let r = self.vertices[this_id].resolved_generator_ref(src);
            self.set_vertex_generator(self_id, g, r);
        }
        self.set_vertex_generator(self_id, 0, GeneratorRef::new(GeneratorKind::Point, involved_front_id));
        self.points[involved_front_id].my_vertices_ids.push(self_id);

        let endpoint_id = self.vertices[self_id].resolved_endpoint_id(0);
        if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
            self.insertion_error_scale = self.power_err;
            return false;
        }
        let slot = self.vertices[endpoint_id].endpoint_slot_to(this_id);
        self.set_vertex_endpoint_deferred(endpoint_id, slot, self_id);
        if self.vertices[self_id].power_value.abs() <= self.power_err {
            self.insertion_error_scale = self.power_err;
            return false;
        }
        true
    }

    /// For each replaced vertex, walks its edges to persisting neighbors and interpolates a new vertex at the power-zero crossing on each such edge.
    fn create_finite_vertices_from_replaced(&mut self) -> bool {
        for ridx in 0..self.replaced_ids.len() {
            let replaced_id = self.replaced_ids[ridx];
            if replaced_id == GeneratorRef::INVALID_ID || replaced_id >= self.vertices.len() {
                continue;
            }
            let start = if self.vertices[replaced_id].is_corner() { 1 } else { 0 };
            for g in (start..=3).rev() {
                let endpoint_id = self.vertices[replaced_id].resolved_endpoint_id(g);
                if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id >= self.vertices.len() {
                    continue;
                }
                if self.vertices[endpoint_id].rrv <= Scalar::zero() && !self.try_to_build_vertex_on_edge(replaced_id, g) {
                    return false;
                }
            }
        }
        true
    }

    /// Pairs newly created frontier vertices that share the same two involved generators (lie on the same new facet) and links them as each other's endpoints.
    fn connect_new_finites_among_themselves_3d(&mut self) -> bool {
        let involved_count = self.involved_refs.len();
        if involved_count == 0 {
            return true;
        }
        let needed = involved_count * involved_count;
        if needed > self.planes.len() {
            self.planes.resize(needed, EdgeEnds::default());
        }

        let involved_front_id = self.involved_id_at(0);
        if involved_front_id == GeneratorRef::INVALID_ID {
            return false;
        }
        for vertex_idx in 0..8.min(self.vertices.len()) {
            if self.vertices[vertex_idx].rrv > Scalar::zero() {
                self.vertices[vertex_idx].power_value = self.points[involved_front_id].power(self.vertices[vertex_idx].position);
                self.set_vertex_generator(vertex_idx, 0, GeneratorRef::new(GeneratorKind::Point, involved_front_id));
            }
        }

        for vidx in 0..self.points[involved_front_id].my_vertices_ids.len() {
            let vid = self.points[involved_front_id].my_vertices_ids[vidx];
            if vid == GeneratorRef::INVALID_ID || vid >= self.vertices.len() {
                continue;
            }
            self.register_vertex_for_connection_3d(vid);
        }
        true
    }

    /// Records a new vertex in the edge-pair table for each of the three generator pairs it participates in, enabling `connect_new_finites_among_themselves_3d` to find and link matching vertices.
    fn register_vertex_for_connection_3d(&mut self, self_id: usize) {
        if self_id >= self.vertices.len() || self.involved_refs.is_empty() {
            return;
        }
        let g1_ref = self.vertices[self_id].resolved_generator_ref(1);
        let g2_ref = self.vertices[self_id].resolved_generator_ref(2);
        let g3_ref = self.vertices[self_id].resolved_generator_ref(3);
        let g1 = self.generator_visited_as(g1_ref);
        let g2 = self.generator_visited_as(g2_ref);
        let g3 = self.generator_visited_as(g3_ref);
        if g1 < 0 || g2 < 0 || g3 < 0 {
            return;
        }
        let n = self.involved_refs.len();
        let i1 = g1 as usize;
        let i2 = g2 as usize;
        let i3 = g3 as usize;
        if i1 >= n || i2 >= n || i3 >= n {
            return;
        }
        let idx_a = i2 * n + i1;
        let idx_b = i3 * n + i1;
        let idx_c = i3 * n + i2;
        if idx_a < self.planes.len() {
            let mut edge = std::mem::take(&mut self.planes[idx_a]);
            edge.store_or_connect(self, self_id, 3);
            self.planes[idx_a] = edge;
        }
        if idx_b < self.planes.len() {
            let mut edge = std::mem::take(&mut self.planes[idx_b]);
            edge.store_or_connect(self, self_id, 2);
            self.planes[idx_b] = edge;
        }
        if idx_c < self.planes.len() {
            let mut edge = std::mem::take(&mut self.planes[idx_c]);
            edge.store_or_connect(self, self_id, 1);
            self.planes[idx_c] = edge;
        }
    }

    /// Removes a vertex id from a cell's `my_vertices_ids` cache, used when a vertex is invalidated during insertion.
    fn erase_cell_my_vertex_by_id(&mut self, r#ref: GeneratorRef, id: usize) {
        if !self.valid_generator_ref(r#ref) {
            return;
        }
        let list = match r#ref.kind {
            GeneratorKind::Point => &mut self.points[r#ref.index].my_vertices_ids,
            GeneratorKind::Side => &mut self.side_generators[r#ref.index].my_vertices_ids,
        };
        if let Some(pos) = list.iter().position(|&x| x == id) {
            list.remove(pos);
        }
    }

    /// Transfers replaced vertices to the unused free-list for recycling. Vertices below the revert snapshot boundary go to `invalids` instead and cannot be recycled across a revert.
    fn update_unused(&mut self) {
        self.unused.truncate(self.n_unused);
        let mut idx = 0usize;
        while idx < self.replaced_ids.len() {
            let replaced_id = self.replaced_ids[idx];
            if replaced_id == GeneratorRef::INVALID_ID || replaced_id >= self.vertices.len() {
                self.replaced_ids.swap_remove(idx);
                continue;
            }
            if self.vertices[replaced_id].is_corner() {
                self.vertices[replaced_id].rrv = Scalar::zero();
                self.replaced_ids.swap_remove(idx);
                continue;
            }
            if replaced_id < self.n_revert_vertices {
                self.invalids.push(replaced_id);
            }
            self.vertices[replaced_id].invalid = true;
            idx += 1;
        }

        if self.n_revert_vertices == 0 {
            for &replaced_id in &self.replaced_ids {
                if replaced_id != GeneratorRef::INVALID_ID && replaced_id < self.vertices.len() {
                    self.unused.push(replaced_id);
                }
            }
        } else {
            for ridx in 0..self.replaced_ids.len() {
                let replaced_id = self.replaced_ids[ridx];
                if replaced_id == GeneratorRef::INVALID_ID || replaced_id >= self.vertices.len() {
                    continue;
                }
                if replaced_id >= self.n_revert_vertices {
                    self.unused.push(replaced_id);
                } else {
                    let refs = self.vertices[replaced_id].generator_refs;
                    for refg in refs {
                        if self.valid_generator_ref(refg) {
                            self.erase_cell_my_vertex_by_id(refg, replaced_id);
                        }
                    }
                }
            }
        }
        self.n_unused = self.unused.len();
    }

    /// For each cell touched by the insertion, records one of its vertices as a starting hint for the next insertion's graph traversal.
    fn assign_representative_vertices_to_cells(&mut self, default_id: usize) {
        let involved_front_id = self.involved_id_at(0);
        if involved_front_id == GeneratorRef::INVALID_ID {
            return;
        }
        if self.points[involved_front_id].my_vertices_ids.is_empty() {
            if default_id != GeneratorRef::INVALID_ID && default_id < self.vertices.len() {
                self.points[involved_front_id].my_vertices_ids.push(default_id);
            }
            return;
        }
        let new_representative_id = self.points[involved_front_id].my_vertices_ids[0];
        if new_representative_id == GeneratorRef::INVALID_ID || new_representative_id >= self.vertices.len() {
            return;
        }
        for involved_idx in 1..self.involved_refs.len() {
            let involved_id = self.involved_id_at(involved_idx);
            if involved_id == GeneratorRef::INVALID_ID || self.points[involved_id].my_vertices_ids.is_empty() {
                continue;
            }
            let representative_id = self.points[involved_id].my_vertices_ids[0];
            if representative_id == GeneratorRef::INVALID_ID || representative_id >= self.vertices.len() {
                continue;
            }
            if !self.vertices[representative_id].is_connected() {
                self.points[involved_id].my_vertices_ids[0] = new_representative_id;
            }
        }
    }

    #[inline(always)]
    /// Returns `true` if this vertex has a Side (bounding-plane) generator in slot 3, meaning it lies on the cube boundary.
    fn has_virtual_generators(&self, vtx: &Vertex<Scalar>) -> bool {
        let refg = vtx.resolved_generator_ref(3);
        !self.ref_is_real_point(refg)
    }

    /// Appends a zero-crossing record to `self.zeros` and registers its index in the owning point cell's `my_zero_points` list.
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

    /// Clears and fully rebuilds the diagram from new atom positions and radii, reusing existing buffer allocations.
    pub(crate) fn recalculate(
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
        self.build_vertices(self.points.len(), 0);

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

    /// Rebuilds the `my_vertices_ids` cache for all cells from scratch, processing all non-corner vertices.
    pub(crate) fn fill_all_my_vertices(&mut self) {
        self.fill_all_my_vertices_from(0, 8);
    }

    /// Assigns each non-corner finite vertex to the point cell(s) whose generator it belongs to. Vertices adjacent to a Side generator are assigned to whichever of their Point generators owns the opposite side of that boundary plane.
    pub(crate) fn fill_all_my_vertices_from(&mut self, from_point: usize, from_vertex: usize) {
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
            let is_corner = vtx.is_corner();
            if !self.has_virtual_generators(vtx) {
                for &refg in &refs {
                    if refg.kind == GeneratorKind::Point && refg.index < self.points.len() {
                        let idx = refg.index;
                        self.points[idx].my_vertices_ids.push(vi);
                        if from_point > 0 && self.points[idx].visited_as == 0 {
                            self.points[idx].visited_as = -1;
                            self.push_involved(GeneratorRef::new(GeneratorKind::Point, idx));
                        }
                    }
                }
            } else if !is_corner {
                for &refg in &refs {
                    if refg.kind == GeneratorKind::Side {
                        break;
                    }
                    if refg.kind == GeneratorKind::Point && refg.index < self.points.len() {
                        let idx = refg.index;
                        self.points[idx].my_vertices_ids.push(vi);
                        if from_point > 0 && self.points[idx].visited_as == 0 {
                            self.points[idx].visited_as = -1;
                            self.push_involved(GeneratorRef::new(GeneratorKind::Point, idx));
                        }
                    }
                }
            } else if !self.params.without_check {
                // Mirrors C++ assert path: corners should not be propagated to SASA ownership sets.
                return;
            }
        }
    }

    /// Populates `neighbours_ids` for every cell by walking its vertices and collecting distinct adjacent Point generators, using `visited_as` marks to avoid duplicates.
    pub(crate) fn fill_all_neighbours(&mut self) {
        for p in &mut self.points {
            p.neighbours_ids.clear();
            p.visited_as = -1;
        }

        for i in 0..self.points.len() {
            for vidx in 0..self.points[i].my_vertices_ids.len() {
                let vid = self.points[i].my_vertices_ids[vidx];
                if vid >= self.vertices.len() {
                    continue;
                }
                let v = self.vertices[vid];
                if v.is_corner() {
                    continue;
                }
                for g in (0..=3).rev() {
                    let gref = v.resolved_generator_ref(g);
                    if !gref.is_valid() || gref.kind != GeneratorKind::Point {
                        continue;
                    }
                    let gidx = gref.index;
                    if gidx == i || gidx >= self.points.len() {
                        continue;
                    }
                    if self.points[gidx].visited_as < i as i32 {
                        self.points[i].neighbours_ids.push(gidx);
                        self.points[gidx].visited_as = i as i32;
                    }
                }
            }
        }
        for p in &mut self.points {
            p.visited_as = 0;
        }
    }

    /// Recomputes all power-zero crossings on vertex edges from scratch, starting after the 8 cube-corner vertices.
    pub(crate) fn fill_all_zero_points(&mut self) {
        self.fill_all_zero_points_from(8, 0);
    }

    /// Scans every vertex edge and solves for the parameter(s) where the power value crosses zero, recording each crossing as a `ZeroPoint`. Only considers edges where the boundary generator is a real point cell and at least one endpoint has positive power.
    pub(crate) fn fill_all_zero_points_from(&mut self, from_vertex: usize, from_zero: usize) {
        if from_zero <= self.zeros.len() {
            self.zeros.truncate(from_zero);
        } else {
            self.zeros.clear();
        }
        for p in &mut self.points {
            p.my_zero_points.clear();
        }
        for vertex_index in from_vertex..self.n_vertices.min(self.vertices.len()) {
            if self.vertices[vertex_index].invalid {
                continue;
            }
            let current_boundary_ref = self.vertices[vertex_index].generator_refs[2];
            if !(current_boundary_ref.kind == GeneratorKind::Point && current_boundary_ref.index < self.points.len()) {
                continue;
            }
            let current_gen3 = self.vertices[vertex_index].generator_refs[3];
            let endpoint_start = if !self.ref_is_real_point(current_gen3) { 3 } else { 0 };
            let current_power_value = self.vertices[vertex_index].power_value;
            let current_position    = self.vertices[vertex_index].position;
            for endpoint_idx in endpoint_start..=3 {
                let endpoint_id = self.vertices[vertex_index].end_point_ids[endpoint_idx];
                if endpoint_id == GeneratorRef::INVALID_ID || endpoint_id <= vertex_index || endpoint_id >= self.n_vertices {
                    continue;
                }
                if current_power_value > Scalar::zero() {
                    let branch = endpoint_idx as i32;
                    let endpoint = &self.vertices[endpoint_id];
                    // Solve for t ∈ (0,1) where power crosses zero along this edge.
                    // Power is quadratic in t; coefficients use values at three points:
                    //   v2 = power at the current vertex (t=0),
                    //   v3 = power at the far endpoint (t=1),
                    //   v1 = power of the third generator evaluated at the reflection 2·current − endpoint.
                    // quot = 2*(v1 + v3 − 2*v2) is the second-order coefficient.
                    // The two roots sol1, sol2 are the candidate crossing parameters.
                    let v3 = endpoint.power_value;
                    let v2 = current_power_value;
                    let refg = self.vertices[vertex_index].generator_refs[if branch == 0 { 1 } else { 0 }];
                    let v1 = self.get_generator(refg).power(current_position * Scalar::from_f64(2.0).unwrap() - endpoint.position);
                    let quot = Scalar::from_f64(2.0).unwrap() * (v1 + v3 - Scalar::from_f64(2.0).unwrap() * v2);
                    let rootsq = (v1 - v3) * (v1 - v3) - Scalar::from_f64(4.0).unwrap() * quot * v2;
                    if rootsq <= Scalar::zero() {
                        continue;
                    }
                    if quot < self.power_err && v1 >= Scalar::zero() && v2 >= Scalar::zero() && v3 >= Scalar::zero() {
                        continue;
                    }
                    let rootquot = rootsq.sqrt() / quot;
                    let min = (v1 - v3) / quot;
                    let sol1 = min + rootquot;
                    let sol2 = min - rootquot;
                    let ep_positive = v3 > Scalar::zero();
                    if sol1 > Scalar::zero() && sol1 < Scalar::one() {
                        if ep_positive {
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
                } else {
                    let endpoint = &self.vertices[endpoint_id];
                    let ep_power = endpoint.power_value;
                    let ep_position = endpoint.position;
                    if ep_power > Scalar::zero() {
                        let branch = endpoint_idx as i32;
                        // Solve for t ∈ (0,1) where power crosses zero along this edge.
                        // Power is quadratic in t; coefficients use values at three points:
                        //   v2 = power at the current vertex (t=0),
                        //   v3 = power at the far endpoint (t=1),
                        //   v1 = power of the third generator evaluated at the reflection 2·current − endpoint.
                        // quot = 2*(v1 + v3 − 2*v2) is the second-order coefficient.
                        // The two roots sol1, sol2 are the candidate crossing parameters.
                        let v3 = ep_power;
                        let v2 = current_power_value;
                        let refg = self.vertices[vertex_index].generator_refs[if branch == 0 { 1 } else { 0 }];
                        let v1 = self.get_generator(refg).power(current_position * Scalar::from_f64(2.0).unwrap() - ep_position);
                        let quot = Scalar::from_f64(2.0).unwrap() * (v1 + v3 - Scalar::from_f64(2.0).unwrap() * v2);
                        let rootsq = (v1 - v3) * (v1 - v3) - Scalar::from_f64(4.0).unwrap() * quot * v2;
                        if rootsq <= Scalar::zero() {
                            continue;
                        }
                        if quot < self.power_err && v1 >= Scalar::zero() && v2 >= Scalar::zero() && v3 >= Scalar::zero() {
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

    #[inline(always)]
    /// Returns a mutable reference to the point cell at the given index.
    pub(crate) fn get_cell_mut(&mut self, id: usize) -> &mut Cell<Scalar> {
        &mut self.points[id]
    }

    /// Looks up a generator by reference, returning the corresponding point cell or side-plane cell.
    pub(crate) fn get_generator(&self, aref: GeneratorRef) -> &Cell<Scalar> {
        match aref.kind {
            GeneratorKind::Point => &self.points[aref.index],
            GeneratorKind::Side => &self.side_generators[aref.index],
        }
    }

    #[inline(always)]
    /// Returns `true` if the reference is a `Point` kind with an index within the point table.
    pub(crate) fn ref_is_real_point(&self, aref: GeneratorRef) -> bool {
        aref.kind == GeneratorKind::Point && aref.index < self.points.len()
    }

    /// Returns `true` if the reference's index is in bounds for its kind.
    pub(crate) fn valid_generator_ref(&self, aref: GeneratorRef) -> bool {
        match aref.kind {
            GeneratorKind::Point => aref.index < self.points.len(),
            GeneratorKind::Side => aref.index < self.side_generators.len(),
        }
    }

    /// Returns `true` if the zero-crossing record refers to in-bounds, mutually connected vertices.
    pub(crate) fn zero_point_valid(&self, zp: &ZeroPoint<Scalar>) -> bool {
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
        let to_id = from.resolved_endpoint_id(zp.branch as usize);
        if to_id == GeneratorRef::INVALID_ID || to_id >= self.n_vertices || to_id >= self.vertices.len() {
            return false;
        }
        !self.vertices[to_id].invalid
    }

    #[inline(always)]
    /// Linearly interpolates between the two endpoint positions of a zero-crossing edge: `from*(1-pos) + to*pos` where `pos ∈ [0,1]`.
    pub(crate) fn zero_point_pos(&self, zp: &ZeroPoint<Scalar>) -> Vector3<Scalar> {
        if !self.zero_point_valid(zp) {
            return Vector3::zeros();
        }
        let from = &self.vertices[zp.from_id];
        let to = &self.vertices[from.resolved_endpoint_id(zp.branch as usize)];
        to.position * zp.pos - from.position * (zp.pos - Scalar::one())
    }

    #[inline(always)]
    /// Executes the `error` step of the power-diagram algorithm/state machine.
    pub(crate) fn error(f: Scalar) -> Scalar {
        // C++ semantics rely on std::numeric_limits<T>::min() (minimum positive normal),
        // making error() strictly non-negative and approximately |f|*eps for practical ranges.
        f.abs() * Scalar::default_epsilon()
    }
}

#[inline(always)]
/// Returns the `n`-th element of the sequence 0..4 with `without` removed: used to select 3 of the 4 generator slots when one is reserved for the new generator.
pub(crate) fn nth(n: i32, without: i32) -> i32 {
    // Enumerates {0,1,2,3} \ {without} in order:
    // indices below `without` are unchanged; indices at or above shift up by 1.
    n + i32::from(without <= n)
}

/// Computes a radius-aware bounding box and applies default padding used by initial clipping.
pub(crate) fn get_bounding_box<Scalar>(
    pos: &[Vector3<Scalar>],
    strength: &[Scalar],
) -> (Vector3<Scalar>, Vector3<Scalar>)
where
    Scalar: RealField + Copy,
{
    let additional_cube_size = Scalar::from_f64(2.0f64.powf(1.0 / 3.0) - 1.0).unwrap();
    get_bounding_box_with_padding(pos, strength, additional_cube_size)
}

/// Computes a radius-aware bounding box with caller-specified extra padding.
pub(crate) fn get_bounding_box_with_padding<Scalar>(
    pos: &[Vector3<Scalar>],
    strength: &[Scalar],
    additional_cube_size: Scalar,
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
    low += (low - center) * additional_cube_size;
    high += (high - center) * additional_cube_size;
    (low, high)
}
