use thiserror::Error;

#[derive(Debug, Error, Clone)]
pub enum SasaError {
    #[error("power diagram exception")]
    MyException,
    #[error("identical point exception")]
    IdenticalPoint,
    #[error("vertices-full exception")]
    VerticesFull,
    #[error("axis too short for stable orientation test")]
    AxisTooShort,
    #[error("invalid get_ang input")]
    InvalidGetAngInput,
    #[error("invalid get_next input")]
    InvalidGetNextInput,
    #[error("surface vertex index out of bounds")]
    SurfaceVertexIndexOutOfBounds,
    #[error("degenerate circle ordering")]
    DegenerateCircleOrdering,
    #[error("rank out of range")]
    RankOutOfRange,
    #[error("atom index out of bounds")]
    AtomIndexOutOfBounds,
    #[error("ambiguous vertex power near tolerance")]
    AmbiguousVertexPower,
    #[error("neighbour index out of bounds")]
    NeighbourIndexOutOfBounds,
    #[error("invalid zero-point partner count")]
    InvalidZeroPointPartnerCount,
    #[error("partner indices out of range")]
    PartnerOutOfRange,
    #[error("invalid partner point counts")]
    InvalidPartnerPointCounts,
    #[error("volume node out of bounds")]
    VolumeNodeOutOfBounds,
    #[error("zero denominator in volume correction")]
    ZeroVolumeDenominator,
    #[error("invalid volume sign configuration")]
    InvalidVolumeSignConfiguration,
    #[error("shared-partner overflow")]
    SharedPartnerOverflow,
}
