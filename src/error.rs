use thiserror::Error;

#[derive(Error, Debug)]
pub enum SlampFssError {
    #[error("Matrix rank mismatch: expected {expected}, got {actual}")]
    RankMismatch { expected: usize, actual: usize },

    #[error("No solution found for linear system")]
    NoSolution,

    #[error("Matrix is not full rank")]
    NotFullRank,

    #[error("Invalid state: T_rb = 0")]
    InvalidState,

    #[error("Bad solution: insufficient rank")]
    BadSolution,

    #[error("Key generation failed after {attempts} attempts")]
    KeyGenFailed { attempts: usize },
}

pub type Result<T> = std::result::Result<T, SlampFssError>;