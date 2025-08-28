use crate::field_arithmetic::*;
use ndarray::{Array1, Array2};

fn gaussian_elimination(
    a: &Array2<u128>,
    b: &Array1<u128>,
    field_size: usize,
) -> Option<EliminationResult> {
    let (rows, cols) = a.dim();

    let mut binary_a = Array2::zeros((rows, cols));
    let mut field_b = b.clone();

    for i in 0..rows {
        for j in 0..cols {
            binary_a[[i, j]] = if a[[i, j]] == 0 { 0 } else { 1 };
        }
    }

    let mut pivot_cols = Vec::new();
    let mut current_row = 0;

    for col in 0..cols {
        let pivot_row = find_pivot(&binary_a, current_row, rows, col)?;

        if let Some(pivot) = pivot_row {
            if pivot != current_row {
                swap_rows(&mut binary_a, current_row, pivot, cols);
                field_b.swap(current_row, pivot);
            }

            pivot_cols.push(col);

            eliminate_col(
                &mut binary_a,
                &mut field_b,
                current_row,
                col,
                rows,
                cols,
                field_size,
            );

            current_row += 1;
        }
    }

    for row in current_row..rows {
        if field_b[row] != 0 {
            return None;
        }
    }

    let free_cols: Vec<usize> = (0..cols)
        .filter(|&col| !pivot_cols.contains(&col))
        .collect();

    Some(EliminationResult {
        pivot_cols,
        free_cols,
        binary_matrix: binary_a,
        field_vector: field_b,
    })
}

fn eliminate_col(
    matrix: &mut Array2<u128>,
    field_vector: &mut Array1<u128>,
    pivot_row: usize,
    col: usize,
    num_rows: usize,
    num_cols: usize,
    field_size: usize,
) {
    let rows_to_eliminate: Vec<usize> = (0..num_rows)
        .filter(|&row| row != pivot_row && matrix[[row, col]] == 1)
        .collect();

    if rows_to_eliminate.is_empty() {
        return;
    }
    let pivot_field_value = field_vector[pivot_row];
    const BATCH_SIZE: usize = 8;
    for batch in rows_to_eliminate.chunks(BATCH_SIZE) {
        for &row in batch {
            for j in (0..num_cols).step_by(8) {
                matrix[[row, j]] ^= matrix[[pivot_row, j]];
                matrix[[row, j + 1]] ^= matrix[[pivot_row, j + 1]];
                matrix[[row, j + 2]] ^= matrix[[pivot_row, j + 2]];
                matrix[[row, j + 3]] ^= matrix[[pivot_row, j + 3]];
                matrix[[row, j + 4]] ^= matrix[[pivot_row, j + 4]];
                matrix[[row, j + 5]] ^= matrix[[pivot_row, j + 5]];
                matrix[[row, j + 6]] ^= matrix[[pivot_row, j + 6]];
                matrix[[row, j + 7]] ^= matrix[[pivot_row, j + 7]];
            }
            field_vector[row] = gf_add(field_vector[row], pivot_field_value, field_size);
        }
    }
}

pub fn solve_binary_system(
    a: &Array2<u128>,
    b: &Array1<u128>,
    field_size: usize,
) -> Option<Vec<u128>> {
    let elimination_result = gaussian_elimination(a, b, field_size)?;
    construct_solution(&elimination_result, a.ncols(), field_size)
}

#[derive(Debug)]
pub struct EliminationResult {
    pub pivot_cols: Vec<usize>,
    pub free_cols: Vec<usize>,
    pub binary_matrix: Array2<u128>,
    pub field_vector: Array1<u128>,
}

fn find_pivot(
    matrix: &Array2<u128>,
    start_row: usize,
    end_row: usize,
    col: usize,
) -> Option<Option<usize>> {
    for row in start_row..end_row {
        if matrix[[row, col]] == 1 {
            return Some(Some(row));
        }
    }
    Some(None)
}

fn swap_rows(matrix: &mut Array2<u128>, row1: usize, row2: usize, num_cols: usize) {
    for j in 0..num_cols {
        let temp = matrix[[row1, j]];
        matrix[[row1, j]] = matrix[[row2, j]];
        matrix[[row2, j]] = temp;
    }
}

fn construct_solution(
    result: &EliminationResult,
    num_vars: usize,
    field_size: usize,
) -> Option<Vec<u128>> {
    let mut solution = vec![0u128; num_vars];

    for &col in &result.free_cols {
        solution[col] = (rand::random::<u8>() & 1) as u128;
    }

    for (i, &pivot_col) in result.pivot_cols.iter().enumerate() {
        let mut value = result.field_vector[i];
        for j in (pivot_col + 1)..num_vars {
            if result.binary_matrix[[i, j]] == 1 {
                value = gf_add(value, solution[j], field_size);
            }
        }
        solution[pivot_col] = value;
    }

    Some(solution)
}

pub fn gf2_matrix_rank(matrix: &Array2<u128>) -> usize {
    let mut a = matrix.clone();
    let (rows, cols) = a.dim();
    let mut rank = 0;

    for col in 0..cols {
        if rank >= rows {
            break;
        }

        let mut pivot_row = None;
        for row in rank..rows {
            if a[[row, col]] == 1 {
                pivot_row = Some(row);
                break;
            }
        }

        if let Some(pivot) = pivot_row {
            if pivot != rank {
                for j in 0..cols {
                    let temp = a[[rank, j]];
                    a[[rank, j]] = a[[pivot, j]];
                    a[[pivot, j]] = temp;
                }
            }

            for row in 0..rows {
                if row != rank && a[[row, col]] == 1 {
                    for j in 0..cols {
                        a[[row, j]] ^= a[[rank, j]];
                    }
                }
            }
            rank += 1;
        }
    }

    rank
}
