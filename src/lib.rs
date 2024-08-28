extern crate alloc;
use ark_poly::Polynomial;
pub mod air;
pub mod arkfield;
mod symbolic_builder;
mod symbolic_expression;
mod symbolic_variable;
mod utils;

use ark_ff::PrimeField;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};
use arkfield::ArkField;
use core::panic;
use folding_schemes::utils::vec::SparseMatrix;
use p3_air::Air;
use p3_field::{Field, PrimeField as p3PrimeField};
use p3_matrix::dense::RowMajorMatrix;
use std::{collections::HashMap, hash::Hash, usize};
pub use symbolic_builder::*;
use symbolic_expression::{Location, SymbolicExpression};
use thiserror::Error;
use utils::naive_arkpoly_mul;

use crate::symbolic_variable::{Public, Query, SymbolicVariable, Var};

pub fn extract_public_values<F: Field>(
    e: &SymbolicExpression<F>,
) -> Option<((usize, Location), usize)> {
    use SymbolicExpression as SE;
    use SymbolicVariable as SV;
    let (mul_lhs, mul_rhs) = match e {
        SE::Mul(lhs, rhs) => (&**lhs, &**rhs),
        _ => return None,
    };
    let (cell_location, (sub_lhs, sub_rhs)) = match (mul_lhs, mul_rhs) {
        (SE::Location(location @ (Location::FirstRow | Location::LastRow)), SE::Sub(lhs, rhs)) => {
            (*location, (&**lhs, &**rhs))
        }
        _ => return None,
    };
    let (cell_column, public) = match (sub_lhs, sub_rhs) {
        (
            SE::Variable(SV(
                Var::Query(Query {
                    is_next: false,
                    column,
                }),
                _,
            )),
            SE::Variable(SV(Var::Public(Public { index }), _)),
        ) => (*column, *index),
        _ => return None,
    };
    Some(((cell_column, cell_location), public))
}

pub fn extract_constants<F: Field>(e: &SymbolicExpression<F>) -> Option<((usize, Location), F)> {
    use SymbolicExpression as SE;
    use SymbolicVariable as SV;
    let (mul_lhs, mul_rhs) = match e {
        SE::Mul(lhs, rhs) => (&**lhs, &**rhs),
        _ => return None,
    };

    // 
    let (cell_location, (sub_lhs, sub_rhs)) = match (mul_lhs, mul_rhs) {
        (SE::Location(location @ (Location::FirstRow | Location::LastRow)), SE::Sub(lhs, rhs)) => {
            (*location, (&**lhs, &**rhs))
        }
        // match arm handling the assert_zero() constant constraint
        (SE::Location(location @ (Location::FirstRow | Location::LastRow)), SE::Variable(SV(
                Var::Query(Query {
                    is_next: false,
                    column,
                }),
                _,
            ))) => return Some(((*column, *location), F::zero())), 
        _ => return None,
    };
    let (cell_column, constant) = match (sub_lhs, sub_rhs) {
        (
            SE::Variable(SV(
                Var::Query(Query {
                    is_next: false,
                    column,
                }),
                _,
            )),
            SE::Constant(value),
        ) => (*column, *value),
        _ => return None,
    };
    Some(((cell_column, cell_location), constant))
}

pub fn compile_circuit_cs<F, A>(
    air: &A,
    num_public_values: usize,
) -> SymbolicAirBuilder<ArkField<F>>
where
    F: PrimeField + Hash,
    A: Air<SymbolicAirBuilder<ArkField<F>>>,
{
    let mut builder = SymbolicAirBuilder::<ArkField<F>>::new(air.width(), num_public_values);
    air.eval(&mut builder);
    builder
}

/// Instantiate a single term, sparse multivariate in `n_variables` polynomial from a symbolic
/// variable
fn sym_variable_to_sparse_poly<F: PrimeField>(
    variable: &SymbolicVariable<ArkField<F>>,
    n_variables: usize,
) -> SparsePolynomial<F, SparseTerm> {
    let n_cols = n_variables / 2;
    match variable.0 {
        Var::Query(query) => {
            let var = query.is_next as usize * n_cols + query.column;
            let poly = SparsePolynomial::from_coefficients_vec(
                n_variables,
                vec![(F::one(), SparseTerm::new(vec![(var, 1)]))],
            );
            debug_assert!(poly.num_vars() == n_variables);
            return poly;
        }
        Var::Public(_) => {
            panic!("Public variables are not supported here!")
        }
    }
}

fn p3_field_element_to_ark_field_element<F: p3PrimeField, AF: PrimeField>(element: &F) -> AF {
    AF::from_be_bytes_mod_order(&element.as_canonical_biguint().to_bytes_be())
}

/// Instantiate a constant sparse multivariate in `n_variables` polynomial from a symbolic
/// constant
fn sym_constant_to_sparse_poly<F: p3PrimeField, AF: PrimeField>(
    constant: &F,
    n_variables_in_polynomial: usize,
) -> SparsePolynomial<AF, SparseTerm> {
    let element = p3_field_element_to_ark_field_element::<F, AF>(constant);
    SparsePolynomial::from_coefficients_vec(
        n_variables_in_polynomial,
        vec![(element, SparseTerm::new(vec![]))],
    )
}

/// Turns a symbolic expression into a polynomial. This is used notably to collect transition
/// constraints.
fn expr_to_polynomial<F: PrimeField + Hash>(
    e: &SymbolicExpression<ArkField<F>>,
    n_variables_in_polynomial: usize,
) -> SparsePolynomial<F, SparseTerm> {
    use SymbolicExpression as SE;
    // Number of variables in the polynomial that we will build from the provided constraint
    match e {
        SE::Variable(var) => sym_variable_to_sparse_poly(var, n_variables_in_polynomial),
        SE::Neg(e) => -expr_to_polynomial(e, n_variables_in_polynomial),
        SE::Sub(lhs, rhs) => {
            &expr_to_polynomial(lhs, n_variables_in_polynomial)
                - &expr_to_polynomial(rhs, n_variables_in_polynomial)
        }
        SE::Location(_) => {
            // assuming that location is always multiplying some other polynomial
            // e.g. air constraint is of the form: Mul(Location(FirstRow), Sub(Variable(...)))
            sym_constant_to_sparse_poly(&ArkField(F::one()), n_variables_in_polynomial)
        }
        SE::Constant(val) => sym_constant_to_sparse_poly(val, n_variables_in_polynomial),
        SE::Add(lhs, rhs) => {
            &expr_to_polynomial(lhs, n_variables_in_polynomial)
                + &expr_to_polynomial(rhs, n_variables_in_polynomial)
        }
        SE::Mul(lhs, rhs) => naive_arkpoly_mul(
            &expr_to_polynomial(lhs, n_variables_in_polynomial),
            &expr_to_polynomial(rhs, n_variables_in_polynomial),
        ),
    }
}

// Build polynomials from the transition constraints, expressed as symbolic expressions in plonky3
// air circuits
pub fn build_polynomials_from_constraints<F: Field, AF: PrimeField>(
    constraints: &Vec<SymbolicExpression<ArkField<AF>>>,
    n_variables_in_polynomial: usize,
) -> Vec<SparsePolynomial<AF, SparseTerm>> {
    let mut polys = Vec::<SparsePolynomial<AF, SparseTerm>>::with_capacity(constraints.len());
    for constraint in constraints {
        let poly = expr_to_polynomial(constraint, n_variables_in_polynomial);
        polys.push(poly);
    }
    polys
}

/// Air polynomials come in two different forms: transition polynomials and boundary polynomials
/// We do not support periodic constraints (i.e. polynomials applying to a subset of rows) since
/// this may vary following the specific downstream SNARK used.
/// `transition_polynomials` are used to check the correctness of rows AIR transitions, `boundary_polynomials`
/// are used to check the correctness of public input values.
/// `constant_polynomials` are similar to `boundary_polynomials` since they define constants
/// against which trace values should be equal to. 
/// TODO: merge together constant and boundary polynomials
pub struct AirPolynomials<F: PrimeField> {
    pub transition_polynomials: Vec<SparsePolynomial<F, SparseTerm>>,
    pub boundary_polynomials: HashMap<usize, F>,
    pub constant_polynomials: HashMap<usize, F>
}

/// Given a set of AIR constraints, returns an `AirPolynomials` struct, containing both transition
/// and boundary constraints. `n` is the number of rows in the air trace
/// TODO: see how to support constant polynomials with transition constraints
pub fn air_constraints_to_air_polynomials<F: PrimeField>(
    constraints: &Vec<SymbolicExpression<ArkField<F>>>,
    trace: &RowMajorMatrix<ArkField<F>>,
    n: usize,
    n_cols: usize,
) -> AirPolynomials<F> {
    // We assume that boundary polynomials are of degree one, and we represent those
    // using a pair `(wtns_idx, value)`. This will later create the constraint that
    // a specific position in the `z` vector should equal `value`
    let mut boundary_polynomials = HashMap::new();
    let mut constant_polynomials = HashMap::new();
    let mut transition_polynomials = vec![];
    for constraint in constraints {
        if let Some(((cell_col, cell_loc), public)) = extract_public_values(constraint) {
            match cell_loc {
                Location::FirstRow => {
                    let trace_pos = cell_col;
                    let value = trace.values[trace_pos];
                    boundary_polynomials.insert(trace_pos, value.0);
                }
                Location::LastRow => {
                    let trace_pos = n_cols * (n - 1) + cell_col;
                    let value = trace.values[trace_pos];
                    boundary_polynomials.insert(trace_pos, value.0);
                }
                Location::Transition => {
                    panic!(
                            "[NOT SUPPORTED] Found a public value in a transition constraint, cell_loc: {}, cell_col: {}, public: {}",
                            cell_loc, cell_loc, public
                        );
                }
            };
            continue;
        }

        if let Some(((cell_col, cell_loc), constant)) = extract_constants(constraint) {
            match cell_loc {
                Location::FirstRow => {
                    let trace_pos = cell_col;
                    constant_polynomials.insert(trace_pos, constant.0);
                }
                Location::LastRow => {
                    let trace_pos = n_cols * (n - 1) + cell_col;
                    constant_polynomials.insert(trace_pos, constant.0);
                }
                Location::Transition => {
                    panic!(
                            "[NOT SUPPORTED] Found a constant in a transition constraint, cell_loc: {}, cell_col: {}, constant: {}",
                            cell_loc, cell_loc, constant 
                        );
                }
            };
            continue;
        } else {
            let poly = expr_to_polynomial(constraint, n_cols * 2);
            transition_polynomials.push(poly);
        }
    }

    AirPolynomials {
        transition_polynomials,
        boundary_polynomials,
        constant_polynomials
    }
}

/// Builds the `z` vector which will satisfy a CCS
pub fn air_trace_to_z<F: PrimeField>(
    trace: &RowMajorMatrix<ArkField<F>>,
    boundary_polynomials: &HashMap<usize, F>,
) -> Vec<F> {
    let mut z = vec![];
    let mut pubs = vec![];

    for (trace_pos, value) in trace.values.iter().enumerate() {
        if boundary_polynomials.contains_key(&trace_pos) {
            pubs.push(value.0);
        } else {
            z.push(value.0);
        }
    }
    z.extend(pubs);
    z.push(F::one());
    z
}

/// This is used to turn several air polynomials into a single one.
/// The random element `r` is provided by a verifier. Described in
/// https://eprint.iacr.org/2023/552.pdf, p. 7
pub fn air_transition_polynomials_to_unique_polynomial<F: PrimeField>(
    r: &F,
    transition_polynomials: &Vec<SparsePolynomial<F, SparseTerm>>,
) -> SparsePolynomial<F, SparseTerm> {
    let n_variables = transition_polynomials[0].num_vars;
    let mut mul = F::one();
    let zero_poly = SparsePolynomial::from_coefficients_vec(n_variables, vec![]);
    transition_polynomials.iter().fold(zero_poly, |acc, poly| {
        let r_factor = SparsePolynomial::from_coefficients_vec(
            n_variables,
            vec![(mul, SparseTerm::new(vec![]))],
        );
        mul *= r;
        acc + naive_arkpoly_mul(&r_factor, &poly)
    })
}

#[derive(Debug)]
pub struct AirCCSConstants {
    pub m: usize,
    pub t: usize,
    pub d: usize,
    pub q: usize,
    pub n: usize,
}

#[derive(Debug, Error)]
pub enum AirToCCSError {
    #[error("Incorrect length for public values. It should be of the same length as `number_of_cols * 2`. Got length: {0}")]
    IncorrectPublicValuesLength(usize),
}

pub fn derive_ccs_constants<F: PrimeField>(
    trace: &RowMajorMatrix<ArkField<F>>,
    n_cols: usize,
    air_final_poly: &SparsePolynomial<F, SparseTerm>,
) -> Result<AirCCSConstants, AirToCCSError> {
    let (m, t, d, q) = (
        (trace.values.len() / trace.width) - 1,
        trace.width * 2,
        air_final_poly.degree(),
        air_final_poly.terms().len(),
    );
    let n = (m + 1) * t / 2 + 1;
    Ok(AirCCSConstants { m, t, d, q, n })
}

// Sticking to notation, allowing a non snake case for N, the number of non zero entries in the
// matrices
#[allow(non_snake_case)]
pub struct AirCCSMatrices<F: PrimeField> {
    pub matrices: Vec<SparseMatrix<F>>,
    pub N: usize,
    pub n_boundary_and_constants_constraints_rows: usize,
}

/// Builds the CCS matrices out of both boundary and transition constraints
pub fn build_ccs_matrices<F: PrimeField>(
    ccs_constants: &AirCCSConstants,
    air_polynomials: &AirPolynomials<F>,
    num_public_values: usize,
    z: &Vec<F>,
    n_cols: usize,
) -> AirCCSMatrices<F> {
    let n_boundary_constraints = air_polynomials.boundary_polynomials.keys().len();
    let n_constants_constraints = air_polynomials.constant_polynomials.keys().len();
    let  n_boundary_and_constants_constraints_rows= (n_boundary_constraints + n_constants_constraints).div_ceil(ccs_constants.t);
    let mut matrices: Vec<SparseMatrix<F>> = vec![];
    for _ in 0..ccs_constants.t {
        let mut coeffs = vec![];
        for _ in 0..ccs_constants.m +  n_boundary_and_constants_constraints_rows{
            coeffs.push(vec![]);
        }
        matrices.push(SparseMatrix {
            n_rows: ccs_constants.m + n_boundary_and_constants_constraints_rows,
            n_cols: ccs_constants.n,
            coeffs,
        })
    }

    // This is the algorithm described in lemma 3, with a few modifications
    // We modify the calculation of k_j and handle the case of an air trace with no public inputs
    let len_w_air = ccs_constants.n - num_public_values - 1;
    #[allow(non_snake_case)]
    let mut N = 0;
    for i in 0..ccs_constants.m {
        for j in 0..ccs_constants.t {
            if i == 0 && j < ccs_constants.t / 2 {
                if num_public_values > 0 {
                    matrices[j].coeffs[i].push((F::one(), j + len_w_air));
                } else {
                    matrices[j].coeffs[i].push((F::one(), i * ccs_constants.t / 2 + j));
                }
            } else if i == ccs_constants.m - 1 && j >= ccs_constants.t / 2 {
                if num_public_values > 0 {
                    matrices[j].coeffs[i].push((F::one(), j + len_w_air));
                } else {
                    matrices[j].coeffs[i].push((F::one(), i * ccs_constants.t / 2 + j));
                }
            } else {
                if num_public_values > 0 {
                    let k_j = ((ccs_constants.t as i128) * ((i as i128) - 1) / 2) + (j as i128);
                    matrices[j].coeffs[i].push((F::one(), k_j as usize));
                } else {
                    matrices[j].coeffs[i].push((F::one(), i * ccs_constants.t / 2 + j));
                }
            }
            N += 1
        }
    }

    // Add boundary constraints
    // Assumes that there is a single boundary constraint per idx
    // Assumes that boundary constraints are not applied to intermediary witness values
    for (i, (trace_idx, value)) in air_polynomials.boundary_polynomials.iter().enumerate() {
        let z_idx = if trace_idx < &n_cols {
            // converts trace idx to idx in z
            // vector
            *trace_idx + len_w_air
        } else {
            *trace_idx
        };
        let matrix_idx = i % ccs_constants.t;
        let row_idx = ccs_constants.m + i / ccs_constants.t;
        matrices[matrix_idx].coeffs[row_idx].push((F::one(), z_idx));
        matrices[matrix_idx].coeffs[row_idx].push((-*value, z.len() - 1));
    }

    // Add constant constraints
    for (i, (trace_idx, value)) in air_polynomials.constant_polynomials.iter().enumerate() {
        let z_idx = *trace_idx;
        let matrix_idx = i % ccs_constants.t;
        let row_idx = ccs_constants.m + i / ccs_constants.t;
        matrices[matrix_idx].coeffs[row_idx].push((F::one(), z_idx));
        matrices[matrix_idx].coeffs[row_idx].push((-*value, z.len() - 1));
    }

    AirCCSMatrices {
        matrices,
        N,
        n_boundary_and_constants_constraints_rows,
    }
}

pub fn build_multisets_and_c_coefficients<F: PrimeField>(
    ccs_constants: &AirCCSConstants,
    unique_polynomial: &SparsePolynomial<F, SparseTerm>,
) -> (Vec<Vec<usize>>, Vec<F>) {
    // deriving multisets
    let mut multisets: Vec<Vec<usize>> = vec![];
    for _ in 0..ccs_constants.q {
        multisets.push(vec![]);
    }

    let mut c_coeffs = vec![];
    for (i, term) in unique_polynomial.terms.iter().rev().enumerate() {
        c_coeffs.push(term.0);
        for (variable, power) in term.1.iter() {
            if *power > 0 {
                for _ in 0..*power {
                    multisets[i].push(*variable)
                }
            }
        }
    }

    (multisets, c_coeffs)
}

#[cfg(test)]
pub mod tests {
    use ark_ff::{PrimeField, Zero};
    use ark_poly::{
        multivariate::{SparsePolynomial, SparseTerm},
        DenseMVPolynomial, Polynomial,
    };
    use std::rc::Rc;

    use ark_bn254::Fr;

    use crate::{
        arkfield::ArkField, expr_to_polynomial, SymbolicExpression as SE, SymbolicVariable as SA,
    };

    pub(crate) fn check_expr_to_polynomial<F: PrimeField>(
        expr: &SE<ArkField<F>>,
        poly: SparsePolynomial<F, SparseTerm>,
        expected_str_repr: &str,
        expected_n_variables: usize,
        expected_n_terms: usize,
        eval_point: &Vec<F>,
        expected_eval_result: u128,
    ) {
        assert_eq!(expr.to_string(), expected_str_repr);
        assert_eq!(poly.num_vars(), expected_n_variables);
        assert_eq!(poly.terms().len(), expected_n_terms);
        assert_eq!(poly.evaluate(eval_point), F::from(expected_eval_result));
    }

    #[test]
    fn test_expr_to_polynomial() {
        let var_1 = Rc::new(SE::Variable::<ArkField<Fr>>(SA::new_query(false, 1)));
        let var_0_prime = Rc::new(SE::Variable(SA::new_query(true, 0)));
        let var_1_prime = Rc::new(SE::Variable::<ArkField<Fr>>(SA::new_query(true, 1)));

        let eval_point = &vec![Fr::from(2), Fr::from(10), Fr::from(2), Fr::from(10)];
        let num_cols = 2;
        let num_variables_poly = num_cols * 2;

        let expr = SE::Sub(var_1, var_0_prime.clone());
        let poly = expr_to_polynomial(&expr, num_variables_poly);
        check_expr_to_polynomial(&expr, poly, "(w1 - w0')", 4, 2, &eval_point, 8);

        let expr_2 = expr.clone() + expr.clone();
        let poly = expr_to_polynomial(&expr_2, num_variables_poly);
        check_expr_to_polynomial(
            &expr_2,
            poly,
            "((w1 - w0') + (w1 - w0'))",
            4,
            2,
            &eval_point,
            16,
        );

        let expr_3 = expr_2.clone() * expr.clone();
        let poly = expr_to_polynomial(&expr_3, num_variables_poly);
        check_expr_to_polynomial(
            &expr_3,
            poly,
            "((w1 - w0') + (w1 - w0')) * (w1 - w0')",
            4,
            3,
            &eval_point,
            128,
        );

        let expr_4 = SE::Add(var_1_prime, var_0_prime);
        let expr_5 = expr_3.clone() * expr_4;
        let poly = expr_to_polynomial(&expr_5, num_variables_poly);
        check_expr_to_polynomial(
            &expr_5,
            poly,
            "((w1 - w0') + (w1 - w0')) * (w1 - w0') * (w1' + w0')",
            4,
            6,
            &eval_point,
            1536,
        );

        let expr_zero = expr_5.clone() - expr_5.clone();
        let poly = expr_to_polynomial(&expr_zero, num_variables_poly);
        assert!(poly.is_zero());
    }
}