// fibonacci circuit with no public inputs but checks on the first (and last?) row.
use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::PrimeField;
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::MatrixRowSlices;
use std::borrow::Borrow;

pub struct FibonacciAir {}

impl<F> BaseAir<F> for FibonacciAir {
    fn width(&self) -> usize {
        NUM_FIBONACCI_COLS
    }
}

impl<AB: AirBuilder> Air<AB> for FibonacciAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();

        let local: &FibonacciRow<AB::Var> = main.row_slice(0).borrow();
        let next: &FibonacciRow<AB::Var> = main.row_slice(1).borrow();

        let mut when_first_row = builder.when_first_row();

        when_first_row.assert_zero(local.left);
        when_first_row.assert_one(local.right);

        let mut when_transition = builder.when_transition();

        // a' <- b
        when_transition.assert_eq(next.left, local.right + local.left);

        // b' <- a + b
        when_transition.assert_eq(next.right, next.left + local.right);

        //        builder.when_last_row().assert_eq(local.left, x);
        //       builder.when_last_row().assert_eq(local.right, y);
    }
}

pub fn generate_trace_rows<F: PrimeField>(a: u64, b: u64, n: usize) -> RowMajorMatrix<F> {
    let mut trace =
        RowMajorMatrix::new(vec![F::zero(); n * NUM_FIBONACCI_COLS], NUM_FIBONACCI_COLS);

    let (prefix, rows, suffix) = unsafe { trace.values.align_to_mut::<FibonacciRow<F>>() };
    assert!(prefix.is_empty(), "Alignment should match");
    assert!(suffix.is_empty(), "Alignment should match");
    assert_eq!(rows.len(), n);

    rows[0] = FibonacciRow::new(F::from_canonical_u64(a), F::from_canonical_u64(b));

    for i in 1..n {
        rows[i].left = rows[i - 1].right + rows[i - 1].left;
        rows[i].right = rows[i].left + rows[i - 1].right;
    }

    trace
}

const NUM_FIBONACCI_COLS: usize = 2;

pub struct FibonacciRow<F> {
    pub left: F,
    pub right: F,
}

impl<F> FibonacciRow<F> {
    fn new(left: F, right: F) -> FibonacciRow<F> {
        FibonacciRow { left, right }
    }
}

impl<F> Borrow<FibonacciRow<F>> for [F] {
    fn borrow(&self) -> &FibonacciRow<F> {
        debug_assert_eq!(self.len(), NUM_FIBONACCI_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<FibonacciRow<F>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

#[cfg(test)]
pub mod tests {
    use ark_std::UniformRand;
    use folding_schemes::arith::Arith;
    use rand::SeedableRng;

    use super::{generate_trace_rows, FibonacciAir, NUM_FIBONACCI_COLS};
    use ark_bn254::Fr;
    use plonky3_ccs::arkfield::ArkField;
    use plonky3_ccs::{
        air_constraints_to_air_polynomials, air_trace_to_z,
        air_transition_polynomials_to_unique_polynomial, build_ccs_matrices,
        build_multisets_and_c_coefficients, compile_circuit_cs, derive_ccs_constants,
    };

    use folding_schemes::arith::ccs::CCS;

    #[test]
    fn test_fib_no_public() {
        let k = 2;
        let n = 2usize.pow(k);
        let mut air = FibonacciAir {};
        let num_public_values = 0;
        let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(42);
        let trace = generate_trace_rows::<ArkField<Fr>>(0, 1, n);
        let builder = compile_circuit_cs::<Fr, _>(&mut air, num_public_values);
        let air_polynomials =
            air_constraints_to_air_polynomials(&builder.constraints, &trace, n, NUM_FIBONACCI_COLS);
        let z = air_trace_to_z(&trace, &air_polynomials.boundary_polynomials);

        let r = Fr::rand(&mut rng);
        let final_poly = air_transition_polynomials_to_unique_polynomial(
            &r,
            &air_polynomials.transition_polynomials,
        );
        let ccs_constants = derive_ccs_constants(&trace, NUM_FIBONACCI_COLS, &final_poly).unwrap();
        let air_ccs_matrices = build_ccs_matrices(
            &ccs_constants,
            &air_polynomials,
            num_public_values,
            &z,
            NUM_FIBONACCI_COLS,
        );
        let (multisets, c_coeffs) = build_multisets_and_c_coefficients(&ccs_constants, &final_poly);

        let ccs = CCS {
            m: ccs_constants.m + air_ccs_matrices.n_boundary_and_constants_constraints_rows,
            n,
            l: 2,
            t: ccs_constants.t,
            q: multisets.len(),
            d: 1, // max degree in each variable
            s: 3,
            s_prime: 3,
            M: air_ccs_matrices.matrices,
            S: multisets,
            c: c_coeffs,
        };

        let res = ccs.check_relation(&z);
        assert!(res.is_ok());
    }
}
