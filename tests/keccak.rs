use ark_bn254::Fr;
use ark_poly::polynomial::Polynomial;
use ark_std::UniformRand;
use folding_schemes::arith::ccs::CCS;
use folding_schemes::arith::Arith;
use p3_keccak_air::NUM_KECCAK_COLS;
use p3_keccak_air::{generate_trace_rows, KeccakAir};
use plonky3_ccs::{
    air_constraints_to_air_polynomials, air_trace_to_z,
    air_transition_polynomials_to_unique_polynomial, arkfield::ArkField, build_ccs_matrices,
    build_multisets_and_c_coefficients, compile_circuit_cs, derive_ccs_constants,
};
use rand::random;
use rand::SeedableRng;

#[test]
fn test_keccak() {
    let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(42);
    let num_hashes = 4;
    let inputs = (0..num_hashes).map(|_| random()).collect::<Vec<_>>();
    let mut air = KeccakAir {};
    let trace = generate_trace_rows::<ArkField<Fr>>(inputs);

    let num_public_values = 0;
    let n = trace.values.len() / NUM_KECCAK_COLS;
    let builder = compile_circuit_cs::<Fr, _>(&mut air, num_public_values);
    let air_polynomials =
        air_constraints_to_air_polynomials(&builder.constraints, &trace, n, trace.width);

    let z = air_trace_to_z(&trace, &air_polynomials.boundary_polynomials);
    let r = Fr::rand(&mut rng);
    let final_poly = air_transition_polynomials_to_unique_polynomial(
        &r,
        &air_polynomials.transition_polynomials,
    );

    // change how constants are derived to account for the case where there is no public inputs
    let ccs_constants = derive_ccs_constants(&trace, trace.width, &final_poly).unwrap();
    let air_ccs_matrices = build_ccs_matrices(
        &ccs_constants,
        &air_polynomials,
        num_public_values,
        &z,
        trace.width,
    );
    let (multisets, c_coeffs) = build_multisets_and_c_coefficients(&ccs_constants, &final_poly);

    let m = ccs_constants.m + air_ccs_matrices.n_boundary_and_constants_constraints_rows;
    let ccs = CCS {
        m,
        n,
        l: num_public_values,
        t: ccs_constants.t,
        q: final_poly.terms.len(),
        d: final_poly.degree(),
        s: (m.ilog2() + 1) as usize,
        s_prime: (n.ilog2() + 1) as usize,
        M: air_ccs_matrices.matrices,
        S: multisets,
        c: c_coeffs,
    };

    let res = ccs.check_relation(&z);
    assert!(res.is_ok());
}
