use ark_ff::{PrimeField, Zero};
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};

fn extract_vars_powers(term: &SparseTerm) -> Vec<(usize, usize)> {
    let vars = term.vars();
    let powers = term.powers();
    vars.into_iter()
        .zip(powers.into_iter())
        .collect::<Vec<(usize, usize)>>()
}

// Adapted from `ark_poly`, see:
// https://github.com/arkworks-rs/algebra/blob/dcf73a5f9610ba9d16a3c8e0de0b3835e5e5d5e4/poly/src/polynomial/multivariate/sparse.rs#L327
pub fn naive_arkpoly_mul<F: PrimeField>(
    cur: &SparsePolynomial<F, SparseTerm>,
    other: &SparsePolynomial<F, SparseTerm>,
) -> SparsePolynomial<F, SparseTerm> {
    if cur.is_zero() || other.is_zero() {
        SparsePolynomial::zero()
    } else {
        let mut result_terms = Vec::new();
        for (cur_coeff, cur_term) in cur.terms.iter() {
            for (other_coeff, other_term) in other.terms.iter() {
                let mut term = extract_vars_powers(cur_term);
                term.extend(extract_vars_powers(other_term));
                result_terms.push((*cur_coeff * *other_coeff, SparseTerm::new(term)));
            }
        }
        SparsePolynomial::from_coefficients_slice(cur.num_vars, result_terms.as_slice())
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::Fr;
    use ark_ff::One;
    use ark_poly::multivariate::Term;
    use ark_poly::Polynomial;
    use ark_poly::{
        multivariate::{SparsePolynomial, SparseTerm},
        DenseMVPolynomial,
    };

    use super::naive_arkpoly_mul;

    #[test]
    pub fn test_naive_arkpoly_mul() {
        let x0 = SparsePolynomial::from_coefficients_vec(
            3,
            vec![(Fr::one(), SparseTerm::new(vec![(0, 1)]))],
        );
        let x1 = SparsePolynomial::from_coefficients_vec(
            3,
            vec![(Fr::one(), SparseTerm::new(vec![(1, 1)]))],
        );
        let x2 = SparsePolynomial::from_coefficients_vec(
            3,
            vec![(Fr::one(), SparseTerm::new(vec![(2, 1)]))],
        );

        // x0 * x1
        let x0_x1 = naive_arkpoly_mul(&x0, &x1);
        assert_eq!(x0_x1.degree(), 2);
        assert_eq!(x0_x1.terms().len(), 1);

        // x0^3
        let x0_pow3 = naive_arkpoly_mul(&naive_arkpoly_mul(&x0, &x0), &x0);
        assert_eq!(x0_pow3.degree(), 3);
        assert_eq!(x0_pow3.terms().len(), 1);

        // (x0 * x1 + x0^3) * (x1 + x0^3)
        let poly_1 = &x0_x1 + &x0_pow3;
        let poly_2 = &x1 + &x0_pow3;
        let res = naive_arkpoly_mul(&poly_1, &poly_2);
        assert_eq!(res.degree(), 6);
        assert_eq!(res.terms().len(), 4);

        let res_2 = naive_arkpoly_mul(&res, &x2);
        assert_eq!(res_2.degree(), 7);
        assert_eq!(res.terms().len(), 4);
    }
}
