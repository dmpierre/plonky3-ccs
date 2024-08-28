/// `ArkField` wraps an `arkworks` field element, for which we implement both plonky3 `Field` and `AbstractField` traits
mod serde;
use ark_ff::BigInteger;
use std::{
    fmt::Display,
    hash::Hash,
    iter::{Product, Sum},
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

use ark_ff::PrimeField;
use num_bigint::BigUint;
use p3_field::{AbstractField, Field, Packable, PrimeField as p3PrimeField, PrimeField64};

#[derive(Default, Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ArkField<F: PrimeField>(pub F);

impl<F: PrimeField> Display for ArkField<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0.to_string())
    }
}

impl<F: PrimeField> MulAssign for ArkField<F> {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}

impl<F: PrimeField> Product for ArkField<F> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|val| val.0).product())
    }
}

impl<F: PrimeField> Sum for ArkField<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|value| value.0).sum())
    }
}

impl<F: PrimeField> Add for ArkField<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<F: PrimeField> AddAssign for ArkField<F> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl<F: PrimeField> Sub for ArkField<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<F: PrimeField> SubAssign for ArkField<F> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl<F: PrimeField> Neg for ArkField<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl<F: PrimeField> Mul for ArkField<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl<F: PrimeField> Div for ArkField<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0 / rhs.0)
    }
}

impl<F: PrimeField> Field for ArkField<F> {
    type Packing = Self;

    fn try_inverse(&self) -> Option<Self> {
        let inv = self.0.inverse();
        inv.map(|val| Self(val))
    }

    fn order() -> BigUint {
        BigUint::from_bytes_be(&F::MODULUS.to_bytes_be())
    }
}

impl<F: PrimeField> Packable for ArkField<F> {}

impl<F: PrimeField> p3PrimeField for ArkField<F> {
    fn as_canonical_biguint(&self) -> BigUint {
        let bg = self.0.into_bigint();
        BigUint::from_bytes_be(&bg.to_bytes_be())
    }
}

impl<F: PrimeField> AbstractField for ArkField<F> {
    type F = Self;

    fn zero() -> Self {
        Self(F::zero())
    }

    fn one() -> Self {
        Self(F::one())
    }

    fn two() -> Self {
        Self(F::from(2 as u8))
    }

    fn neg_one() -> Self {
        Self(-F::one())
    }

    fn from_f(f: Self::F) -> Self {
        Self(f.0)
    }

    fn from_bool(b: bool) -> Self {
        Self(F::from(b))
    }

    fn from_canonical_u8(n: u8) -> Self {
        Self(F::from(n))
    }

    fn from_canonical_u16(n: u16) -> Self {
        Self(F::from(n))
    }

    fn from_canonical_u32(n: u32) -> Self {
        Self(F::from(n))
    }

    fn from_canonical_u64(n: u64) -> Self {
        Self(F::from(n))
    }

    fn from_canonical_usize(n: usize) -> Self {
        Self(F::from(n as u64))
    }

    fn from_wrapped_u32(n: u32) -> Self {
        Self(F::from(n))
    }

    fn from_wrapped_u64(n: u64) -> Self {
        Self(F::from(n))
    }

    fn generator() -> Self {
        Self(F::GENERATOR)
    }
}

// From https://github.com/privacy-scaling-explorations/halo2/blob/26fc7c279b7863a4673325517c38079e84c0f2d1/p3_frontend/src/fwrap.rs#L241C1-L254C2
// HACK: In general an `FWrap` will need more than 64 bits.  This trait is only implemented in
// order to use `FWrap` with witness generation from plonky3 that requries this trait but doesn't
// use the order.  Do not use an `ff::PrimeField` on a circuit that requires a 64 bit prime field
// (i.e. relies on the `ORDER_U64` value), only use it on circuits that always assign less than 64
// bit values on the field elements.
impl<F: PrimeField + Hash + Ord> PrimeField64 for ArkField<F> {
    const ORDER_U64: u64 = 0;

    fn as_canonical_u64(&self) -> u64 {
        self.as_canonical_biguint()
            .try_into()
            .expect("field fits in u64")
    }
}
