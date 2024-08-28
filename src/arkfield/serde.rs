use crate::arkfield::ArkField;
use ark_ff::BigInteger;
use ark_ff::PrimeField;
use serde::{de::Visitor, Deserialize, Serialize};
use std::marker::PhantomData;

impl<F: PrimeField> Serialize for ArkField<F> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let bg = self.0.into_bigint();
        serializer.serialize_bytes(&bg.to_bytes_be())
    }
}

struct ArkFieldVisitor<F: PrimeField> {
    _f: PhantomData<F>,
}

impl<'de, F: PrimeField> Visitor<'de> for ArkFieldVisitor<F> {
    type Value = ArkField<F>;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("A field element, represented as bytes")
    }
}

impl<'de, F: PrimeField> Deserialize<'de> for ArkField<F> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let val = PhantomData::<F> {};
        deserializer.deserialize_bytes(ArkFieldVisitor::<F> { _f: val })
    }
}
