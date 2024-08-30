A [plonky3](https://github.com/Plonky3/Plonky3) to [CCS](https://eprint.iacr.org/2023/552.pdf) converter. 

## Notes

[plonky3](https://github.com/Plonky3/Plonky3) is a popular toolkit for building SNARK and STARK proofs. Plonky3's underlying arithmetization uses an Algebraic Intermediate Representation (AIR) of constraints. In the context of a proliferation of such arithmetization schemes, Setty et al. (2023) proposed a customizable constraint system (CCS) that encompasses all three major arithmetization techniques (R1CS, plonkish, AIR) under a single framework.

This library allows you to port your Plonky3 AIR circuits to CCS.

Note that this library does not support periodic constraints (i.e. constraints applied to specific subsets of rows). Indeed, the way a CCS built from an AIR containing periodic constraints is processed may vary depending on the specific downstream SNARK used.

This library does not provide any SNARK backend for the constructed CCS.

## Acknowledgements

- Most of the code used for extracting polynomials from Plonky3 AIR circuits was inspired by (and in some parts, directly copied from) [ed255](https://github.com/ed255)'s Plonky3 frontend for Halo2, available [here](https://github.com/privacy-scaling-explorations/halo2/tree/main/p3_frontend).
- [ed255](https://github.com/ed255) provided valuable advice on how to proceed. Thanks to him!
