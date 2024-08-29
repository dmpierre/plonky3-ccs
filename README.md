A [plonky3](https://github.com/Plonky3/Plonky3) to [CCS](https://eprint.iacr.org/2023/552.pdf) converter. 

## Notes

[plonky3](https://github.com/Plonky3/Plonky3) is a popular toolkit for building SNARK and STARK proofs. Plonky3's underlying arithmetization uses an algebraic intermediate representation (AIR) of constraints, another well-known technique. In the context of a proliferation of such arithmetization schemes, Setty et al. (2023) proposed a customizable constraint system (CCS), encompassing all three major arithmetization techniques (R1CS, plonkish, AIR) under a single framework. 

This library lets you port your plonky3 AIR circuits to CCS. 

## Acknowledgements 

- Most of the code used for extracting polynomials out of plonky3 air circuits has been inspired (if not copied) from [ed255](https://github.com/ed255)'s plonky3 frontend for halo2. It is available [here](https://github.com/privacy-scaling-explorations/halo2/tree/main/p3_frontend). 
- [ed255](https://github.com/ed255) gave me very valuable advices on how to proceed, thanks to him!


