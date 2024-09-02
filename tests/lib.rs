pub mod fib;

#[cfg(test)]
mod tests {

    use crate::fib::FibonacciRow;
    use ark_bn254::Fr;
    use p3_air::AirBuilder;
    use p3_field::Field;
    use plonky3_ccs::{arkfield::ArkField, SymbolicAirBuilder};
    use plonky3_ccs::{
        extract_constants, AirConstant, Location, SymbolicExpression as SE, SymbolicVariable as SV,
    };
    use std::u64;

    fn initialize_builder_row_and_variable(
        variable_is_next: bool,
        variable_col: usize,
        constant_value: u64,
    ) -> (
        SymbolicAirBuilder<ArkField<Fr>>,
        FibonacciRow<ArkField<Fr>>,
        SE<ArkField<Fr>>,
        SE<ArkField<Fr>>,
    ) {
        let mut builder = SymbolicAirBuilder::<ArkField<Fr>>::new(2, 0);
        let fib_row =
            FibonacciRow::new(ArkField(Fr::from(10 as u64)), ArkField(Fr::from(10 as u64)));
        let var = SE::Variable::<ArkField<Fr>>(SV::new_query(variable_is_next, variable_col));
        let constant = SE::Constant(ArkField(Fr::from(constant_value)));
        (builder, fib_row, var, constant)
    }

    // TODO:assert_eq! on the location as well
    fn check_air_constraint_to_constant(
        extracted_constant: Option<AirConstant<ArkField<Fr>>>,
        expected_col: usize,
        expected_value: u64,
    ) {
        let constant = extracted_constant.unwrap();
        assert_eq!(expected_col, constant.column);
        assert_eq!(ArkField(Fr::from(expected_value)), constant.value);
    }

    #[test]
    fn test_extract_zero_constant() {
        // extract non zero constant not next row, 1st column
        let (mut builder, row, var, constant) = initialize_builder_row_and_variable(false, 1, 2);
        let mut when_first_row = builder.when_first_row();
        when_first_row.assert_zero(row.right);
        let zero_constant = extract_constants(&builder.constraints[0]);
    }

    #[test]
    fn test_extract_non_zero_constant() {
        // extract non zero constant not next row, 1st columwn, when_first_row
        let (mut builder, _, var, constant) = initialize_builder_row_and_variable(false, 1, 1);
        let mut when_first_row = builder.when_first_row();
        when_first_row.assert_eq(constant.clone(), var.clone());
        when_first_row.assert_eq(var.clone(), constant.clone());
        let constant_left = extract_constants(&builder.constraints[0]);
        let constant_right = extract_constants(&builder.constraints[1]);
        check_air_constraint_to_constant(constant_left, 1, 1);
        check_air_constraint_to_constant(constant_right, 1, 1);

        // extract non zero constant next row, 1st column, when_first_row
        let (mut builder, row, var, constant) = initialize_builder_row_and_variable(true, 1, 2);
        let mut when_first_row = builder.when_first_row();
        when_first_row.assert_eq(row.right, var.clone());
        when_first_row.assert_eq(var.clone(), row.right);
        let constant_left = extract_constants(&builder.constraints[0]);
        let constant_right = extract_constants(&builder.constraints[1]);

        //check_air_constraint_to_constant(constant_left, 1, 2);
        //check_air_constraint_to_constant(constant_right, 1, 2);

        // extract non zero constant not next row, 1st columwn, when_last_row
        let (mut builder, row, var, constant) = initialize_builder_row_and_variable(false, 1, 3);
        let mut when_last_row = builder.when_last_row();
        when_last_row.assert_eq(row.right, var.clone());
        when_last_row.assert_eq(var.clone(), row.right);
        let constant_left = extract_constants(&builder.constraints[0]);
        let constant_right = extract_constants(&builder.constraints[1]);

        // extract non zero constant next row, 1st column, when_last_row
        let (mut builder, row, var, constant) = initialize_builder_row_and_variable(true, 1, 4);
        let mut when_last_row = builder.when_last_row();
        when_last_row.assert_eq(row.right, var.clone());
        when_last_row.assert_eq(var.clone(), row.right);
        let constant_left = extract_constants(&builder.constraints[0]);
        let constant_right = extract_constants(&builder.constraints[1]);
    }

    #[test]
    fn test_extract_public_values() {}
}
