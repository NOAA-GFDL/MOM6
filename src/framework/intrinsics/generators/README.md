# Intrinsic constant generators

These tools generate lookup tables or polynomial coefficients for the Fortran
implementations in the parent directory. They are not needed to build MOM6.

Run the following commands from `src/framework/intrinsics/`.

## Logarithm tables

```sh
python3 generators/generate_MOM_log_data_n128.py > MOM_log_data_n128.F90
```

The default reciprocal-first search requires NumPy with extended-precision
`longdouble`. Final candidates are evaluated with Python's `Decimal` arithmetic.
Use `--help` for the alternative log-first strategy and search options.

## Table-reduced logarithm polynomial

```sh
sollya generators/generate_MOM_log_remez.sollya
```

Prints the coefficients used by `log1p_tail_remez_4` in `MOM_log.F90`.
The recipe approximates `(log(1+x)-x)/x^2` on `[-0.004, 0.004]`.

## Near-one logarithm polynomial

```sh
sollya generators/generate_MOM_log_near.sollya
```

Prints candidate coefficients for `log1p_compensated_near`. This reconstructs
the fitting problem but does **not** reproduce the current coefficients exactly.
It does not update the Fortran source. Validate the complete compensated
evaluation before adopting new coefficients.
