Forcing
=======
Data Override
-------
When running MOM6 with the Flexible Modelling System (FMS) coupler, forcing can be specified by a `data_table` file. This is particularly useful when running MOM6 with a data atmosphere, as paths to the relevent atmospheric forcing products (eg. JRA55-do or ERA5) can be provided here. Each item in the data table must be separated by a new line, and contains the following information:

| ``gridname``: The component of the model this data applies to. eg. `atm` `ocn` `lnd` `ice`.
| ``fieldname_code``: The field name according to the model component. eg. `salt`
| ``fieldname_file``: The name of the field within the source file. 
| ``file_name``: Path to the source file.
| ``interpol_method``: Interpolation method eg. `bilinear`
| ``factor``: A scalar by which to multiply the field ahead of passing it onto the model.

| 
The data table can be written in two formats: "legacy" or as a standard yaml file.

Legacy Format:
    "ATM", "t_bot",  "t2m", "./INPUT/2t_ERA5.nc", "bilinear", 1.0

Speficying a constant value:
    Rather than overriding with data from a file, one can also set a field to constant. To do this, pass empty strings to `fieldname_file` and `file_name`. The `factor` now corresponds to the override value. For example, the following sets the temperature at the bottom of the atmosphere to 290 Kelvin: 

.. code-block:: rst
    "ATM", "t_bot",  "", "", "bilinear", 290.0

More information can be found in the FMS data_override `readme file <https://github.com/NOAA-GFDL/FMS/tree/main/data_override>`_. 

.. toctree::
    :maxdepth: 2

    api/generated/pages/Solar_Radiation
    api/generated/pages/Tracer_Fluxes
