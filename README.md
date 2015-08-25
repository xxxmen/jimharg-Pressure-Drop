README - PRESSURE DROP CALCULATOR V0.02
=======================================

## What is it?

Pressure Drop Calculator (PDC) is an R script for calculating the pressure drop 
in piping systems transporting incompressible fluids (i.e. liquids). PDC also
includes functions for pump sizing and calculating NPSHa.


## How do I use it?

Configurable parameters are defined in section 1. This includes process
material properties and system geometries. Enter the values for your system.

The input data to PDC can be vectorised, allowing for large numbers of pressure
drop calculations to be done at once. E.g, setting pipe diameter as such:

pipe.dia <- c(0.025, 0.050, 0.075, 0.100)

Will determine the pressure drop at the 4 different pipe sizes.

Section 3 contains executed statements. These should be modified or added to
to suit your application and desired output. The default is to display a summary
of results in a table.


## The Latest Version

The latest version is available at: https://github.com/jimharg/Pressure-Drop


## Licensing

Pressure drop Calculator.R is free software licensed under the GNU GPL v3. 
For full details see license.md.


## Contacts

Bug reports to 'jamespaulhargreaves at gmail dot com'


## File Manifest

Pressure Drop Calculator.R

## Changelog

###[0.2] - 25/08/15

- CalcDarcyFactor() now accepts vector arguments. This means the script
  can now perform multiple pressure drop calculations in one go.
 
- KFitting is now manually input. Automated calculation is too 'messy' to
  implement cleanly. This will hopefully be revisted in a future release.

- CalcPumpPower() variable nomenclature is now consistent with the rest of
  the script.

- ConvertH2P() and ConvertP2H() functions have been added for calculating pressure
  from fluid head, and fluid head from pressure, respectively.

