# Changelog
All notable changes to HVCWISE_TEA will be documented in this file.


## Unreleased
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [add_tests] - 2025-02-07
### Added
- Number of parallel lines in AC branches

## [add_tests] - 2025-01-23
### Added
- Length in Macro-scenario template
### Changed
- Reliability data template
- Refactor build_availability_series.m
### Deprecated
- Remove mentions to Matlab: HVDCWISE_TEA is only compatible with Octave
### Fix
- Right line lengths considered in availability series generation

## [add_tests] - 2024-11-26
### Fix
- Space in macro-scenario

## [add_tests] - 2024-11-25
### Changed
- Replace "" by "_" in macro-scenario
- Update main.jl

## [add_tests] - 2024-11-22
### Changed
- No temporary file written in the matlab_tools folder during the process

## [add_tests] - 2024-11-15
### Added
- File main.jl as launching script for the user
### Fixed
- Reliability parameters for DC branch are the ones for DC cables (not DC OHL)
- Matlab/Octave launcher detection

## [add_tests] - 2024-10-31
### Changed
- Automatised running of the Matlab scripts

## [add_tests] - 2024-10-29
### Added
- Sheet 'gen_res' in the user inputs (model, cost, series)
### Changed
- Matrices sorted by alphabetic order in the .m file
### Fixed
- Emission factors matrix (gen ids ordering different from the user inputs ones)

## [add_tests] - 2024-10-10
### Changed
- Several economic indicators instead of SEW: Generation_cost, RES_curtailment_cost, Load_shedding_cost, Producers_surplus, Consumers_surplus

## [add_tests] - 2024-10-09
### Changed
- Storage considered in SEW computation

## [add_tests] - 2024-10-08
### Changed
- Template model.xlsx: removal of 'connection type'
- Template default_values.xlsx: update comments and values (some default values are now default models)
### Fixed
- Resistance pu conversion in DC

## [add_tests] - 2024-10-03
### Changed
- RES generators (ndgen) are the ones with generation_cost = 0
### Fixed
- KPI computation
- build_grid_model. convdc power rating in MW and branch_currents in MW

## [add_tests] - 2024-09-26
### Added
- Emission factors (only CO2) in KPI
- Units in KPI results
### Changed
- matlab tool files refactoring
### Fixed
- KPI computation when empty folders in simulation results 

## [add_tests] - 2024-09-03
### Added
- KPI computation with Matlab code
### Changed
- main_attributes in build_outputs_from_json()
### Fixed
- Solver double definition in example_tea.jl and run_study.jl
- test_load_case.jl

## [add_tests] - 2024-08-27
### Added
- Reliability data taken from the user inputs
- First post-processing file
- template folder (for user inputs template files)
### Fixed
- Flexible load attributes missing in build_grid_model()
- ndgen building in build_grid_model()

## [add_tests] - 2024-07-25
### Changed
- run_study() in run_study.jl enables to run a study from user inputs to simulation results. A user intervention is needed to run a Matlab code to generate the availability series.

## [add_tests] - 2024-07-23
### Added
- contingencies_generation (for N-1 only)

## [add_tests] - 2024-07-22
### Added
- build_simulation_inputs.jl
### Changed
- run_study.jl uses now the simulation inputs built by build_simulation_inputs() and splits the yearly problem in short subproblems (ex: 1 week)

## [add_tests] - 2024-05-28
### Added
- build_outputs_from_csv()
- build_outputs_from_json()
### Fixed
- parse_data() used in load_case() was deleted in the multi-period branch

## [add_tests] - 2024-05-23
### Added
- test05 (1 gen & 1 load & 1 storage) implemented & validated
- run_study.jl
- Raw results saving in JSON

## [add_tests] - 2024-04-30
### Added
- Excel user interface for inputs & outputs
### Changed
- Test results saved in the untracked folder "output"
### Fixed
- test01 (1 gen & 1 load)
- test02 (1 ndgen & 1 flexible load)
- test04 (1 DC line & 2 converters)
- test05 (1 gen & 1 load & 1 storage) implemented but not validated

## [add_tests] - 2024-04-11
### Changed
- load_case.jl: conversion from MW (csv inputs) to per unit thanks to baseMVA
- test01: expected cost updated and use of MW values
### Added
- test05: storage with dispatchable generator & fixed load

## [add_tests] - 2024-01-24
### Added
- load_case.jl
- test_load_case.jl
- test folders