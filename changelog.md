# Changelog
All notable changes to HVCWISE_TEA will be documented in this file.


## Unreleased
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [add_tests] - 2024-09-03
### Added
- KPI computation with Matlab code

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