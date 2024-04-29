# Changelog
All notable changes to HVCWISE_TEA will be documented in this file.


## Unreleased
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [add_tests] - 2024-04-29
### Added
- Excel user interface for inputs & outputs
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