# Changelog
All notable changes to HVCWISE_TEA will be documented in this file.


## Unreleased
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [add_tests] - 2024-04-25
### Added
- Excel user interface for inputs & outputs
### Fixed
- test01 (1 gen & 1 load)
- almost test05 (1 gen & 1 load & 1 storage)

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