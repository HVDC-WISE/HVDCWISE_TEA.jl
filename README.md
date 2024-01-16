# HVDCWISE_TEA.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HVDC-WISE.github.io/HVDCWISE_TEA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HVDC-WISE.github.io/HVDCWISE_TEA.jl/dev/)
[![Build Status](https://github.com/HVDC-WISE/HVDCWISE_TEA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HVDC-WISE/HVDCWISE_TEA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HVDC-WISE/HVDCWISE_TEA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HVDC-WISE/HVDCWISE_TEA.jl)

## Overview

HVDCWISE_TEA.jl is a Julia package to perform techno-economic analyses of HVDC expansion candidates for large scale, hybrid AC/DC power transmission networks.
A multi-period, hybrid AC/DC Optimal Power Flow problem is built in order to quantify the change is operational costs related to the HVDC expansion candidates.
The package builds upon the [PowerModels](https://github.com/lanl-ansi/PowerModels.jl), [PowerModelsMCDC](https://github.com/Electa-Git/PowerModelsMCDC.jl), [FlexPlan](https://github.com/Electa-Git/FlexPlan.jl) and [CbaOPF](https://github.com/Electa-Git/CbaOPF.jl) packages, using a similar structure.

To activate the virtual environment of this project (in a julia terminal): 

]
activate my\personnal\path\to\HVDCWISE_TEA.jl
instantiate

## Development

HVDCWISE_TEA.jl is early-stage, research-grade software under current development.
If you have suggestions for improvement, please contact us via the Issues page on the repository.

## Acknowledgements

This code is being developed as part of the HVDC-WISE project  under the European Union’s Horizon 2020 research and innovation programme (grant agreement no. 101075424).

## License

This code is provided under a [BSD 3-Clause License](/LICENSE.md).