# FaEB

Fast Energy Barrier (FaEB) evaluation tool to assess diffusion in nanoporous materials based on Gemmi symmetrical energy grids.

## Installation

Check if you have `c++11` compiler installed (may work with other compiler but has been tested and mainly used with this compiler).

Tested for c++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0

Compilation that follows the rules set in the `Makefile`:
```
make all
```

## How to use it

```
./faeb structure/KAXQIL_clean_14.cif forcefield/UFF.def 298 12 Xe 0.12
```

You can replace the structure and forcefield files and the other parameters.

The results are printed in a comma separated (csv) format: 

```
structure, adsorption enthalpy (kJ/mol), Henry coefficient (mol/kg/Pa), Basin energy (kJ/mol), Transition state energy (kJ/mol), Time (s)
```
```
KAXQIL_clean_14,-44.626,0.0300725,-44.3181,-25.8181,1.78605
```

To get Barrier energy, you simply need to substract the TS energy by the basin energy. 

## Grid visualization

Please refer to a previous work to print and visualize the energy grids: https://github.com/coudertlab/GrAED

## Acknowledgement

This code has been developed during a PhD thesis co-financed by the CEA and Orano under the supervision of François-Xavier Coudert: https://github.com/fxcoudert

This code includes the library developed in Gemmi: 
https://github.com/project-gemmi/gemmi.git

## License

The MIT License (MIT)

Copyright (c) 2023 Emmanuel Ren

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.