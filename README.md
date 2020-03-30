## Global Scalable Paleo Landscape Evolution  / gospl


[![PyPI](https://img.shields.io/pypi/v/gospl)](https://pypi.org/project/gospl/)  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![Build Status](https://travis-ci.org/Geodels/gospl.svg?branch=master)](https://travis-ci.org/Geodels/gospl) [![Updates](https://pyup.io/repos/github/Geodels/gospl/shield.svg)](https://pyup.io/repos/github/Geodels/gospl/) [![Coverage Status](https://coveralls.io/repos/github/Geodels/gospl/badge.svg?branch=master)](https://coveralls.io/github/Geodels/gospl?branch=master)

[![Documentation Status](https://readthedocs.org/projects/gospl/badge/?version=latest)](https://gospl.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


### What's in the box?

**gospl** (pronounced: /ˈɡospel/, [ˈɡo̞s̠pe̞l]) is a scalable parallelised Python-based numerical model built to simulate paleo-landscapes and basins reconstruction at global scale.


**gospl** is a parallel TIN-based landscape evolution model, built to simulate topography dynamic of the Earth over millions of years.

The model accounts for hillslope processes (soil creep using linear diffusion), fluvial incision (stream power law), spatially and temporally varying tectonics (horizontal and vertical displacements) and climatic forces (temporal and spatial precipitation changes and/or sea-level fluctuations).


### Specs


The model is based on the following approaches:

+ an adaptation of the implicit, parallelisable method for calculating drainage area for both single (D8) and multiple flow direction (Dinf) from [Richardson et al., 2014],
+ the methods developed in _badlands_ [Salles et al., 2018] for marine sediment distribution,
+ a _PETSc_ layer similar to the one in _eSCAPE_ [Salles, 2019] for mash partitioning and solvers.


### License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see _http://www.gnu.org/licenses/lgpl-3.0.en.html_.

### Documentation & Installation

_https://gospl.readthedocs.io/_


#### Artemis HPC Sydney

```bash
module load gcc/4.9.3  python/3.6.5 petsc-gcc-mpich/3.11.1

export F77FLAGS=-fPIC
export FCFLAGS=-fPIC

python3 setup.py install --user
```

#### Performance tests

+ Resolution distribution for 10 million points (min:7.6 km, max:10.3 km, mean: 9.1 km)

| PTS | CPUS | WALLTIME | NODES |
| --- | --- | --- | --- |
| 10612062 | 8 | 17:08:58 | 1 |
| 10612062 | 16 | 08:58:32 | 2 |
| 10612062 | 32 | 04:57:58 | 4 |
| 10612062 | 64 | 03:32:15 | 6 |
| 10612062 | 96 | 02:41:17 | 6 |
| 10612062 | 128 | 02:40:57 | 7 |

+ Resolution distribution for 17 million points (min:4.8 km, max:7.6 km, mean: 6.0 km)

| PTS | CPUS | WALLTIME | NODES |
| --- | --- | --- | --- |
| 17004184 | 64 | 07:28:41 | 4 |
| 17004184 | 96 | 06:38:11 | 4 |
| 17004184 | 128 | 05:29:51 | 6 |
| 17004184 | 144 | 04:59:49 | 6 |
| 17004184 | 168 | 04:31:14 | 7 |
| 17004184 | 192 | 06:23:02 | 7 |
