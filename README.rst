## Global Scalable Paleo Landscape Evolution  / gospl


|    |    |
| --- | --- |
| Build Status | [![Linux/MacOS Build Status](https://travis-ci.org/Geodels/gospl.svg?branch=master)](https://travis-ci.org/Geodels/gospl) [![Coverage Status](https://coveralls.io/repos/github/Geodels/gospl/badge.svg?branch=master)](https://coveralls.io/github/Geodels/gospl?branch=master) |
| Latest release | [![Github release](https://img.shields.io/github/release/Geodels/gospl.svg?label=tag&colorB=11ccbb)](https://github.com/Geodels/gospl/releases) [![PyPI version](https://badge.fury.io/py/gospl.svg?colorB=cc77dd)](https://pypi.org/project/gospl)  [![Total alerts](https://img.shields.io/lgtm/alerts/g/Geodels/gospl.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Geodels/gospl/alerts/)  [![Updates](https://pyup.io/repos/github/Geodels/gospl/shield.svg)](https://pyup.io/repos/github/Geodels/gospl/) |
| Documentation | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Documentation Status](https://readthedocs.org/projects/gospl/badge/?version=latest)](https://gospl.readthedocs.io/en/latest/?badge=latest) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)  [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Geodels/gospl.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Geodels/gospl/context:python)|

### What's in the box?

**gospl** (pronounced: /ˈɡospel/, [ˈɡo̞s̠pe̞l]) is a scalable parallelised Python-based numerical model built to simulate paleo-landscapes and basins reconstruction at global scale.


The model accounts for hillslope processes (soil creep using linear diffusion), fluvial incision (stream power law), spatially and temporally varying tectonics (horizontal and vertical displacements) and climatic forces (temporal and spatial precipitation changes and/or sea-level fluctuations).


### Specs


The model is based on the following approaches:

+ an adaptation of the implicit, parallelisable method for calculating drainage area for both single (D8) and multiple flow direction (Dinf) from [Richardson et al., 2014],
+ the methods developed in _badlands_ [Salles et al., 2018] for marine sediment distribution,
+ a _PETSc_ layer similar to the one in _eSCAPE_ [Salles, 2019] for mesh partitioning and solvers.


### License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU  General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see _http://www.gnu.org/licenses/gpl-3.0.en.html_.

### Documentation & Installation

_https://gospl.readthedocs.io/_