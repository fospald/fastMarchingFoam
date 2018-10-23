<p align="center">
  <a href="LICENSE" alt="GPLv3 license"><img src="https://img.shields.io/badge/license-GPLv3-brightgreen.svg" /></a>
  <a href="#" alt="no warranty"><img src="https://img.shields.io/badge/warranty-no-red.svg" /></a>
</p>

# fibergen

A FFT-based homogenization tool.

* FFT-based homogenization based on Lippmann-Schwinger equation with staggered grid approach \cite{SchneiderOspaldKabel2015:1}
* homogenization for linear elasticity, large deformations, Stokes flow and heat equation
* C++, OpenMP multiprocessing, XML + Python scripting interface
* laminate mixing for interfaces \cite{KabelMerkertSchneider2014,SchneiderOspaldKabel2015:2}
* mixed boundary conditions \cite{Kabel2016}
* generation of fibers distributions
* use of tetrahedrical mesh geometries
* arbitrarily many materials
* reading of raw CT data (gzip compressed)
* identification of homogenized material parameters
* ...


## Requirements

The code was developed using [OpenFOAM](https://www.openfoam.com/)-2.3.x.
It might or might not work with newer versions of OpenFOAM (please test and let me know).


## Installation

1. download source
```
git clone https://github.com/fospald/fastMarchingFoam.git
```
2. run Allwmake, which uses the wmake tool of OpenFOAM
```
sh Allwmake
```


## Run Demo




## Acknowledgements

[Felix Ospald](https://www.tu-chemnitz.de/mathematik/part_dgl/people/ospald) gratefully acknowledges financial support by the [German Research Foundation](http://www.dfg.de/en/) (DFG), [Federal Cluster of Excellence EXC 1075](https://www.tu-chemnitz.de/MERGE/) "MERGE Technologies for Multifunctional Lightweight Structures". Many thanks to [Matti Schneider](https://www.itm.kit.edu/cm/287_3957.php) for his helpful introduction to FFT-based homogenization and ideas regarding the ACG distribution.

