<p align="center">
  <a href="LICENSE" alt="GPLv3 license"><img src="https://img.shields.io/badge/license-GPLv3-brightgreen.svg" /></a>
  <a href="#" alt="no warranty"><img src="https://img.shields.io/badge/warranty-no-red.svg" /></a>
</p>

# fastMarchingFoam

The shortest path algorithm on the surface of OpenFOAM meshes. 

Note: the code from the geodesic subdirectory is a clone of the [geodesic](https://code.google.com/archive/p/geodesic) project (the downloadable zip file "geodesic_cpp_03_02_2008.zip"). This project is licensend under the MIT License. See geodesic/readme.txt for further information.


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

