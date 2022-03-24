# dbdModels
 Miscellaneous DBD models implemted in OpenFOAM
 
 ## Shyy model
 
Shyy's model was implemented in OpenFOAM by means of an fvOptions script, to use this implementation just place the script in the ```/system``` directory of the case. This code can be used with the solvers:

* simpleFoam
* pisoFoam
* pimpleFoam
* Or any other solver that can make use of fvOptions

The code was validated by comparing the velocity profiles with the measurements reported by Shyy at the positions denoted ST1 to ST4.

 ## Dörr & Kloker model
 
 Upcoming
 
## Requirements

To use these codes an installation of OpenFOAM is required. For the post processing ParaView is required https://www.paraview.org/. 

The codes were tested with OpenFOAM version 8.0 from the https://openfoam.org/ branch. But it should also work with the legacy versions and newer versions. 

## References

* Shyy model:
SHYY, Wei; JAYARAMAN, B.; ANDERSSON, A. Modeling of glow discharge-induced fluid dynamics. Journal of applied physics, 2002, vol. 92, no 11, p. 6434-6443 https://aip.scitation.org/doi/abs/10.1063/1.1515103
