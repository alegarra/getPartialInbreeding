# getPartialInbreeding and TtimesImP.jl

These programs allow computing of partial inbreeding coefficients (Garcia Cortes et al. 2010) and its posterior use in BLUPto estimates individual Inbreeding Depression Load as in Varona et al. 2019; Martinez-Castillero et al. 2021 

The programs need a renumbered (parents before offspring) pedigree with missing parents set to 0. A pedigree example is provided: the genealogy of bull Comet (Garcia-Cortes et al., 2010).



## getPartialInbreeding.f90

This fortran code can be compiled e.g. as `ifort -O3 -heap-arrays -fopenmp getPartialInbreeding.f90 -o getPartialInbreeding`.

And run using e.g.

```
./getPartialInbreeding pedigree_file [first last]
```

the first and last are optional integers and indicate that we want to compute the partial coefficients from animal "first" to animal "last"



García-Cortés, L. A., Martínez-Ávila, J. C., & Toro, M. A. (2010). Fine decomposition of the inbreeding and the coancestry coefficients by using the tabular method. Conservation Genetics, 11, 1945-1952.
Martinez-Castillero, M., Varona, L., Pegolo, S., Rossoni, A., & Cecchinato, A. (2021). Bayesian inference of the inbreeding load variance for fertility traits in Brown Swiss cattle. Journal of Dairy Science, 104(9), 10040-10048.
Varona, L., Altarriba, J., Moreno, C., Martínez-Castillero, M., & Casellas, J. (2019). A multivariate analysis with direct additive and inbreeding depression load effects. Genetics Selection Evolution, 51(1), 1-12.


