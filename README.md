# getPartialInbreeding and TtimesImP.jl

These programs allow computing of partial inbreeding coefficients (Garcia Cortes et al. 2010) and its posterior use in BLUPto estimates individual Inbreeding Depression Load as in Varona et al. 2019; Martinez-Castillero et al. 2021 

The programs need a renumbered (parents before offspring) pedigree with missing parents set to 0. A pedigree example is provided: the genealogy of bull Comet (Garcia-Cortes et al., 2010).



## getPartialInbreeding.f90

This fortran code can be compiled e.g. as `ifort -O3 -heap-arrays -fopenmp getPartialInbreeding.f90 -o getPartialInbreeding`.

And run using e.g.

```
./getPartialInbreeding pedigree_file [first last]
```

the first and last are optional integers and indicate that we want to compute the partial coefficients from animal "first" to animal "last".

For instance

```
$ ./getPartialInbreeding geneal_comet.txt 10 12
 pedfile with          12 lines
 first=           10 last=           12
    Calculating Inbreeding by M&L, elapsed time  4.9998052E-06
 found           21  non zero Partial F
 time for 2nd approach    1.73147599515505  
```

The file `PartialInbreeding.txt` (+first to last if included) is as follows
```
$ cat PartialInbreeding.txt10to12 
 ancestor individual Fpartial numberOfAncestors
        3       10  0.04687500                  9
        4       10  0.03125000                  9
        5       10  0.03125000                  9
        7       10  0.06250000                  9
        1       10  0.01562500                  9
        1       11  0.02343750                 10 
```

so individual 10 has 5 partial inbreeding coefficients due to individuals 1,3,4,5,7. The column `numberOfAncestors` indicates how many ancestors it has (including itself). The order of ancestor within individual is unpredictable as several threads are potentially used, but this does not affect the final result.

## TtimesI-P.jl

This program computes matrix K=T(I-P). The coefficients in K are the ones used to estimated EBVs for inbreeding depression load.   It is called as e.g.

```
$ julia TtimesImP.jl geneal_comet.txt PartialInbreeding.txt1to12 0.01
```

where the last number is a filter so that values in K with absolute value smaller than, in this ecample, 0.01 are not written. This is to avoid many small values that are actually not much useful for estimation (Martinez Castillero et al.). If one wants to see all values, then use a value of 0. This is the output:

```
$ cat K.txt 
       10        3        2        0.0468750000        3        0.0625000000        7
       11        3        2       -0.0117187500        1        0.1250000000        9
       12        6        5       -0.0102539100        1        0.0107421900        3       -0.0312500000        8        0.0312500000        9        0.1250000000       10
```

the first column is the animal, the second column indicates how many animals in column 10 of K get a non-zero value, the third column indicates how many animals have values higher than the thresholds, and then there are pairs of (value, individual) which are the corresponding values in K. E.g. the element (11,9) of K is 0.125.

# References

García-Cortés, L. A., Martí-Avila, J. C., & Toro, M. A. (2010). Fine decomposition of the inbreeding and the coancestry coefficients by using the tabular method. Conservation Genetics, 11, 1945-1952.
Martinez-Castillero, M., Varona, L., Pegolo, S., Rossoni, A., & Cecchinato, A. (2021). Bayesian inference of the inbreeding load variance for fertility traits in Brown Swiss cattle. Journal of Dairy Science, 104(9), 10040-10048.
Varona, L., Altarriba, J., Moreno, C., Martínez-Castillero, M., & Casellas, J. (2019). A multivariate analysis with direct additive and inbreeding depression load effects. Genetics Selection Evolution, 51(1), 1-12.


