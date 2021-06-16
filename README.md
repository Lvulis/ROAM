# ROAM
R implementation of the Opening Angle Method for delta shoreline extraction.  This is currently a set of functions, but may eventually be re-hashed into an R package.

## Requirements:
Currently requires `pracma` and `EBimage` to handle some computations. Built on `R 4.1.0`

## R implementation and improvements
The Opening Angle Method is a method proposed by Shaw et al. 2008 to define the shoreline of river mouth sections. The method relies on computing how much of the open ocean (i.e. the test set) can be seen from a point in the near shore (i.e. the query set), and then defining a continuous line at a critical viewing angle. John Shaw (personel communication) provided a slightly modified version of the algorithm which speeds up OAM by computing how is blocked by land. Here, we document an additional algorithmic change that speeds up computation time. Please note it before running this script or using this implementation.

The OAM algorithm runtime is proportional to the size of the query set *Q* and the size of the test set *T*. A natural ways to speed up the process involves reducing either the size of *Q* or of *T*. Here we reduce the size of *Q* by starting with an initial query set *Q<sub>o</sub>*, which consists of all water pixels adjacent to land. Then, compute the viewing angle *a* for every *q* in *Q<sub>o</sub>*. For each pixel *q<sub>o</sub>* in *Q<sub>o</sub>* with *a<sub>o</sub>* below some stopping criteria *a<sub>max</sub>*,  add to *Q<sub>o</sub>* the neighboring pixels that aren't already in *Q<sub>o</sub>* to *Q<sub>o</sub>*. Compute *a* for the neigbors and continue to expand the *Q<sub>o</sub>* until very few pixels are left with *a* less than *a<sub>max</max>*.

An additional improvement is parallelization of the computation of `a`. 

I have not run extensive benchmarking to compare this with the OAM version provided by Shaw et al., but expect it to be quite a bit faster. 
