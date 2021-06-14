# ROAM
R implementation of the Opening Angle Method with some improvements. May eventually be re-hashed into an R package, but available as a set of functions & an example for now.

# Requirements:
Currently requires `pracma` to handle some computations. Built on `R 4.1.0`

# R implementation and improvements
The Opening Angle Method is a method proposed by Shaw et al. 2008 to define the shoreline of river mouth sections. The method relies on computing how much of the open ocean (i.e. the test set) can be seen from a point in the near shore (i.e. the query set), and then defining a continuous line at a critical viewing angle. John Shaw (personel communication) provided a slightly modified version of the algorithm which speeds up OAM by computing how is *blocked* by land. Here, we document an additional algorithmic change that speeds up computation time.

The OAM algorithm runtime is proportional to the size of the query set `Q` and the size of the test set `T`. A natural ways to speed up the process involves reducing either the size of `Q` or of `T`. Here we reduce the size of `Q` by starting with an initial query set `Q<sub>o</sub>`, which consists of all water pixels adjacent to land. Then, compute the viewing angle `&theta` for every `q` in `Q<sub>o</sub>`. For each pixel `q<sub>o</sub>` in `Q<sub>o</sub>` with &theta<sub>o</sub> below some stopping criteria `&theta[max]`, then add to `Q<sub>o</sub>` add the neighboring pixels that aren't in `Q<sub>o</sub>` to `Q<sub>o</sub>`. Iteratively compute `&theta` for the neighbors until most pixels have been evaluated.




