#' Opening Angle Method
#'
#' Compute the opening angle method
#'
#' @details This function computes the opening angle method based on Shaw et al., 2008 and through an a modification approximation used in Liang et al., (2015a). This R based version uses the rayshader package to efficiently cast rays from an angle phi to every pixel in the query set Q.
#'
#' The parallel parameter takes integer values from 0 to 3.  0 will indicate no parallelization. Three parallelization schemes are possible: 1 use rayshaders-native parallelization; 2 use socket-based parallelization, neceessary on windows machines; 3 use fork-based parallelization, possible on *nix machines. Rayshaders native parallelization splits the input image into blocks and computes them in parallel, while the other 2 options will compute chunks of opening angles in parallel.
#'
#'
#' @param watermask Numeric or stars: binary matrix where 1 indicates water and 0 indicates land. The georeferenced image as a stars object can be passed.
#' @param precision Integer: Number of discrete angles to split circle into (angles tested given by 360/precision)
#' @param save_im Logical: Whether or not the output should be written to a geoTIFF via stars
#' @param fn_r Character: Output file name if save_im is true.
#' @param ext st_crs: Geographic extent information from the stars object.
#' @param no_cores Integer: Number of cores to use if running in parallel. Parallel is assumed by default.
#' @param parallel Integer: Parallelization scheme (see details)
#' @return Numeric: matrix with same dimensions as watermask
#' @export
#'
#' @examples
#' # See example in Github repo.
OAM <- function(watermask, precision = 360,  save_im = T,
                fn_r = NULL, ext = NULL, no_cores = 3,
                parallel = 0) {
  # Handle a stars-type input
  if(methods::is(watermask, "stars")) {
    ext <- stars::st_dimensions(watermask)$x$refsys
    watermask <- as.matrix.stars(watermask)
  }
  if((save_im) & ((is.null(fn_r)) | (is.null(ext)))) stop("Provide write information")

  Slopes <- pracma::gradient(watermask)
  G <- sqrt(Slopes$X^2 + Slopes$Y^2)
  rm(Slopes)
  # extract coordinates of edge pixels (edge set to test against)
  allshore <- idExBase(which((G > 0) & (watermask == 0)),
                       nrow(watermask),
                       ncol(watermask))
  # Extract coordinates of water laying on edge (initial query set)
  edges <- (G > 0) & (watermask > 0)
  rm(G)
  # Get shallow sea pixels ( pixels where threshold image is 1)
  water_ind <- which(watermask>0.5)
  sea <- idExBase(water_ind, nrow(watermask), ncol(watermask))
  K   <- grDevices::chull(allshore) #convex hull of the edge set

  # Test what points are *in* the convex hull. This is used to set up
  # the maximum possible query set size
  In <- pracma::inpolygon(sea[, 1], sea[, 2], allshore[K, 1], allshore[K, 2]) # double check that row/column aren't fflippe
  tm = sweepimage = matrix(0, nrow = nrow(edges), ncol = ncol(edges))
  tm[water_ind[In]] <- 1

  if(parallel <= 1) {
    if(parallel == 0) {no_cores = 1}
    print("Using rayshader's image breakdown parallelization")
    for(k in seq(360/precision, 360, length.out = precision)) {
      print(k)
      sweepimage = sweepimage + floor(rayshader::ray_shade((!watermask)*1L, anglebreaks = seq(1e-4, 1e-3, by = 1e-4), sunangle = k,
                                                zscale = .1, lambert = F, cachemask = tm,
                                                multicore=T, no_cores = no_cores))
      gc()
    }
    sweepimage2 <- sweepimage[nrow(sweepimage):1, ]*360/precision
  }

  nb = precision/no_cores

  if(parallel == 3) {
    print("Parallelizing angles into chunks using fork (optimal on *nix)")
    sweepimage_par =  parallel::mclapply(0:(no_cores-1), \(j) {
      print(j)
      for(k in (j*nb+(1:nb))*360/precision) {
        sweepimage = sweepimage + floor(rayshader::ray_shade((!watermask)*1L, anglebreaks = seq(1e-4, 1e-3, by = 1e-4), sunangle = k,
                                                  zscale = .1, lambert = F, cachemask = tm,
                                                  multicore=F, no_cores = no_cores))
        gc()
      }
      gc()
      return(sweepimage)
    })
  }
  if(parallel == 2) {
    print("Parallelizing angles into chunks using socket")
    cl <- parallel::makeCluster(no_cores)

    # Export objects from local environment:
    parallel::clusterExport(cl, varlist=c("nb", "precision",
                                          "watermask", "no_cores", "sweepimage"),
                            envir=environment())
    parallel::clusterEvalQ(cl, library(rayshader))

    sweepimage_par = parallel::parLapply(cl, 0:(no_cores-1), \(j) {
      print(j)
      for(k in (j*nb+(1:nb))*360/precision) {
        print(k)
        sweepimage = sweepimage + floor(rayshader::ray_shade((!watermask)*1L, anglebreaks = seq(1e-4, 1e-3, by = 1e-4), sunangle = k,
                                                  zscale = .1, lambert = F, cachemask = tm,
                                                  multicore=F, no_cores = no_cores))
        gc()
      }
      gc()
      return(sweepimage)
    })
    parallel::stopCluster(cl)

  }

  if(parallel >= 2) {
    sweepimage2 = matrix(0, nrow = nrow(sweepimage), ncol = ncol(sweepimage))
    for (j in 1:(no_cores)) {
      sweepimage2 = sweepimage2 + sweepimage_par[[j]]
    }
    sweepimage2 <- sweepimage2[nrow(sweepimage2):1, ]*360/precision

  }

  sweepimage2[watermask==0] <- 0

  sweepimage2[water_ind[-match(water_ind[In], water_ind)]] <- 200L

  if(save_im == T) {
    OA_r <- stars::st_as_stars(sweepimage2)
    stars::st_dimensions(OA_r) <- ext

    stars::write_stars(OA_r, fn_r,
                       options = c("COMPRESS=LZW", "TFW=NO"))

  }

  return(sweepimage2)
}
