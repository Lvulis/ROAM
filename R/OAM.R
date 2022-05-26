#' Opening Angle Method
#'
#' Compute the opening angle method
#'
#' @details This function implements the Opening Angle Method. This R based version uses an approximation of the OAM proposed by John Shaw and commonly used by the delta morphology community. See the `ROAM` package vignette for details.
#'
#' The parallel parameter takes integer values from 0 to 3.  0 indicates no parallelization. Three parallelization schemes are possible: 1 use rayshaders-native parallelization which dissects an input scene into chunks before raytracing for a given angle; 2 and 3 split up the computation of the raytracing into `no_cores` blocks and then tallies the results. `2` indicates the use socket-based parallelization, necessary on windows machines; 3 use fork-based parallelization, possible on *nix machines.
#'
#' @param watermask Numeric or stars: binary n x m matrix where 1 indicates water and 0 indicates land. The input image should ideally be passed as a georeferenced stars object.
#' @param precision Integer: Number (`p`) of discrete sectors to split circle into (angles tested given by 360/precision)
#' @param save_im Logical: Whether or not the output should be written to a geoTIFF via `stars``
#' @param fn_r Character: Output file name if save_im is true.
#' @param geoInfo st_dimensions: Georaphic information including coordinate reference system, resolution, etc. of the georeferenced output. If a stars object is passed to `map`, its st_dimensions properties will be used.
#' @param no_cores Integer: Number of cores to use if running in parallel. No parallelization is assumed by default.
#' @param parallel Integer: Parallelization scheme (see details)
#'
#' @return Numeric: matrix with same dimensions as watermask
#' @export
#'
#' @examples
#' \dontrun{
#' watermask = stars::read_stars("SaoFrancisco.tif")
#' plot(watermask, main = "Sao Francisco")
#' OAM_map = OAM(watermask, precision = 720, save_im = T,
#' fn_r = "SaoFrancisco_OAM.tif", geoInfo = NULL, no_cores = 3,
#' parallel = 2)
#' }
OAM <- function(watermask, precision = 360,  save_im = T,
                fn_r = NULL, geoInfo = NULL, no_cores = 3,
                parallel = 0) {
  # Handle a stars-type input
  if(methods::is(watermask, "stars")) {
    geoInfo <- stars::st_dimensions(watermask)
    watermask <- as.matrix.stars(watermask)
  }
  if((save_im) & ((is.null(fn_r)) | (is.null(geoInfo)))) stop("Provide write information")

  Slopes <- pracma::gradient(watermask)
  G <- sqrt(Slopes$X^2 + Slopes$Y^2)
  rm(Slopes)
  # extract coordinates of edge pixels (edge set to test against)
  allshore <- id_ex_base(which((G > 0) & (watermask == 0)),
                       nrow(watermask),
                       ncol(watermask))
  # Extract coordinates of water laying on edge (initial query set)
  edges <- (G > 0) & (watermask > 0)
  rm(G)
  # Get shallow sea pixels ( pixels where threshold image is 1)
  water_ind <- which(watermask>0.5)
  sea <- id_ex_base(water_ind, nrow(watermask), ncol(watermask))
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
    stars::st_dimensions(OA_r) <- geoInfo

    stars::write_stars(OA_r, fn_r,
                       options = c("COMPRESS=LZW", "TFW=NO"))

  }

  return(sweepimage2)
}
