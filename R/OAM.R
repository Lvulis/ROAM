#' Compute the Opening Angle Method
#'
#' test 1
#'
#' This function computes the opening angle method based on Shaw et al., 2008 and through an a modification approximation used in Liang et al., (2015a). This approximation uses the rayshader package to efficiently cast rays from an angle phi to every pixel in the query set Q.
#'
#' @param watermask Numeric: binary matrix where 1 indicates water and 0 indicates land. Should be passed with georeferenced info.
#' @param precision Integer: Number of discrete angles to split circle into (angles tested given by 360/precision)
#' @param save_im Logical: Whether or not the output should be written to a geoTIFF via stars
#' @param fn_r Character: Output file name if save_im is true.
#' @param ext st_crs: Geographic extent information from the stars object.
#' @param no_cores Integer: Number of cores to use if running in parallel. Parallel is assumed by default.
#' @param fork Logical: whether to use fork-type parallelization. See details. Only one of fork, socket, and native should be true.
#' @param socket Logical: whether to use socket-type parallelization.
#' @param native Logical: whether to use native-type parallelization.
#'
#' @return Numeric: matrix with same dimensions as watermask
#' @export
#'
#' @examples
#' # See example in Github repo.
OAM <- function(watermask, precision = 360,  save_im = T,
                       fn_r = NULL, ext = NULL, no_cores = 3,
                       fork = F, socket = F, native = F) {
  # first get shoreline and pertinent border
  if((save_im == T) & ((is.null(fn_r)) | (is.null(ext)))) stop("Provide write information")
  partest = sum(c(fork, socket, native))
  if(partest > 1) stop("Choose parallelization scheme")
  if(partest <1) native <- T
  # Arguments:
  #  watermask: Binary watermask (0 = land, 1 = water)
  #  precision: Number of circular sectors to divide the viewshed of a pixel into
  #  no_cores: Number of cores
  # Returns:
  #  Opening Angle map (OA_map)
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

  if(native==T) {
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


  if(fork) {
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
  if(socket) {
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

  if(socket | fork) {
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
