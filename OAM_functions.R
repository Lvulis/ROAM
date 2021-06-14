library(pracma)
library(parallel)
library(compiler)

idExBase <- cmpfun(function(id, nr, nc) {
  # Retrieve x and y indices of an N by M matrix
  # Arguments:
  #   A: Rectangular matrix
  # Returns:
  #   NxM by 2 matrix containing indicse of every point
  idx <- id %% nr
  if(length(idx)==1) {
    if(idx==0) {
      idx <- nr
    }
  } else {
    idx[idx==0] <- nr
  }
  idy <- (id-idx)/nr+1
  return(cbind(idx, idy))
})

idIn <- cmpfun(function(idx, idy, nro, nco) {
  # Turns x and y indices into running index 
  # Arguments:
  #   idx: x coordinates
  #   idy: y coordinates
  #   nro: nrow of matrix
  #   nco: ncol of matrix
  # Returns:
  #   running index equivalents of (idx, idy)
  # function to convert x, y indices to their running index value
  ifelse(((0 < idx) & (idx <= nro)) & ((0 < idy) & (idy <= nco)),
         (idy-1)*nro + idx,
         NA)
})

cRotate <- cmpfun(function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
})

deg2Rad <- cmpfun(function(deg) {
  deg*pi/180
})

openingAngle <- function(testset, queryset, precision = precision, NC = no_cores) {
  # casts Rays and computes opening angle of the queryset against the test set
  # Uses the approximation of seeangles2 from Shaw et al.
  # Arguments:
  #  testset: nx2 matrix containing x, y coordinates of test set pixel centers
  #  queryset: nx2 matrix containing x, y coordinates of query set pixelc enters
  #  precision: integer determing number of viewing angles to use
  # Returns:
  #  Viewing angle of each pixel in the queryset against the testset
  theta = seq(360/precision, 360, length.out = precision)
  rot_array <- array(unlist(lapply(theta, function(x) cRotate(deg2Rad(x)))), dim = c(2, 2, precision))
  r1 <- nrow(queryset)
  hypo <- sqrt(2)/2
  cl <- makeCluster(NC)
  
  # Export objects from local environment:
  clusterExport(cl, varlist=c("r1", "testset", "queryset", "precision",
                              "theta", "hypo", "rot_array"),
                envir=environment())
  # # Export objects from global environment
  clusterExport(cl, varlist=c("deg2Rad","cRotate"))
  
  # rotate all of the differences by theta(j) around the center point
  
  oAset <- parSapply(cl, seq_len(r1), function(i) {
    d <- sweep(testset, 2,  queryset[i, ])
    oa <- apply(rot_array, 3, function(ROT) {
      spin <- ROT %*% t(d)
      examine <- spin[, spin[1, ]>=0, drop = F]
      (min(abs(examine[2, ]))  > hypo)*1L
    })
    sum(oa)
  }) * 360/precision
  
  stopCluster(cl)
  proc.time() - ptm
  return(oAset)
}


Seeangles2 <- function(watermask, precision = 360, no_cores = 3) {
  # Compute opening angle method with modifications for accuracy
  # Primary modification is an iterative expansion of the
  # queryset. Requires pracma library 
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
  # Extract coordinates of water laying on edge (initial set to use)
  edges <- (G > 0) & (watermask > 0)
  Shallowsea_iter <- idExBase(which(edges>0), nrow(edges), ncol(edges))
  rm(G)
  # Get shallow sea pixels ( pixels where threshold image is 1)
  water_ind <- which(watermask>0.5)
  sea <- idExBase(water_ind, nrow(watermask), ncol(watermask))
  K   <- chull(allshore) #convex hull of the edge set
  In <- pracma::inpolygon(sea[, 1], sea[, 2], allshore[K, 1], allshore[K, 2]) # double check that row/column aren't fflippe
  ## Sea outside the convex hull
  ## Turn on if you want to show area outside CHULL but...not necessary

  ## Sea within the convex hull originally from OAM code

  Shallowsea <- matrix(NA, nrow = sum(In), ncol = 3)
  ## The "In" part is basically only to ensure you grow up to 
  ## the convex hull, and isn't a bottleneck.
  ######## INSERT DISPLAY CODE HERE

  ### Compute opening angle of the first set of pixels
  oA <- openingAngle(testset = allshore, queryset = Shallowsea_iter, precision = precision, NC = no_cores)
  
  ## We typically don't care about shorelines extracted at angles
  # greater than 90 degrees, so only compute angles up to that.
  # This can be adjusted (not passed as an argument but could be)
  notSL <- which(oA <= 90)
  N_notSL <- length(notSL)
  iter_bot <- nrow(Shallowsea_iter)
  
  Shallowsea[1:nrow(Shallowsea_iter), 1:2] <- Shallowsea_iter
  Shallowsea[1:nrow(Shallowsea_iter), 3] <- oA
  
  gc()
  bufferind <- which(watermask==0)
  # while N of small viewshed objects is more than or equal to 1
  # dilate the image around ONLY THOSE and look more
  ptm <- proc.time()
  while(N_notSL >= 100) {
    
    # construct empty matrix to grow
    tobuild <- matrix(0, nrow = nrow(watermask), ncol = ncol(watermask))
    tobuild[idIn(Shallowsea_iter[notSL, 1], Shallowsea_iter[notSL, 2], nrow(watermask), ncol(watermask))] <- 1L
    # grow the matrix
    tobuild <- dilate(tobuild, makeBrush(3, 'box'))
    # drop any land pixels
    tobuild[bufferind] <- 0L
    # new possible pixels to get
    Shallowsea_iter <- idExBase(which((tobuild-edges) == 1), nrow(watermask), ncol(watermask))
    edges[tobuild==1] <- 1L
    # compute their oA
    oA <- openingAngle(testset = allshore, queryset = Shallowsea_iter, precision = precision, NC = no_cores)
    
    notSL <- which(oA <= 90)
    N_notSL <- length(notSL)
    
    ### Update running list of preallocated vector. Preallocation is
    ### important to prevent out of control RAM usage
    Shallowsea[(iter_bot)+seq_len(length(oA)), 1:2] <- Shallowsea_iter
    Shallowsea[(iter_bot)+seq_len(length(oA)), 3]   <- oA
    iter_bot <- iter_bot + length(oA)
    ### Self-modify
  
    print(N_notSL)
    gc()
    
  }
  proc.time()- ptm
  
  Shallowsea_short <- Shallowsea[!is.na(Shallowsea[, 1]), ]
  Shallowsea_ind <- idIn(Shallowsea_short[, 1], Shallowsea_short[, 2], nrow(watermask), ncol(watermask))
  
  # Identify which pixels are in the shallowsea but haven't had OA computed
  setabove <- water_ind[-match(Shallowsea_ind, water_ind)]
  # Matrix for final map that'll be output
  OA_map <- matrix(0, nrow = nrow(watermask), ncol = ncol(watermask))
  # Add in OA values for shallowsea pixels that were computed
  OA_map[Shallowsea_ind] <- Shallowsea_short[, 3]
  # Remaining shallowsea pixels are 200
  OA_map[setabove] <- 200L
  
  return(OA_map)
}


