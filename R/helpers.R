id_ex_base <- function(id, nr, nc) {
  # Conv x and y indices of an N by M matrix
  # Arguments:
  #   A: Rectangular matrix
  # Returns:
  #   NxM by 2 matrix containing indices of every point
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
}

as.matrix.stars <- function(x, i) {
  # Convert the first layer of a stars object into a matrix
  # Necessary for some kind of javascript error that arises sometimes
  # Arguments:
  #  x: stars object on a regular grid
  #  i: band
  # Returns:
  #  matrix containing the gridded values
  return(matrix(as.matrix(x[[1]][, , i]), nrow = dim(x)[1], ncol = dim(x)[2]))
}

fill_holes <- function(scene, maxholesize = 0) {
  # Translated from Jon Schwenk, RivGraph.
  # Fills holes in binary mask of size <= maxholesize
  if (maxholesize == 0) {
    scene = EBImage::fillHull(scene)
  } else {
    scenecomp <- (!scene)*1L
    mode(scenecomp) <- 'integer'
    scenecomp <- EBImage::bwlabel(scenecomp)
    mode(scenecomp) <- 'integer'
    ox <- lengths(split_objects(scenecomp))

    scenecomp <- EBImage::rmObjects(scenecomp, names(which(ox <= maxholesize)), reenumerate = F)
    scenecomp[scenecomp > 0] <- 1L
    scene <- (!scenecomp)*1L
  }
  scene
}

split_objects = function(x) {
  # From EBImage: get R indices of each object in x matrix
  # Arguments:
  #  x: Integer valued labelled image
  # Returns:
  #  list w/ each entry corresponding to indices of a unique label in x
  z = which(as.integer(x)>0L)
  split(z, x[z])
}

keep_largest <- function(img) {
  # Keep only largest element of an image. Relies on python
  # Arguments:
  #  img: integer-valued binary matrix
  # Returns:
  #  Same sized matrix only containing largest cluster (label id still there)
  cncmp <- round(as.matrix(imager::label(imager::as.cimg(img), high_connectivity = T)))

  mode(cncmp) <- 'integer'
  ox <- lengths(split_objects(cncmp))

  if(length(ox) > 1) {
    tokp <- names(which.max(ox))
    clean_mask <- EBImage::rmObjects(cncmp, setdiff(names(ox), tokp), reenumerate = F)
    return(clean_mask)
  } else {
    return(cncmp)
  }

}

bin_thresh <- function(x, thresh) {
  # Binary threshold a numeric/integer vector or matrix
  # Arguments:
  #  x: Matrix
  #  thresh: threshold value
  # Returns:
  #  x: Binary matrix
  x[x < thresh] <- 0L
  x[x >= thresh] <- 1L
  mode(x) <- 'integer'
  x
}

rename_geometry <- function(g, name){
  # Rename the geometry column of sf `g` into `name``
  # Arguments:
  #  g: sf object
  #  name: new geometry column name (string)
  # Returns:
  #  x: sf object
  current = attr(g, "sf_column")
  names(g)[names(g)==current] = name
  st_geometry(g)=name
  g
}

card_select <- function(newmap, card) {
  # Obtain the coordinate of in the "card" position along edge of a matrix
  # Arguments:
  #  newmap: n x m matrix
  #  card: cardinal location, see guideline below
  # Returns:
  #  2-vector with cx, cy.
  if(card == 1) {
    cx = 1; cy = 1
  } else if(card == 2) {
    cx = round(nrow(newmap)/4); cy = 1
  } else if(card == 3) {
    cx = round(nrow(newmap)/2); cy = 1
  } else if(card == 4) {
    cx = round(nrow(newmap)*3/4); cy = 1
  } else if(card == 5) {
    cx = nrow(newmap); cy = 1
  } else if(card == 6) {
    cx = nrow(newmap); cy = round(ncol(newmap)/4)
  } else if(card == 7) {
    cx = nrow(newmap); cy = round(ncol(newmap)/2)
  } else if(card == 8) {
    cx = nrow(newmap); cy = round(ncol(newmap)*3/4)
  } else if(card == 9) {
    cx = nrow(newmap); cy = ncol(newmap)
  } else if(card == 10) {
    cx = round(nrow(newmap)*3/4); cy = ncol(newmap)
  } else if(card == 11) {
    cx = round(nrow(newmap)/2); cy = ncol(newmap)
  } else if(card == 12) {
    cx = round(nrow(newmap)/4); cy = ncol(newmap)
  } else if(card == 13) {
    cx = 1; cy = ncol(newmap)
  } else if(card == 14) {
    cx = 1; cy = round(ncol(newmap)/4)
  } else if(card == 15) {
    cx = 1; cy = round(ncol(newmap)/2)
  } else if(card == 16) {
    cx = 1; cy = round(ncol(newmap)*3/4)
  }
  matrix(c(cx, cy), ncol = 2)
}


