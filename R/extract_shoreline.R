#' Extract Shoreline vector
#'
#' Extract a shoreline vector
#'
#' This function is used to define a continuous shoreline vector (i.e. linestring geometry) separating the subaerial delta from the subaqueous delta & the ocean.
#' By definition, a continuous shoreline vector has an ordered set of coordinates, which is useful for conducting certain analyses, e.g. the spectral analyses described in Vulis et al., ($DATE). A shoreline vector is also necessary for the the Python package RivGraph which is used to analyze delta channel networks.
#' As the shoreline is defined at the pixel scale, it is often rough and ragged. This causes the typical vectorization (polygonization) approaches implemented in the `stars` and `sf` space to fail. To overcome this, a "walker" approach, wherein the pixels are sorted by walking along the shoreline.
#'
#' At
#' @param map numeric: n x m matrix corresponding to an OAM image
#' @param theta numeric: critical value at which to define a shoreline
#' @param card integer: Cardinal direction at which to start walking (see details)
#' @param NA_buff integer: indices to remove before walking (see details)
#' @param pixres numeric: Pixel resolution (see details)
#'
#' @return `sf` object storing the shoreline vector.
#' @export
#'
#' @examples
extract_shoreline <- function(map, theta = 45, card = NULL, NA_buff = NULL, pixres = NULL) {
  # Error handling
  if(is.null(card)) stop("You need to specify a starting point to walk from")
  if(methods::is(map, "stars")) {
    if(is.null(pixres)) pixres <- dim(map)$x$delta
    map <- as.matrix.stars(map)
  }
  if(is.null(pixres)) stop("You need to specify an input pixel resolution")

  ## CHECK THIS
  newmap <- bin_thresh(theta)
  # newmap <- map
  # newmap[map < th] <- 0
  # newmap[newmap>0] <- 1

  newmap <- fill_holes(newmap, 4)

  newmap <- (!fill_holes((!newmap)*1L, 4)) * 1L
  mode(newmap) <- 'integer'
  kern <- makeBrush(3, 'diamond')
  kern[2, 2] <- 0

  ## Identify boundary pixels
  new_sl <- newmap - EBImage::erode(newmap, kern)
  if(!is.null(NA_buff)) {
    new_sl[NA_buff==1] <- 0
  }

  new_sl <- keep_largest(new_sl)
  new_sl <- bin_thresh(new_sl, 1)

  ####WALKER search each time...
  COORD <- COORD_REM <- id_ex_base(which(new_sl>=1), nrow(new_sl), ncol(new_sl))
  colnames(COORD_REM) <- c("X", "Y")

  ### Convert COORDS_REM from pixel number to geographic representation
  # check to make sure these are actually aligned to the real thing and not offset by a pixel
  COORD_REM[, 1] <- (COORD_REM[, 1] - 0.5)*stars::st_dimensions(OA_r)$x$delta + stars::st_dimensions(OA_r)$x$offset
  COORD_REM[, 2] <- (COORD_REM[, 2] - 0.5)*stars::st_dimensions(OA_r)$y$delta + stars::st_dimensions(OA_r)$y$offset
  # COORD_REM <- st_coordinates(s_r_obj2)[, 1:2]
  nr <- nrow(COORD_REM)

  # Find starting point as closest COORD to the reference starter = (cx, cy)
  starter = card_select(newmap, card)

  dCOORD <- sqrt(rowSums(sweep(COORD, 2, starter, "-")^2))
  allmin <- ind <- which.min(dCOORD)
  gc()

  # COORD_new: Preallocate matrix to store the sorted coordinates
  cnt <-1
  COORD_NEW <- matrix(0, nrow = nrow(COORD_REM), ncol = ncol(COORD_REM))
  COORD_NEW[1, ] <- COORD_REM[ind, ] # First item is our start
  #
  while(cnt < nr) {
    # Reference coord is the i_ind value
    i_ind <- COORD_REM[ind, , drop = F]
    COORD_REM <- COORD_REM[-allmin, , drop = F] #Drop "allmin"

    # Smallest neighbor L2 distance from remainder of points.
    d_v <- sqrt(rowSums(sweep(COORD_REM, 2, i_ind, "-")^2))

    # When the distances are too far- you've walked off the trailand/or done!
    if(min(abs(d_v)) >= (30*pixres)) {
      print("TOOBIG")
      cnt <- nr
    } else if(min(abs(d_v)) < (30*pixres)) {
      ## INTRODUCE A "NUKE" -->
      # If there are multiple objects with distance < sqrt(2),
      # then remove from COORD_REM all the guys that aren't ind
      # This isn't a very good fix, the best would be to keep track and if
      # you go too far down the wrong path then turn around,
      # but that could involve a lot of memory. This is a very basic
      # way to keep searching the space
      ind <- which.min(abs(d_v)) # identify closest neighbor

      cnt <- cnt+1
      COORD_NEW[cnt, ] <- COORD_REM[ind, ] # count up and add a number
      allmin <- which(abs(d_v)<=(pixres*sqrt(2)*1.01)) # Identify next nearest neighbor? that's at most a diag-pix away?
      # That part is there to keep going towards the smallest neighb?
      if(length(allmin)==0) {
        allmin<-ind
      }
    }
  }
  # if any COORD_New are at 0,0 drop 'em
  if(any(rowSums(COORD_NEW)<1e-10)) {
    COORD_NEW <- COORD_NEW[-which(rowSums(COORD_NEW)<1e-10), ]
  }


  ### Loops can pop up if there's possibly an extra pixel
  ### Remove loops of any size by trimming all the pixels in between
  ### the looped index


  # First create lines from every consecutive pair of points
  temp_LINE_package <- lapply((1:(nrow(COORD_NEW)-1)), \(i) {
    sf::st_linestring(COORD_NEW[i:(i+1), ])
  }) |>
    sf::st_sfc() |>
    sf::st_sf()

  # Identify where all the lines intersect:

  intersection_mat <- sf::st_intersects(temp_LINE_package)
  # Some lines may intersect MULTIPLE times -- suggestive of a
  # loop. We test if the 2nd element in the intersection is the self-intersection
  intersect_test <- sapply(seq_along(intersection_mat[-1]), \(i) {
    intersection_mat[[i+1]][2]==i+1
  })
  toclip <- which(intersect_test==F)+1

  nodes_rm <- unlist(lapply(toclip, \(ii) {
    # first, find the smallest link #
    idmin <- min(intersection_mat[[ii]])
    # second, find the largest link
    # (will basically route in between these)
    idmax <- max(intersection_mat[[ii]])
    # Interesection point:
    st_intersection(temp_LINE_package[idmin, ], temp_LINE_package[ii, ])
    # Keep the earleist line, and shortcut from its endpoint to the start
    # of the next line. This'll clip out what was there.
    noderm <- setdiff(idmin:idmax, c(idmin, idmin+1, idmax))
  }))
  ### REMOVE these nodes
  if(length(nodes_rm)>0){
    COORD_NEW <- COORD_NEW[-nodes_rm, ]
  }
  s_r_obj3 <- sf::st_linestring(COORD_NEW) |>
    sf::st_sfc(crs = sf::st_crs(OA_r)) |>
    sf::st_sf() |>
    rename_geometry("geometry")
  s_r_obj3$id <- NA

  return(s_r_obj3)
}
