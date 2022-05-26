#' Extract Shoreline vector
#'
#' Extract a shoreline vector
#'
#' This function is used to define a continuous shoreline vector (i.e. linestring geometry) separating the subaerial delta from the subaqueous delta & the ocean.
#' By definition, a continuous shoreline vector has an ordered set of coordinates, which is useful for conducting certain analyses, e.g. the spectral analyses described in Vulis et al., ($DATE). A shoreline vector is also necessary for the the Python package RivGraph which is used to analyze delta channel networks.
#' As the shoreline is defined at the pixel scale, it is often rough and ragged. This causes the typical vectorization (polygonization) approaches implemented in the `stars` and `sf` space to fail. To overcome this, a "walker" approach, wherein the pixels are sorted by walking along the shoreline.
#' The `card` argument specifies what position along the edge of the image to start walking from, and takes on values from 1 to 16 which correspond to the following locations:
#' ## NW = 1, N = 2.... NW-W = 16
#1  0  0 2  0 0  3 0 0  4 0 5
#0  0  0 0  0 0  0 0 0  0 0 0
#0  0  0 0  0 0  0 0 0  0 0 0
#16 0  0 0  0 0  0 0 0  0 0 6
#0  0  0 0  0 0  0 0 0  0 0 0
#0  0  0 0  0 0  0 0 0  0 0 0
#15 0  0 0  0 0  0 0 0  0 0 7
#0  0  0 0  0 0  0 0 0  0 0 0
#0  0  0 0  0 0  0 0 0  0 0 0
#14 0  0 0  0 0  0 0 0  0 0 8
#0  0  0 0  0 0  0 0 0  0 0 0
#0  0  0 0  0 0  0 0 0  0 0 0
#13 0  12 0 0 11 0 0 10 0 0 9
#' It is highly recommended to pass a georeferenced input image to `map` and avoid any possible issues.
#'
#' @param map numeric or stars: OAM map to be processed
#' @param theta numeric: critical value at which to define a shoreline
#' @param card integer: Cardinal direction at which to start walking (see details)
#' @param NA_buff integer: indices to remove before walking (see details)
#' @param pixres numeric: pixel resolution. Shouldn't be passed if passing georeferenced image
#' @param outCRS st_crs: coordinate reference system of output geometry. Shouldn't be passed if passing georeferenced image
#' @param xoff numeric: offset value to add to x coordinates. Shouldn't be passed if passing georeferenced image
#' @param yoff numeric: offset value to add to y coordinates. Shouldn't be passed if passing georeferenced image
#'
#' @return `sf` object storing the shoreline vector.
#' @export
#'
#' @examples
#' \dontrun{
#' fn_OAM = system.file("extdata", "SaoFrancisco_OAM.tif", package = "ROAM")
#' OAM_map = stars::read_stars(fn_OAM)
#' shoreline = extract_shoreline(OAM_map, theta = 45, card = 13, NA_buff = NULL,
#' pixres= NULL, outCRS = NULL)
#' }
extract_shoreline <- function(map, theta = 45, card = NULL,
                              NA_buff = NULL, pixres = NULL,
                              outCRS = NULL, xoff = NULL, yoff = NULL) {
  # Error handling
  if(is.null(card)) stop("You need to specify a starting point to walk from")
  if(methods::is(map, "stars")) {
    if(is.null(pixres)) pixres <- stars::st_dimensions(map)$x$delta
    if(is.null(outCRS)) outCRS <- sf::st_crs(map)
    if(is.null(xoff))   xoff   <- stars::st_dimensions(map)$x$offset
    if(is.null(yoff))   yoff   <- stars::st_dimensions(map)$y$offset
    deltax = stars::st_dimensions(map)$x$delta
    deltay = stars::st_dimensions(map)$y$delta
    map <- as.matrix.stars(map)
  }
  if(is.null(pixres)) stop("You need to specify an input pixel resolution")
  if(is.null(outCRS)) stop("You need to specify an output coordinate reference system")
  if(is.null(xoff))   stop("You need to specify an input x-offset")
  if(is.null(yoff))   stop("You need to specify an input y-offset")
  if(!exists("deltax")) deltax = deltay = pixres
  ## CHECK THIS
  newmap <- bin_thresh(map, theta)


  newmap <- fill_holes(newmap, 4)

  newmap <- (!fill_holes((!newmap)*1L, 4)) * 1L
  mode(newmap) <- 'integer'
  kern <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)

  ## Identify boundary pixels
  # new_sl <- newmap - EBImage::erode(newmap, kern)
  new_sl <- newmap - round(as.matrix(imager::erode(imager::as.cimg(newmap), mask = imager::as.cimg(kern))))

  if(!is.null(NA_buff)) {
    new_sl[NA_buff==1] <- 0
  }

  new_sl <- keep_largest(new_sl)

  ####WALKER search each time...
  COORD <- COORD_REM <- id_ex_base(which(new_sl>=1), nrow(new_sl), ncol(new_sl))
  colnames(COORD_REM) <- c("X", "Y")

  ### Convert COORDS_REM from pixel number to geographic representation
  # check to make sure these are actually aligned to the real thing and not offset by a pixel
  COORD_REM[, 1] <- (COORD_REM[, 1] - 0.5)*deltax + xoff
  COORD_REM[, 2] <- (COORD_REM[, 2] - 0.5)*deltay + yoff
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
    tf = (suppressWarnings(min(abs(d_v))) >= (30*pixres))
    if(tf) {
      # print("TOOBIG")
      cnt <- nr
    } else if(!tf) {
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
    sf::st_intersection(temp_LINE_package[idmin, ], temp_LINE_package[ii, ])
    # Keep the earleist line, and shortcut from its endpoint to the start
    # of the next line. This'll clip out what was there.
    noderm <- setdiff(idmin:idmax, c(idmin, idmin+1, idmax))
  }))
  ### REMOVE these nodes
  if(length(nodes_rm)>0){
    COORD_NEW <- COORD_NEW[-nodes_rm, ]
  }
  s_r_obj3 <- sf::st_linestring(COORD_NEW) |>
    sf::st_sfc(crs = outCRS) |>
    sf::st_sf() |>
    rename_geometry("geometry")
  s_r_obj3$id <- NA

  return(s_r_obj3)
}
