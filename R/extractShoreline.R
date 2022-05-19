extractShoreline <- function(map, theta_c, card = NULL, NA_buff = NULL, pixres = NULL) {

  if(is.null(card)) stop("You need to specify a starting point to walk from")
  if(methods::is(map, "stars")) {
    if(is.null(pixres)) pixres <- dim(map)$x$delta
    map <- as.matrix.stars(map)
  }


  newmap <- map
  newmap[map < th] <- 0
  newmap[newmap>0] <- 1

  newmap <- fillHoles(newmap, 4)

  newmap <- (!fillHoles((!newmap)*1L, 4)) * 1L
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
  COORD <- COORD_REM <- idExBase(which(new_sl>=1), nrow(new_sl), ncol(new_sl))
  colnames(COORD_REM) <- c("X", "Y")

  ### Convert COORDS_REM from pixel number to geographic representation
  # check to make sure these are actually aligned to the real thing and not offset by a pixel
  COORD_REM[, 1] <- (COORD_REM[, 1] - 0.5)*st_dimensions(OA_r)$x$delta + st_dimensions(OA_r)$x$offset
  COORD_REM[, 2] <- (COORD_REM[, 2] - 0.5)*st_dimensions(OA_r)$y$delta + st_dimensions(OA_r)$y$offset
  # COORD_REM <- st_coordinates(s_r_obj2)[, 1:2]
  nr <- nrow(COORD_REM)


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



  # Find starting point as closest COORD to the reference (cx, cy)
  dCOORD <- sqrt(rowSums(sweep(COORD, 2, cbind(cx, cy), "-")^2))
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
      # That part is there to keep going towards the smallest neighb? not clea
      if(length(allmin)==0) {
        allmin<-ind
      }


    }
    print(cnt)

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
    sf::st_linestring(COORD_NEW[i:(i+1), ]) # |>
    # st_sfc()
    #st_sf()
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
  ### REMOVE these goddamn nodes (hope it worked)
  if(length(nodes_rm)>0){
    COORD_NEW <- COORD_NEW[-nodes_rm, ]
  }

  s_r_obj3 <- sf::st_linestring(COORD_NEW) |>
    sf::st_sfc(crs = sf::st_crs(OA_r)) |>
    sf::st_sf() |>
    rename_geometry("geometry")

  s_r_obj3$id <- NA

  shoreline_L[[which(th_l==th)]] <- s_r_obj3

  rm(intersection_mat, new_sl, newmap)
  gc()
}
