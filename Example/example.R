#### Header ####
# Example script to use the OAM functions

library(stars)
library(pracma)
library(parallel)
library(compiler)

source("~/deltas/ROAM/OAM_functions.R")

watermask_r <- read_stars("~/deltas/ROAM/Example/colville.tif", proxy = F)
# watermask_r[[1]][NA_ind] <- 0 # we need NA's for OAM


watermask <- matrix(as.matrix(watermask_r[[1]]), nrow = dim(watermask_r)[1], ncol = dim(watermask_r)[2])

## Example of using OAM
ptm <- proc.time()
OAM_map <- Seeangles2(watermask, precision = 360, no_cores = 3)
proc.time() - ptm
gc()


### Consider writing the OAM map as a backup
OA_r <- st_as_stars(OAM_map)
st_dimensions(OA_r) <- st_dimensions(filled_R)

write_stars(OA_r, "~/deltas/ROAM/example",
            options = c("COMPRESS=LZW", "TFW=NO"),
            type = "UInt16")
###
