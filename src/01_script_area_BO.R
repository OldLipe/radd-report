library(sits)
library(bench)
library(bayts)
library(raster)
library(parallel)
library(bfast)
library(ggbeeswarm)

#---- Execution of RADD-SITS ----#
#
# Creating cube using Landsat-7 images
#
cube_l7 <- sits_cube(
    source = "MPC",
    collection = "LANDSAT-C2-L2",
    data_dir = "./data/raw/images_radd_exp_l7/",
    parse_info = c("X1", "band", "tile", "date"),
    start_date = "2015-02-08",
    end_date = "2016-05-25"
)

#
# Creating cube using Sentinel-1 images
#
cube_s1 <- sits_cube(
    source = "MPC",
    collection = "SENTINEL-1-RTC",
    data_dir = "./data/raw/images_radd_exp_s1/",
    parse_info = c("X1", "band", "tile", "date"),
    start_date = "2015-02-08",
    end_date = "2016-05-25"
)

#
# Defining mean values for NDVI and VV
#
mean_stats <- tibble::tibble(
    "label" = c("F", "NF"),
    "NDVI" = c(0, -0.5),
    "VV" = c(-1, -4)
)

#
# Defining standard deviation values for NDVI and VV
#
sd_stats <- tibble::tibble(
    "label" = c("F", "NF"),
    "NDVI" = c(0.1, 0.125),
    "VV" = c(0.75, 1)
)

#
# Combining Sentinel-1 and Landsat-7 cubes
#
cube <- sits_merge(cube_l7, cube_s1, irregular = TRUE)

#
# Defining the parameters of RADD execution
#
res_r <- sits:::sits_radd(
    data = cube,
    mean_stats = mean_stats,
    sd_stats = sd_stats,
    deseasonlize = NULL,
    multicores = 1,
    memsize = 4,
    start_date = "2016-01-02",
    end_date = "2016-05-26",
    impute_fn = identity,
    output_dir = "./data/output/",
    version = "RADD-SITS-V1"
)

#
# Compute evaluation metrics for RADD-SITS
#
predict_sits <- bench::mark(
    `RADD-SITS` = res_r(),
    iterations = 10,
    time_unit = "s",
    memory = TRUE
)

#---- Execution of RADD-REICHE ----#
#
# Reading Landsat-7 and Sentinel-1 images as a brick
#
lndviD <- raster::stack("./data/raw/LANDSAT_NDVI_230073_BRICK.tif")
s1vvD <- raster::stack("./data/raw/SENTINEL_VV_230073_BRICK.tif")

#
# Reading the image dates
#
lndvi_date <- readRDS("./data/raw/LNDVI_DATES.rds")
s1vv_date <- readRDS("./data/raw/S1VV_DATES.rds")

#
# Defining mean and standard deviation for each sensor
#
s1vvD_pdf <- c(c("gaussian","gaussian"),c(-1,0.75),c(-4,1))
lndviD_pdf <- c(c("gaussian","gaussian"),c(0,0.1),c(-0.5,0.125))

#
# Defining the parameters of RADD execution
#
chi = 0.9
start = 2016

#
# Running RADD-REICHE implementation
#
res <- baytsSpatial(
    bL       = list(lndviD, s1vvD),
    datesL   = list(lndvi_date, lndvi_date),
    pdfL     = list(lndviD_pdf, s1vvD_pdf),
    chi      = chi,
    start    = start,
    mc.cores = 10
)

#
# Function to convert zoo dates to yday
#
convert_date_to_yday <- function(r, ref_date) {
    v <- raster::values(r)

    yday365 <- function(x) {
        x <- as.POSIXlt(x)
        mdays_sum <- c(0L, 31L, 59L, 90L, 120L, 151L, 181L,
                       212L, 243L, 273L, 304L, 334L, 365)
        mdays_sum[1L + x$mon] + x$mday
    }

    dates_transformed <- 1900 + as.POSIXlt(ref_date)$year + (yday365(ref_date) - 1) / 365
    names(ref_date) <- dates_transformed
    v <- lubridate::date(unname(ref_date[as.character(v)]))
    v <- as.numeric(paste0(lubridate::year(v), lubridate::yday(v)))
    raster::values(r) <- v
    r
}

new_rast <- convert_date_to_yday(res[[3]], lndvi_date)

raster::writeRaster(
    x = new_rast,
    datatype = "INT4U",
    filename = "./data/output/LANDSAT_TM-ETM-OLI_230073_2015-02-08_2016-05-25_radd_radd-reiche-v1.tif",
    options = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES",
                "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512")
)

#
# Compute evaluation metrics for RADD-REICHE (~ 30min)
#
predict_radd <- bench::mark(
    `RADD-REICHE` = baytsSpatial(
        bL       = list(lndviD, s1vvD),
        datesL   = list(lndvi_date, lndvi_date),
        pdfL     = list(lndviD_pdf, s1vvD_pdf),
        chi      = chi,
        start    = start,
        mc.cores = 1
    ),
    iterations = 10,
    time_unit = "s",
    memory = TRUE
)

#
# Combining results in a unique tibble
#
res <- dplyr::bind_rows(predict_sits, predict_radd)

#
# Saving results as RDS (large file)
#
saveRDS(res, "./data/output/results_bench.rds")

#
# Plotting the execution time
#
ggplot2::autoplot(res)

#
# Saving output plot
#
ggplot2::ggsave(
    "./data/output/results_comparing.png",
    width = 1600,
    height = 800,
    units = "px",
    device = "png"
)

#---- Computing the difference map ----#

#
# Reading results maps with terra package
#
radd_sits <- terra::rast("./data/output/LANDSAT_TM-ETM-OLI_230073_2015-02-08_2016-05-25_radd_radd-sits-v1.tif")
radd_reiche <- terra::rast("./data/output/LANDSAT_TM-ETM-OLI_230073_2015-02-08_2016-05-25_radd_radd-reiche-v1.tif")

#
# Computing the difference between them
#
radd_diff <- radd_sits - radd_reiche

#
# Checking for differences
#
radd_diff_freq <- terra::freq(radd_diff)

#
# Expecting that all pixels are equal 0
#
testthat::expect_equal(
    radd_diff_freq["count"], 1542
)
testthat::expect_equal(
    radd_diff_freq["valor"], 0
)

#
# Writing the difference map with raster package
#
terra::writeRaster(
    x = radd_diff,
    datatype = "INT1U",
    filename = "./data/output/DIFFERENCE_MAP.tif",
    gdal = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES",
             "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512")
)

#
# Checking for checksum
#
values_sits <- unname(terra::values(radd_sits))
values_reiche <- unname(terra::values(radd_reiche))

#
# Checking values
#
testthat::expect_equal(
    digest::digest(values_sits, algo = "md5"),
    digest::digest(values_reiche, algo = "md5")
)
