library(sits)
library(bench)
library(bayts)

#devtools::install_github("oldlipe/bayts@feat/experiments")

cube_l7 <- sits_cube(
    source = "MPC",
    collection = "LANDSAT-C2-L2",
    data_dir = "~/radd_report/data/raw/images_radd_exp_l7/",
    parse_info = c("X1", "band", "tile", "date"),
    start_date = "2015-02-08",
    end_date = "2016-05-25"
)

cube_s1 <- sits_cube(
    source = "MPC",
    collection = "SENTINEL-1-RTC",
    data_dir = "~/radd_report/data/raw/images_radd_exp_s1/",
    parse_info = c("X1", "band", "tile", "date"),
    start_date = "2015-02-08",
    end_date = "2016-05-25"
)

mean_stats <- tibble::tibble(
    "label" = c("F", "NF"),
    "NDVI" = c(0, -0.5),
    "VV" = c(-1, -4)
)

sd_stats <- tibble::tibble(
    "label" = c("F", "NF"),
    "NDVI" = c(0.1, 0.125),
    "VV" = c(0.75, 1)
)

cube <- sits_merge(cube_l7, cube_s1, irregular = TRUE)

res_r <- sits_radd(
    data = cube,
    mean_stats = mean_stats,
    sd_stats = sd_stats,
    deseasonlize = NULL,
    multicores = 1,
    memsize = 4,
    start_date = "2016-01-02",
    end_date = "2016-05-26",
    impute_fn = identity,
    output_dir = "/home/sits/",
    version = "RADD-SITS-V1"
)

predict_sits <- bench::mark(
    `RADD-SITS` = res_r(),
    iterations = 10,
    time_unit = "s",
    memory = TRUE

)


lndviD <- raster::stack("~/radd_report/LANDSAT_NDVI_230073_BRICK.tif")
s1vvD <- raster::stack("~/radd_report/SENTINEL_VV_230073_BRICK.tif")

lndvi_date <- readRDS("~/radd_report/LNDVI_DATES.rds")
s1vv_date <- readRDS("~/radd_report/S1VV_DATES.rds")

s1vvD_pdf <- c(c("gaussian","gaussian"),c(-1,0.75),c(-4,1))
lndviD_pdf <- c(c("gaussian","gaussian"),c(0,0.1),c(-0.5,0.125))

chi = 0.9
start = 2016

predict_sits <- bench::mark(
    `RADD-SITS` = res_r(),
    iterations = 10,
    time_unit = "s",
    memory = TRUE

)

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

res <- dplyr::bind_rows(predict_sits, predict_radd)
res$result <- NULL

saveRDS(res, "~/radd_report/results_bench.rds")

ggplot2::autoplot(res)

res |>
    tidyr::unnest(c(time, gc)) |>
    dplyr::filter(gc == "none") |>
    dplyr::mutate(expression = as.character(expression)) |>
    ggplot2::ggplot(ggplot2::aes(x = mem_alloc, y = time, color = expression)) +
    ggplot2::geom_point() +
    scale_color_bench_expr(scales::brewer_pal(type = "qual", palette = 3))

radd_sits <- terra::rast("~/radd_report/LANDSAT_TM-ETM-OLI_230073_2015-02-08_2016-05-25_radd_v1.tif")
radd_reiche <- terra::rast("~/radd_report/RADD_COMBINED_REICHE_DEFOR_YDAY_CORRETLY.tif")
radd_diff <- radd_sits - radd_reiche

terra::writeRaster(
    x = radd_diff,
    datatype = "INT1U",
    filename = "~/radd_report/RADD-DIFF.tif",
    gdal = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES",
             "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512")
)

values_sits <- unname(terra::values(radd_sits))
values_reiche <- unname(terra::values(radd_reiche))

digest::digest(values_sits, algo = "md5")
digest::digest(values_reiche, algo = "md5")

sum(terra::freq(radd_sits)["count"])
sum(terra::freq(radd_diff)["count"])

