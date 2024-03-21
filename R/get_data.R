# Setup -------------------------------------------------------------------
suppressWarnings({
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(rvest)
    library(httr)
    library(sf)
    library(ggplot2)
    library(terra)
    library(fs)
    library(purrr)
    library(rgbif)
    library(auk)
  })
})


# Get nest data -----------------------------------------------------------

# Scrape Osprey Nest Data
if (!file.exists("data/coords.rds")) {
  base.url <- "https://www.osprey-watch.org"
  pages <- 1:293
  urls <- paste0(base.url, "/nests?page=", pages)
  
  nests.main.content <- purrr::map(urls, function(url) {
    cat(format(Sys.time()), "Reading content from", url, "\n")
    page.content <- read_html(url)
    # Find all anchor tags within the "nests" div
    nest.ids <- html_nodes(page.content, ".nests a") %>% 
      html_attr("href") %>%
      unique() 
    nests <- paste0(base.url, nest.ids)
    Sys.sleep(1) # Give the server a rest
    out <- list(
      content=page.content,
      urls=nests
    )
    # Return content & urls
    return(out)
  })
  
  # Do some filtering, since some of the observations are definitely wrong
  correct.coords <- function(lat, lon) {
    if (lat > 90 | lat < -90 | lon > 180 | lon < -180) {
      return(c(NA, NA))
    } else if ((lat > 90 | lat < -90) & (lon >= -90 & lon <= 90)) {
      return(c(lon, lat))
    } else {
      return(c(lat, lon))
    }
  }
  
  coords <- purrr::map_df(nests.main.content, function(n) {
    content <- n$content
    text <- html_text(content, trim=T)
    # Regular expression to match the pattern "plotMyNest(id, latitude, longitude);"
    pattern <- "plotMyNest\\((\\d+),(\\-?\\d+\\.\\d+),(\\-?\\d+\\.\\d+)\\);"
    
    # Extract matches
    matches <- stringr::str_extract_all(text, pattern)[[1]]
    
    # Process each match to extract id, latitude, and longitude
    coord.data <- lapply(matches, function(match) {
      components <- stringr::str_match(match, pattern)
      # crds <- correct.coords(as.numeric(components[4]), as.numeric(components[3]))
      return(data.frame(
        id = components[2],
        lat = as.numeric(components[3]),
        long = as.numeric(components[4])
      ))
    })
    
    # Combine all data frames into a single data frame
    df <- do.call(rbind, coord.data)
  })
  
  urls <- purrr::map(nests.main.content, ~.x$urls) %>% unlist()
  coords$url <- urls
  coords <- coords %>% 
    filter(!(lat > 180 | lat < -180 | long > 180 | long < -180)) %>%
    setDT()
  
  saveRDS(coords, "data/coords.rds")
} else {
  coords <- readRDS("data/coords.rds")
}

if (!file.exists("data/nests_us.shp")) {
  # Convert to an sf object
  sf.df <- st_as_sf(coords, coords = c("long", "lat"), crs = 4326)
  
  # Load CONUS shapefile, get intersection with nests
  us.sf <- sf::read_sf("data/conus.shp") %>%
    dplyr::select(NAME, STATE_ABBR, geometry) %>%
    dplyr::rename(state=NAME, state_abbr=STATE_ABBR)
  nests.us <- st_intersection(sf.df, us.sf) %>%
    .[!st_is_empty(.), ]
  # Save
  sf::write_sf(nests.us, "data/nests_us.shp")
} else {
  us.sf <- sf::read_sf("data/conus.shp") %>%
    dplyr::select(NAME, STATE_ABBR, geometry) %>%
    dplyr::rename(state=NAME, state_abbr=STATE_ABBR)
  nests.us <- sf::read_sf("data/nests_us.shp")
}

# Filter to North Carolina
nc.sf <- sf::read_sf("data/North_Carolina_State_and_County_Boundary_Polygons.shp") %>%
  sf::st_transform(st_crs(us.sf)$wkt) %>%
  dplyr::select(County, geometry)
nests.nc <- nests.us[nests.us$state_abbr == "NC",]


# Raster Data -------------------------------------------------------------

r <- lapply(list.files("D:/data/final_us_rasters_2",  
                       full.names=T, pattern="\\.tif$"), rast) %>%
  purrr::reduce(c)
if (!file.exists("data/nc_rasters.tif")) {
  nc.r <- terra::crop(r, terra::vect(us.sf[us.sf$state_abbr=="NC",]), mask=T)
  nc.r$water_shore_dist <- min(nc.r$dist_to_water, nc.r$dist_to_shore)
  terra::writeRaster(nc.r, "data/nc_rasters.tif", overwrite=T)
} else {
  nc.r <- rast("data/nc_rasters.tif")
}

if (!file.exists("data/osprey_nests.csv")) {
  us.sf <- sf::read_sf("data/conus.shp")
  coords <- readRDS("data/coords.rds")
  .nests <- coords %>% 
    st_as_sf(sf_column_name="geometry", 
             coords=c("long", "lat"), crs=crs(r)) %>%
    sf::st_intersection(us.sf) %>%
    dplyr::select(url, NAME, geometry) %>%
    rename(state=NAME)
  
  raster.values <- terra::extract(r, .nests) %>% 
    setDT() %>%
    .[, ID := NULL]
  
  # raster.values[, na.vals := rowSums(is.na(raster.values))]
  # raster.values[na.vals > 0]
  # raster.values[, na.vals := NULL]
  dt <- cbind(.nests, raster.values, sf::st_coordinates(.nests)) %>%
    as.data.table() %>%
    .[, geometry := NULL] %>%
    setnames(old=c("X", "Y"), new=c("long", "lat"))
  dt.sf <- st_as_sf(dt, coords=c("long", "lat"), crs=crs(r))
  fwrite(dt, "data/osprey_nests.csv", row.names=F)
}

# Osprey Observation data (GBIF)
dwnld.tk <- function(key, extent) {
  tk <- tryCatch({
    occ_download(pred("taxonKey", key), 
                 pred_within(extent),
                 format="SIMPLE_CSV") %>%
      occ_download_wait()
  }, error=function(e) {
    browser()
    print(e)
  })
  d <- occ_download_get(tk$key, path="data", overwrite=T) %>%
    occ_download_import()
  
  utils::unzip(file.path("data", paste0(tk$key, ".zip")), 
               exdir=file.path("data", tk$key), overwrite=T) 
  file.remove(file.path("data", paste0(tk$key, ".zip")))
  
  d <- fread(file.path("data", tk$key, paste0(tk$key, ".csv")))
  unlink(file.path("data", tk$key), recursive=T)
  if (nrow(d) > 0) {
    fname <- gsub(" ", "_", d$species[[1]]) %>%
      gsub("[^[:alnum:]_]", "", .) %>%
      tolower() %>%
      paste0(., "_", key, ".csv.gz")
    f <- file.path("data", "species_obs", fname)
    fwrite(d, f)
  } 
  meta <- file.path("data", "species_obs", paste0(key, ".rds"))
  # Save metadata
  saveRDS(tk, meta)
  return(meta)
}

get.species.obs <- function(.data, extent) {
  if (!is.null(.data)) {
    tks <- .data %>%
      filter(publishingCountry == "US") %>%
      .$taxonKey %>%
      unique()
    if (length(tks) > 0) {
      out.files <- purrr::map2(tks, 1:length(tks), function(tk, i) {
        f <- file.path("data", "species_obs", paste0(tk, ".rds"))
        if (!file.exists(f)) {
          f <- dwnld.tk(tk, extent)
          cat(format(Sys.time()), "Saved metadata for", 
              paste0("[", i, "/", length(tks), "]"), "to", f, "\n")
        } else {
          cat(format(Sys.time()), paste0("[", i, "/", length(tks), "]"), 
              "already exists\n")
        }
        
      })
    } else {
      cat(format(Sys.time()), "No data for",  .data$species[[1]], "\n")
    }
  } else {
    cat(format(Sys.time()), "No data found\n")
  }
}

if (!file.exists("data/gbif_obs.csv.gz")) {
  
  .crs <- sf::st_as_text(st_crs(us.sf))
  .extent <- sf::st_bbox(us.sf)
  wkt.extent <- sprintf("POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s))",
                        .extent["xmin"], .extent["ymin"],
                        .extent["xmin"], .extent["ymax"],
                        .extent["xmax"], .extent["ymax"],
                        .extent["xmax"], .extent["ymin"],
                        .extent["xmin"], .extent["ymin"])
  
  species <- "Pandion haliaetus"
  if (!file.exists("data/osprey.occ.rds")) {
    osp.occ <- rgbif::occ_data(scientificName=species)
    saveRDS(osp.occ, "data/osprey.occ.rds")
  } else {
    osp.occ <- readRDS("data/osprey.occ.rds")
  }
  get.species.obs(osp.occ$data, wkt.extent)
  
  osp.dt <- purrr::map(
    list.files("data/species_obs", full.names=T, 
               pattern="\\.csv\\.gz"), 
    ~{
      fread(.x) %>%
        .[stateProvince == "North Carolina"] %>%
        .[!(basisOfRecord %in% c("PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN"))] %>%
        .[!is.na(decimalLatitude) & !is.na(decimalLongitude)] %>%
        .[, `:=` (eventDate=suppressWarnings(lubridate::as_date(eventDate)),
                  dateIdentified=NULL,
                  typeStatus=NULL)] %>%
        .[!is.na(eventDate) & eventDate >= as.Date("2010-01-01")] 
    }
  ) %>% rbindlist() 
  
  fwrite(osp.dt, "data/gbif_obs.csv.gz", row.names=F)
} else {
  osp.dt <- fread("data/gbif_obs.csv.gz")
}


# Osprey Observation Data (eBird)
# https://ebird.org/data/download
# Specify where your eBird datasets were downloaded to;
ebird.download.dirs <- list.dirs("D:/data/ebird/ebird_downloads")[-1]

# Specify where the outputs should be saved
ebird.output.dir <- "D:/data/ebird/ebird/ebird_outputs"


# Parse eBird downloads
for (dir in ebird.download.dirs) {
  # Set the path to EBD text files
  # auk_set_ebd_path(dir, overwrite = T)
  cat("Processing species observation data at", dir, "\n")
  # List .txt files that start with "ebd_"
  in.files <- list.files(path = dir, pattern = "^ebd_.*\\.txt$", full.names = T)
  
  .sampl <- in.files[grepl("*sampling\\.txt$", in.files)]
  in.file <- in.files[!grepl("*sampling\\.txt$", in.files)]
  
  # Use regular expression to extract state abbreviation
  state_abbreviation <- sub(".*_US-([A-Z]{2})_.*", "\\1", in.file)
  out.file <- file.path(ebird.output.dir, paste0(state_abbreviation, ".txt"))
  
  out.sampling.file <- file.path(ebird.output.dir, 
                                 paste0("sampling_", state_abbreviation, ".txt"))
  if (!file.exists(out.file) & !file.exists(out.sampling.file)) {
    # Read in the filtered data using the `auk` library, saving to `out.file`
    auk_ebd(in.file, .sampl) %>%
      auk_species(species = species) %>% 
      auk_complete() %>% # Add this to keep only complete checklists
      auk_filter(file = out.file, file_sampling = out.sampling.file, 
                 overwrite=T, execute=T)
    
    # Remove "Problem" records
    df <- readr::read_delim(out.file, delim = "\t", 
                            show_col_types = F) %>% 
      suppressWarnings()
    df.samp <- readr::read_delim(out.sampling.file, delim = "\t", 
                                 show_col_types = F) %>%
      suppressWarnings()
    p <- problems(df)
    p.s <- problems(df.samp)
    if (nrow(p) > 0) {
      df <- df %>% slice(-p$row)
      readr::write_delim(df, out.file, delim="\t")
    }
    if (nrow(p.s) > 0) {
      df.samp <- df.samp %>% slice(-p$row)
      readr::write_delim(df.samp, out.sampling.file, delim="\t")
    }
  }
}

# Some basic pre-processing
preprocess.obs <- function(data.path) {
  data <- readr::read_delim(data.path, delim = "\t", show_col_types = F) %>%
    suppressMessages()
  names(data) <- gsub(" ", "\\.", tolower(names(data)))
  data <- data %>%
    filter(observation.date >= as.Date("2016-01-01") & 
             observation.date < as.Date("2024-01-01") &
             approved == 1 & observation.count != "X") %>%
    mutate(observation.date = as.Date(observation.Date)) %>%
    dplyr::select(common.name, observation.count, observation.date, 
                  latitude, longitude) %>% 
    group_by(common.name, observation.date, latitude, longitude) %>% 
    summarize(observation.count = sum(as.numeric(observation.count), na.rm=T),
              .groups="keep") %>%
    ungroup() %>%
    as.data.table()
  return(data)
}

if (!file.exists("data/osprey_observations.csv.gz")) {
  osp.dt <- preprocess.obs(file.path(ebird.output.dir, "NC.txt"))
  fwrite(osp.dt, "data/osprey_observations.csv.gz", row.names=F)
} else {
  osp.dt <- fread("data/osprey_observations.csv.gz")
}

# Bounding Box & Padded Bounding Box

expand.bbox <- function(bbox, x=2, y=2) {
  bbox <- c(xmin=bbox$xmin - x/2, ymin=bbox$ymin - y/2, 
            xmax=bbox$xmax + x/2, ymax=bbox$ymax + y/2)
  class(bbox) <- "bbox"
  return(bbox)
}

bb <- sf::st_bbox(nc.sf) 

bb.sf <- bb %>%
  expand.bbox(1, 1) %>% 
  st_as_sfc() %>%
  st_set_crs(st_crs(nc.sf)) %>%
  sf::st_make_valid()# %>%
# st_intersection(us.sf, .) %>%
# .[!st_is_empty(.), ]
