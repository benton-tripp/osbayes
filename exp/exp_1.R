# Setup -------------------------------------------------------------------
suppressWarnings({
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(sf)
    library(ggplot2)
    library(terra)
    library(purrr)
    library(spatstat)
  })
})

# Load data
source("R/get_data.R")
source("R/ppp_preprocessing.R")


# Initial Exploratory Plots -----------------------------------------

# Plot US Observations
ggplot(us.sf) +
  geom_sf(fill = "white", color = "gray") +
  geom_sf_text(aes(label = state), size = 3, check_overlap = T) +
  geom_sf(data=nests.us, color="darkred", size=1.25) + 
  theme_void()

# Plot NC Observations
ggplot(nc.sf) +
  geom_sf(fill = "white", color = "gray") +
  geom_sf(data=us.sf[us.sf$state_abbr=="NC",], fill="darkblue", alpha=0.1) +
  geom_sf_text(aes(label = County), size = 3, check_overlap = T) +
  geom_sf(data=nests.nc, color="darkred", size=1.25) + 
  theme_void()

# For README
if (!file.exists("figures/nc_osprey_nest_sites.png")) {
  nc.nest.plt <- ggplot() + 
    geom_sf(data=us.sf, fill="#c7cd9d", color="#5a5a48") +
    geom_sf(data=us.sf[us.sf$state_abbr=="NC",], fill="#8f955e",
            linewidth=0.5, color="#7e844d") +
    geom_sf(data=nests.nc, color="#142919", size=1.5) + 
    coord_sf(xlim=c(bb$xmin-.5, bb$xmax+.5), 
             ylim=c(bb$ymin-.5, bb$ymax+.5), expand=F) + 
    xlab("Longitude") + ylab("Latitude") + 
    ggtitle("NC Osprey Nesting Sites") + 
    theme(panel.grid.major=element_line(color="#bcbcb8", 
                                        linetype="dashed", 
                                        linewidth=0.5), 
          panel.background=element_rect(fill="#748c82"),
          panel.border=element_rect(colour="#142919", linewidth=0.5,
                                    fill=NA), 
          plot.background=element_rect(color="#142919", linewidth=1, 
                                       fill="#fbfbf1"))
  ggsave("figures/nc_osprey_nest_sites.png", plot=nc.nest.plt, width=6.12, height=4.6)
}


# Poisson Point Process Modeling ------------------------------------------

# Create Poisson Point Process Object
region.list <- convert.to.list.format(nc.sf)

nest.coords <- sf::st_coordinates(nests.nc) %>%
  as.data.table() %>% 
  unique()

# Convert the data to a ppp objects
locations <- spatstat.geom::ppp(nest.coords$X, nest.coords$Y, 
                                poly=region.list) 

# Fit a basic Poisson Point Process model using just location
baseline.fit <- ppm(locations ~ x + y, rbord=.05, method="mpl", emend=T)


region.dt <- nc.r %>%
  terra::as.data.frame(xy=T) %>%
  as.data.table() 

region.imgs <- list(
  x=as.im(region.dt[, .(x, y, z=x)]),
  y=as.im(region.dt[, .(x, y, z=y)])
)

covariates <- c("x", "y", names(nc.r)) %>%
  purrr::set_names() %>%
  purrr::map(function(rname) {
    if (!(rname) %in% c("x", "y")) {
      # Capture the current raster in the iteration
      rc <- nc.r[[rname]]
      
      # Return a closure (function) that captures 'rc'
      fun <- function(x, y) {
        p <- data.table(x = x, y = y)
        extracted.values <- terra::extract(rc, p)
        return(extracted.values[, 2])
      }
    } else {
      fun <- function(x, y) {
        if (rname == "x") return(x) else return(y)
      }
    }
    return(fun)
  })

selected.covariates <- c("dem", "water_shore_dist", "urban_imperviousness", "canopy") %>%
  covariates[.]
.f <- as.formula("~ dem + water_shore_dist + urban_imperviousness + canopy")
fit <- ppm(locations, trend=.f, covariates=selected.covariates,
           rbord=.05, method="mpl", emend=T, use.gam=T)

pred <- spatstat.model::predict.ppm(
  fit, 
  locations=region.dt,
  type="trend")
region.dt[, intensity := pred]

pred.r <- rast(region.dt)
crs(pred.r) <- crs(nc.r)
pred.r <- terra::crop(pred.r,
                      terra::vect(us.sf[us.sf$state_abbr=="NC",]), 
                      mask=T)

ggplot() + 
  geom_raster(data=region.dt, aes(x=x, y=y, fill=intensity)) + 
  geom_sf(data=us.sf[us.sf$state_abbr != "NC",], 
          fill="#c7cd9d", color="#5a5a48") +
  coord_sf(xlim=c(bb$xmin-.5, bb$xmax+.5), 
           ylim=c(bb$ymin-.5, bb$ymax+.5), expand=F) + 
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("Estimated NC Osprey Nesting Site Intesity") + 
  theme(panel.grid.major=element_line(color="#bcbcb8", 
                                      linetype="dashed", 
                                      linewidth=0.5), 
        panel.background=element_rect(fill="#748c82"),
        panel.border=element_rect(colour="#142919", linewidth=1,
                                  fill=NA), 
        plot.background=element_rect(color="#142919", linewidth=1, 
                                     fill="#fbfbf1"))  +
  viridis::scale_fill_viridis(option="A", direction = -1)


# Osprey Observations
osp.dt <- fread("data/osprey_observations.csv.gz")






