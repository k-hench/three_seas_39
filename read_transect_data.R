# loading all needed libraries
# -----------------------------------
library(tidyverse)
library(sf)
library(rnaturalearth)
library(scatterpie)
library(paletteer)
library(ggforce)
library(ggsn)
library(hypoimg)
library(lubridate)

# setting up custom functions
# -----------------------------------

# function to crop shapefiles to plot extent 
crp <- function (poly, xlim = xlim_boc, ylim = ylim_boc) {
  st_intersection(poly, 
                  st_set_crs(st_as_sf(as(raster::extent(xlim[1],
                                                        xlim[2],
                                                        ylim[1],
                                                        ylim[2]), 
                                         "SpatialPolygons")),
                             st_crs(poly)))
}

# function to turn a data.frame (with the columns
# "Longitude" and "Latitude") into a spatial object
tibble_to_sf <- function(tib, crs = 4326){
  tib %>%
    st_as_sf(., coords = c("Longitude","Latitude")) %>%
    st_set_crs(., crs)
}

# funtion to fill empty cells (NA) with zeros 
replace_all_na <- function(tib, replace = 0, exclude = c("Date", "location", "group",
                                                         "habitat_type", "depth")){
  nm_replace <- names(tib)[!(names(tib) %in% exclude)]
  replace_list <- rep(replace,length(nm_replace)) %>% 
    set_names(nm = nm_replace) %>%
    as.list()
  
  tib %>% replace_na(replace = replace_list)
}

# -----------------------------------
# the actual script
# -----------------------------------

# setting the extent of the map
xlim_boc <- c(-82.5, -82)
ylim_boc <- c(9.1, 9.47)


# reading the panama shape files and
# cropping it to the plot extent
# (downloaded from https://www.gadm.org/download_country_v3.html)
bocas <- read_sf('data/PAN_adm0.shp') %>%
  crp()

# reading in the dive sites for the gps locations
sites <- read_tsv('data/dive_spots - Sheet1.tsv') %>%
  filter(!duplicated(Site)) %>%
  select(Site, Latitude, Longitude)

# reading in the transect data, replacing NAs with zeros
# and merge with gps positions
transects <- read_tsv('data/Fish_surveys - Sheet1.tsv',
                      col_types = str_c(c('ccdc',rep('d', 78)),collapse = '')) %>%
  mutate(Site = location %>% str_to_lower() %>% str_remove('_mangrove'),
         Date = Date %>%
           str_replace(pattern = "([0-9]*)/([0-9]*)/([0-9]*)",
                       replacement = "\\3-\\1-\\2") %>%
           as_date(Date)) %>%
  replace_all_na() %>%
  left_join(sites)

# create a spatial object from the transects
transects_sf <- transects %>%
  tibble_to_sf()

# reading in the fish list for the fish families and create an ID column
# for merging with the transect data
fish_list <- read_tsv('data/Fish_list_2020 - Sheet1.tsv') %>%
  select(Common_name:`Family (latin name)`) %>%
  # format the fish id column which is going to be
  # used for the merging with the transect data
  mutate(fish_id = Common_name )# %>%
           # str_to_lower() %>% 
           # str_replace_all(pattern = ' ', 
           #                 replacement = '_'))

# pie shift is set to define the distance between
# the "mangrove" and "reef" pie for each location 
pie_shift <- .03

# merge transect data and fish list
summary_by_family <- transects %>%
  # we need unique identifiers for each transect,
  # so we are going to assign a transet_nr_within_group 
  group_by(location, group) %>%
  mutate(transet_nr_within_group  = row_number(),
         transect_id = str_c(location, group,
                             transet_nr_within_group,
                             sep = '_')) %>%
  ungroup() %>%
  # then we transform the data into "long format" to be able
  # to merge the family names based on the common name
  pivot_longer(cols = bridled_goby:porkfish,
               names_to = 'common_name') %>% 
  # from the common name we create a fish id that is in the EXACT same
  # format as the fish id column of the fish list data frame
  # (this requires some reformating, eg. dropping '_juvenile')
  mutate(fish_id = common_name %>% 
           str_remove("_juvenile|juvenile_|initial_|intermediate_|juv_"))  %>%# View()
  #          str_replace(pattern = 'brown_chrosmis_damselfish',
  #                      replacement = 'brown_chromis') %>%
  # str_replace(pattern = 'graysby_grouper',
  #             replacement = 'graysby') %>%
  # str_replace(pattern = 'doctor_surgeonfish',
  #             replacement = 'doctorfish')) %>%
  # merge with the fish list
  left_join(fish_list) %>% 
     # filter(!duplicated(fish_id)) %>%
     # View()
  # drop all row that do not have an entry in the
  # fish list (typos & missing)
  filter(!is.na(`Latin name`)) %>%
  # sum the count for each family at each habitat type of each location
  group_by(`Family (latin name)`, Site, habitat_type) %>%
  summarise(n = sum(value),Latitude = Latitude[[1]], Longitude = Longitude[[1]]) %>%
  # transform back into "wide format"
  pivot_wider(names_from = `Family (latin name)`, values_from = n) %>%
  # adjust positions of pies
  mutate(grp = str_c(Site, habitat_type, sep = '_'),
         Longitude = ifelse(habitat_type == 'Reef', Longitude - .5 * pie_shift, Longitude + .5 * pie_shift),
         Latitude = ifelse(habitat_type == 'Reef', Latitude - .5 * pie_shift, Latitude +  .5 *pie_shift))

# select the names of all present fish families
fish_columns <- names(summary_by_family)[!(names(summary_by_family) %in%
                                             c("Date", "Site", "habitat_type", "Latitude",
                                               "Longitude", "grp"))]

# set color coding for habitat type
clr_habitat <- c(rgb(0,0,0), rgb(1,1,1)) %>% 
  set_names(nm = c('Reef', 'Mangrove'))

# import the compass image
compass <- hypo_read_svg('north.svg')

# -----------------------------------
# plotting the map
# -----------------------------------

# initalize the plot
ggplot()+
  # add the panama coastline
  geom_sf(data = bocas)+
  # add the sampling sites (as dots)
  geom_sf(data = transects_sf)+
  # add the pies
  geom_scatterpie(aes(x = Longitude, y = Latitude ,
                      r = .6 * pie_shift , group = grp),
                  data = summary_by_family, color = rgb(1,1,1,0),
                  cols = fish_columns) +
  # add a color coded ring around the pie to indicate habitat type
  geom_circle(data = summary_by_family,
              aes(x0 = Longitude, y0 = Latitude ,
                  r = .6 * pie_shift, color = habitat_type )) +
  # add the scalebar
  scalebar(x.min = -82.5, x.max = -82.3,
           y.min = 9.13, y.max = 9.17,transform = TRUE,
           dist = 10, dist_unit = "km",
           st.dist = .35, st.size = 3.5,
           border.size = .1, height = .15)+
  # add the compass
  annotation_custom(grob = compass,
                    xmin = -82.49, xmax = -82.44,
                    ymin = 9.4, ymax = 9.45)+
  # set the color palette for the pies
  scale_fill_manual(values = paletteer_c(n = length(fish_columns),
                                         palette = "ggthemes::Red-Green-Gold Diverging"))+
  # set the color palette for the habitat types
  scale_color_manual(values = clr_habitat, guide = FALSE
                     )+
  # remove any "extra space" on the x and y axis
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # foramt the legend
  guides(fill = guide_legend(title = 'Family',
                             title.position = 'top')#,
         # color = guide_legend(title = 'Habitat type',
         #                     title.position = 'top',
         #                     ncol = 1,override.aes = list(shape = 2))
         )+
  # twaek the plot appearance
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.key.size = unit(14,'pt'),
        legend.key.height = unit(14,'pt'),
        legend.background = element_rect(fill = 'lightgray', 
                                         color = rgb(1,1,1,0)),
        legend.position = 'bottom')

library(ggstance)
summary_by_family %>%
  pivot_longer(cols = Acanthuridae:Urotrygonidae,
               names_to = 'Family', values_to = 'n') %>%
  # mutate(Family = fct_reorder(Family, n)) %>%
  ggplot(aes(y = fct_reorder(Family, n), x = n, fill = Family))+
  geom_barh(stat = 'identity')+
  facet_grid(habitat_type ~ Site)+
  scale_fill_manual(values = paletteer_c(n = length(fish_columns),
                                         palette = "ggthemes::Red-Green-Gold Diverging"))+
  guides(fill = guide_legend(title = 'Family',
                             title.position = 'top',nrow = 3))+
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.key.size = unit(14,'pt'),
        legend.key.height = unit(14,'pt'),
        legend.background = element_rect(fill = 'lightgray', 
                                         color = rgb(1,1,1,0)),
        legend.position = 'bottom')
