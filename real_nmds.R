# loading all needed libraries
# -----------------------------------
library(tidyverse)
library(ggalt)
library(vegan)
library(lubridate)

# setting up custom functions
# -----------------------------------
# funtion to fill empty cells (NA) with zeros 
replace_all_na <- function(tib, replace = 0,
                           exclude = c("Date", "location", "group",
                                       "habitat_type", "depth")){
  nm_replace <- names(tib)[!(names(tib) %in% exclude)]
  replace_list <- rep(replace,length(nm_replace)) %>% 
    set_names(nm = nm_replace) %>%
    as.list()
  
  tib %>% replace_na(replace = replace_list)
}

# funtion to export nmds results for plotting with ggplot
export_nmds <- function(tib, nmds, score_name = "Spot", species_name = "Species"){
  # doublecheck that rownames are uniqe,
  # otherwise the merge at the end duplicates entries 
  if(any(duplicated(rownames(scores(nmds))))){
    stop(paste("Sorry, the rownames of the distance matrix need to be unique.
               Currently that is not the case:\n\n ", 
               paste(c('!',' ')[2 - (duplicated(rownames(scores(nmds))) | duplicated(rownames(scores(nmds)),
                                                                                     fromLast = TRUE))],
                     rownames(scores(nmds)),
                     collapse = '\n  ')))
  }
  
  # export obesravtions in nMDS space
  data_scores <- as.data.frame(scores(nmds)) %>% 
    as_tibble(rownames = score_name)
  # export species in nMDS space
  data_species <- as.data.frame(scores(nmds, "species"))%>% 
    as_tibble(rownames = species_name)
  
  list(spots = tib %>% left_join(data_scores),
       species = data_species)
}
# -----------------------------------
# the actual script
# -----------------------------------
  
# reading in the transect data and replacing NAs with zeros
transects <- read_tsv('data/Fish_surveys - Sheet1.tsv',
                      col_types = str_c(c('ccdc',rep('d', 78)),collapse = '')) %>%
  # we need unique identifiers for each transect,
  # so we are going to assign a transet_nr_within_group 
  group_by(location, group) %>%
  mutate(Site = location %>% str_to_lower() %>% str_remove('_mangrove'),
         transet_nr_within_group  = row_number(),
         Date = Date %>%
           str_replace(pattern = "([0-9]*)/([0-9]*)/([0-9]*)",
                       replacement = "\\3-\\1-\\2") %>%
           as_date(Date),
         Transect_id = str_c(Site, habitat_type, group,
                             transet_nr_within_group ,sep = '_')) %>%
  ungroup() %>%
  replace_all_na()

# alternative 1:
# removing empty transects
transects <- transects %>%
  mutate(total_count = transects %>% 
           rowwise() %>% 
           select(-(Date:depth),-(Site:Transect_id)) %>% rowSums())

transects <- transects %>%
  filter(total_count != 0)

data_matrix <- transects %>%
  select(bridled_goby:scrawled_cowfish) %>%
  as.matrix()

# alternative 2:
# merge by group and fish id
names(transects)
transects <- transects %>%
  pivot_longer(cols = bridled_goby:porkfish,
               names_to = 'common_name',
               values_to = 'n') %>% 
  mutate(fish_id = common_name %>% 
           str_remove("_juvenile|juvenile_|initial_")) %>% 
  group_by(Site, habitat_type, group, fish_id) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = fish_id,
              values_from = n) %>%
  mutate(Transect_id = str_c(Site, habitat_type, group, sep = '_'))


transects <- transects %>%
  mutate(total_count = transects %>% 
           rowwise() %>% 
           select(-(Site:group),-Transect_id) %>%
           rowSums())

transects <- transects %>%
  filter(total_count != 0)

# names(transects)

data_matrix <- transects %>%
  select(atlantic_spadefish:yellowtail_snapper) %>%
  as.matrix()

rownames(data_matrix) <- transects$Transect_id

data_nmds <- metaMDS(data_matrix, k = 2, trymax = 999, distance='bray')

plot(data_nmds)


data_nmds_results <- export_nmds(transects, data_nmds,score_name = 'Transect_id')

p1 <- ggplot(data_nmds_results$spots,
       aes(x = NMDS1, y = NMDS2,
           group = habitat_type)) + 
  coord_equal() +
  #  geom_text(aes(label = Spot)) +
  geom_hline(yintercept = 0, color = 'lightgray', alpha = .4) +
  geom_vline(xintercept = 0, color = 'lightgray', alpha = .4) +
  geom_text(inherit.aes = FALSE,
            data = data_nmds_results$species,
            aes(x = NMDS1, y = NMDS2, label = Species)) +
  geom_point(aes(colour = habitat_type), size = 2) +
  geom_encircle(aes(colour = habitat_type, fill = habitat_type), 
                s_shape = 1, alpha = 0.2, size = 1, expand = 0)+
  scale_x_continuous(expand = c(.01,.05))+
  theme(legend.position = 'bottom',
        panel.grid = element_blank())



p2 <- ggplot(data_nmds_results$spots,
       aes(x = NMDS1, y = NMDS2,
           group = Site)) + 
  coord_equal() +
  #  geom_text(aes(label = Spot)) +
  geom_hline(yintercept = 0, color = 'lightgray', alpha = .4) +
  geom_vline(xintercept = 0, color = 'lightgray', alpha = .4) +
  geom_text(inherit.aes = FALSE,
            data = data_nmds_results$species,
            aes(x = NMDS1, y = NMDS2, label = Species)) +
  geom_point(aes(colour = Site), size = 2) +
  geom_encircle(aes(colour = Site, fill = Site), 
                s_shape = 1, alpha = 0.2, size = 1, expand = 0)+
  scale_x_continuous(expand = c(.01,.05))+
  theme(legend.position = 'bottom',
        panel.grid = element_blank())


library(patchwork)
p1 + p2
