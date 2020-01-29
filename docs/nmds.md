---
output: html_document
editor_options:
  chunk_output_type: console
---
# nMDS and Permanova



## External pacakges

The script starts by loading all needed libraries.
(Technically, the package lubridate is not important, I only used it to propperly format the date column, but the information is never used... )


```r
# loading all needed libraries
# -----------------------------------
library(tidyverse)
library(ggalt)
library(vegan)
library(lubridate)
```

## Custom functions

To make the code mode easily understandable (hopefully), we define some helper functions upfront (as opposed to "in the middle of the process"):

At some point we will need to fill the empty celly of our data sheet with zeros - by default R will fill empty cells with `NA`.
For this we create a function, that will replace `NA` in all columns of a data frame except for a specific subset of columns (the "non-fish" columns).


```r
# setting up custom functions
# -----------------------------------
# funtion to fill empty cells (NA) with zeros 
replace_all_na <- function(tib, replace = 0,
                           exclude = c("Date", "location", "group",
                                       "habitat_type", "depth")){
  
  # get all fish columns ("non-non-fish-columns") 
  nm_replace <- names(tib)[!(names(tib) %in% exclude)]
  
  # create a replacement list in the form of: list(fish_1 = 0, fish_2 = 0, ...)
  replace_list <- rep(replace,length(nm_replace)) %>% 
    set_names(nm = nm_replace) %>%
    as.list()
  
  # apply replacement list
  tib %>% replace_na(replace = replace_list)
}
```

To make plotting the nMDS results with ggplot a little easier, we create a export function that will take the nMDS results as input and return a list with two data sets (`output_name$spots` - the transect scores and `output_name$species` - the species "influence").


```r
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
```

## Loading data

Now, we are all set up and can start with the actual work.

First, we need to import the survey data into R.

Directly while importing we deal with some issues:

- we extract the *pure* site name from the `location` column (dropping the `_mangrove` suffix)
- we proppely format the date (actually transforming the column from `character` into `date`-type)
- since we have several transects per group, we number the transects within each group
- from the columns `Site`, `habitat_type`, `group` and `transet_nr_within_group` we create a unique identifier for each transect
- finally, we use the previously prepared function to replace all empty fish-cells with zeros

<p style='color:#f0a830'>Beware of the `col_types = str_c(c('ccdc',rep('d', 78)),collapse = '')` part: here we define the column types (c = character, d = "double"/number). So I the original google sheet changes (eg. by a merge og "blue hamlet" and "black hamlet") the `78` needs to be updated to *the total number of columns* - 4!</p>



```r
# ----------------------------------
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
```

```
## Warning: Missing column names filled in: 'X74' [74]
```


After a first check, I realised that we will need to merge the observations within each group (eg. line 2 & 3 or line 4 & 5 of the original google sheet) - here we merge several rows of the data set.

We also need to merge the observations over different stages (eg. juvenile_bluehead_wrasse & bluehead_wrasse) - here we merge several columns.

To do this merging we are going to apply some tidyverse-magic (using functions from the packages [tidyr](https://tidyr.tidyverse.org/) and [dplyr](https://dplyr.tidyverse.org/)).

If you have difficulties following these rather complex steps, it might be helpful to highlight and execute just parts of the code (to see what is going on you will allways need to start with `transects %>%` and NOT highlight the `transects <- ` bit):

![](img/highlight.png)



```r
# merge by group and fish id
transects <- transects %>%
  # transfrom the data from "wide" into "long" format
  # note that we select the fish columns here, so the bridled_goby:porkfish
  # part might need to be updated if the google cheet changes to the
  # form of fist_fish:last_fish
  pivot_longer(cols = bridled_goby:porkfish,
               names_to = 'common_name',
               values_to = 'n') %>%
  # create a new column with the fish species (irrespective of stage)
  mutate(fish_id = common_name %>% 
           str_remove("_juvenile|juvenile_|initial_")) %>% 
  # set the logic to summarise over (we want one entry per species
  # for every site, habitat_type and group)
  group_by(Site, habitat_type, group, fish_id) %>%
  # summ the counts according to our grouping
  summarise(n = sum(n)) %>%
  # remove the grouping (the effect of this is not visible,
  # but if we don't do this the data table might behave strange later on..)
  ungroup() %>%
  # reformat back to wide format (note that now the species are sorted alphabetically)
  pivot_wider(names_from = fish_id,
              values_from = n) %>%
  # Since we now have only a single entry per group,
  # we update our transect identifiers
  mutate(Transect_id = str_c(Site, habitat_type, group, sep = '_'))
```

For some reason, some transects are completely empty (0 for all species).
These are not informative in terms of community composition and need to be dropped from the data table (otherwise the nMDS will fail because there is no way to compute ecological distances with those transects).
To check for this, we create a new column with the sum of all fish for each transect.


```r
transects <- transects %>%
  mutate(total_count = transects %>% 
           rowwise() %>% 
           select(-(Site:group),-Transect_id) %>%
           rowSums())
```

Now we filter all transects that have a total count of zero.


```r
transects <- transects %>%
  filter(total_count != 0)
```

## Run nMDS

So far, we have worked with our data in the format of a `data.frame`, yet the nMDS expects a `matrix` (this just some internal stuff in R that has to do with the data stucture and you don't rally have to worry about...).

To prepare the nMDS input we select only the fish columns (<span style='color:#f0a830'>note that the column order has changed and the fish are now sorted alphabetically</span>) and convert these columns into a matrix.

We store the transect ids in the `rownames` of the matrix (this is no data *within* the matrix and hence not considdered for the analysis).


```r
data_matrix <- transects %>%
  select(atlantic_spadefish:yellowtail_snapper) %>%
  as.matrix()

rownames(data_matrix) <- transects$Transect_id
```

Now we can run the nMDS on our data and export the results for ggplot.


```r
data_nmds <- metaMDS(data_matrix, k = 2, trymax = 999, distance='bray')
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.2054951 
## Run 1 stress 0.2073328 
## Run 2 stress 0.2110472 
## Run 3 stress 0.2089344 
## Run 4 stress 0.2112255 
## Run 5 stress 0.2072524 
## Run 6 stress 0.2076932 
## Run 7 stress 0.2061955 
## Run 8 stress 0.2277037 
## Run 9 stress 0.2251799 
## Run 10 stress 0.2120436 
## Run 11 stress 0.2061956 
## Run 12 stress 0.2158478 
## Run 13 stress 0.2061961 
## Run 14 stress 0.2232418 
## Run 15 stress 0.2089348 
## Run 16 stress 0.212634 
## Run 17 stress 0.2100794 
## Run 18 stress 0.229484 
## Run 19 stress 0.2061973 
## Run 20 stress 0.2226701 
## Run 21 stress 0.2269043 
## Run 22 stress 0.208934 
## Run 23 stress 0.2112546 
## Run 24 stress 0.2089026 
## Run 25 stress 0.2084581 
## Run 26 stress 0.2229818 
## Run 27 stress 0.2101229 
## Run 28 stress 0.225972 
## Run 29 stress 0.2114888 
## Run 30 stress 0.2158514 
## Run 31 stress 0.2115015 
## Run 32 stress 0.2102373 
## Run 33 stress 0.2061962 
## Run 34 stress 0.2084403 
## Run 35 stress 0.2084372 
## Run 36 stress 0.210387 
## Run 37 stress 0.2089371 
## Run 38 stress 0.2084382 
## Run 39 stress 0.2056559 
## ... Procrustes: rmse 0.01064424  max resid 0.05805437 
## Run 40 stress 0.2072531 
## Run 41 stress 0.2230131 
## Run 42 stress 0.2089486 
## Run 43 stress 0.2226318 
## Run 44 stress 0.2124895 
## Run 45 stress 0.2084525 
## Run 46 stress 0.2185491 
## Run 47 stress 0.212617 
## Run 48 stress 0.2116721 
## Run 49 stress 0.2099511 
## Run 50 stress 0.2310483 
## Run 51 stress 0.2125223 
## Run 52 stress 0.2112245 
## Run 53 stress 0.2084565 
## Run 54 stress 0.2061955 
## Run 55 stress 0.2112547 
## Run 56 stress 0.2109999 
## Run 57 stress 0.2126238 
## Run 58 stress 0.2069518 
## Run 59 stress 0.2124165 
## Run 60 stress 0.2101305 
## Run 61 stress 0.2134171 
## Run 62 stress 0.2056535 
## ... Procrustes: rmse 0.01052152  max resid 0.0585433 
## Run 63 stress 0.2073191 
## Run 64 stress 0.2076925 
## Run 65 stress 0.2182903 
## Run 66 stress 0.2159858 
## Run 67 stress 0.2061957 
## Run 68 stress 0.2104501 
## Run 69 stress 0.2123883 
## Run 70 stress 0.223215 
## Run 71 stress 0.2089352 
## Run 72 stress 0.2137481 
## Run 73 stress 0.2061977 
## Run 74 stress 0.2056535 
## ... Procrustes: rmse 0.01055998  max resid 0.05926682 
## Run 75 stress 0.2084524 
## Run 76 stress 0.2114861 
## Run 77 stress 0.2061957 
## Run 78 stress 0.2089639 
## Run 79 stress 0.2069527 
## Run 80 stress 0.206213 
## Run 81 stress 0.2110489 
## Run 82 stress 0.2089347 
## Run 83 stress 0.2123777 
## Run 84 stress 0.2103604 
## Run 85 stress 0.2126263 
## Run 86 stress 0.2123735 
## Run 87 stress 0.2123798 
## Run 88 stress 0.2056549 
## ... Procrustes: rmse 0.01064887  max resid 0.06014341 
## Run 89 stress 0.2056536 
## ... Procrustes: rmse 0.01056363  max resid 0.0588167 
## Run 90 stress 0.2084382 
## Run 91 stress 0.2089348 
## Run 92 stress 0.2104627 
## Run 93 stress 0.2159874 
## Run 94 stress 0.2061963 
## Run 95 stress 0.2110054 
## Run 96 stress 0.2105451 
## Run 97 stress 0.2343238 
## Run 98 stress 0.215314 
## Run 99 stress 0.2250141 
## Run 100 stress 0.2112248 
## Run 101 stress 0.2250305 
## Run 102 stress 0.2076914 
## Run 103 stress 0.2089026 
## Run 104 stress 0.206952 
## Run 105 stress 0.2160433 
## Run 106 stress 0.2110327 
## Run 107 stress 0.2089351 
## Run 108 stress 0.2061959 
## Run 109 stress 0.2310673 
## Run 110 stress 0.2100788 
## Run 111 stress 0.2061956 
## Run 112 stress 0.2167566 
## Run 113 stress 0.206197 
## Run 114 stress 0.2152961 
## Run 115 stress 0.2171266 
## Run 116 stress 0.206196 
## Run 117 stress 0.2143276 
## Run 118 stress 0.2159875 
## Run 119 stress 0.2153079 
## Run 120 stress 0.211664 
## Run 121 stress 0.2061955 
## Run 122 stress 0.2111632 
## Run 123 stress 0.2228563 
## Run 124 stress 0.2061956 
## Run 125 stress 0.2234943 
## Run 126 stress 0.2076915 
## Run 127 stress 0.2079455 
## Run 128 stress 0.2260427 
## Run 129 stress 0.2108587 
## Run 130 stress 0.2247557 
## Run 131 stress 0.2260429 
## Run 132 stress 0.2224012 
## Run 133 stress 0.2125277 
## Run 134 stress 0.2123009 
## Run 135 stress 0.2072541 
## Run 136 stress 0.2107816 
## Run 137 stress 0.2269081 
## Run 138 stress 0.2267632 
## Run 139 stress 0.2238753 
## Run 140 stress 0.2089608 
## Run 141 stress 0.2101249 
## Run 142 stress 0.2110071 
## Run 143 stress 0.2062129 
## Run 144 stress 0.2271357 
## Run 145 stress 0.2089359 
## Run 146 stress 0.2061963 
## Run 147 stress 0.2072647 
## Run 148 stress 0.2167583 
## Run 149 stress 0.208935 
## Run 150 stress 0.2161062 
## Run 151 stress 0.20566 
## ... Procrustes: rmse 0.005987957  max resid 0.02943882 
## Run 152 stress 0.2089053 
## Run 153 stress 0.2100831 
## Run 154 stress 0.2163064 
## Run 155 stress 0.2076934 
## Run 156 stress 0.2076075 
## Run 157 stress 0.2111609 
## Run 158 stress 0.2110051 
## Run 159 stress 0.206196 
## Run 160 stress 0.2115747 
## Run 161 stress 0.2114879 
## Run 162 stress 0.2108576 
## Run 163 stress 0.2272283 
## Run 164 stress 0.2061958 
## Run 165 stress 0.2281372 
## Run 166 stress 0.2112237 
## Run 167 stress 0.2124408 
## Run 168 stress 0.2084526 
## Run 169 stress 0.2137899 
## Run 170 stress 0.2061963 
## Run 171 stress 0.2056566 
## ... Procrustes: rmse 0.0105751  max resid 0.06038153 
## Run 172 stress 0.2056534 
## ... Procrustes: rmse 0.01054811  max resid 0.05924621 
## Run 173 stress 0.2105566 
## Run 174 stress 0.2101229 
## Run 175 stress 0.2084381 
## Run 176 stress 0.2128589 
## Run 177 stress 0.2076207 
## Run 178 stress 0.2121746 
## Run 179 stress 0.2061955 
## Run 180 stress 0.2160432 
## Run 181 stress 0.223066 
## Run 182 stress 0.2056537 
## ... Procrustes: rmse 0.01058809  max resid 0.05958554 
## Run 183 stress 0.222882 
## Run 184 stress 0.2084391 
## Run 185 stress 0.2265079 
## Run 186 stress 0.2132956 
## Run 187 stress 0.2185497 
## Run 188 stress 0.2105417 
## Run 189 stress 0.2133624 
## Run 190 stress 0.2124931 
## Run 191 stress 0.2123818 
## Run 192 stress 0.2226805 
## Run 193 stress 0.2061958 
## Run 194 stress 0.2115045 
## Run 195 stress 0.206196 
## Run 196 stress 0.2062142 
## Run 197 stress 0.2089391 
## Run 198 stress 0.2281036 
## Run 199 stress 0.2162323 
## Run 200 stress 0.2123821 
## Run 201 stress 0.2155047 
## Run 202 stress 0.2089485 
## Run 203 stress 0.2076913 
## Run 204 stress 0.2089355 
## Run 205 stress 0.2109996 
## Run 206 stress 0.2084579 
## Run 207 stress 0.2111613 
## Run 208 stress 0.2112233 
## Run 209 stress 0.2089345 
## Run 210 stress 0.212189 
## Run 211 stress 0.2177879 
## Run 212 stress 0.2295741 
## Run 213 stress 0.2061992 
## Run 214 stress 0.2247632 
## Run 215 stress 0.2260427 
## Run 216 stress 0.2250272 
## Run 217 stress 0.2267634 
## Run 218 stress 0.2115603 
## Run 219 stress 0.2125271 
## Run 220 stress 0.2123807 
## Run 221 stress 0.211908 
## Run 222 stress 0.2105442 
## Run 223 stress 0.2123835 
## Run 224 stress 0.2115023 
## Run 225 stress 0.2264694 
## Run 226 stress 0.2054933 
## ... New best solution
## ... Procrustes: rmse 0.001084309  max resid 0.004229843 
## ... Similar to previous best
## *** Solution reached
```

```r
# you can also try plot(data_nmds) at this point
data_nmds_results <- export_nmds(transects, data_nmds, score_name = 'Transect_id')
```

```
## Joining, by = "Transect_id"
```

To look at the results we plot the exported data and color by habitat type:


```r
p1 <- ggplot(data_nmds_results$spots,
       aes(x = NMDS1, y = NMDS2,
           group = habitat_type)) + 
  # nMDS needs equal axis
  coord_equal() +
  # add lines for x = 0 and y = 0
  geom_hline(yintercept = 0, color = 'lightgray', alpha = .4) +
  geom_vline(xintercept = 0, color = 'lightgray', alpha = .4) +
  # add indication of species influence
  geom_text(inherit.aes = FALSE,
            data = data_nmds_results$species,
            aes(x = NMDS1, y = NMDS2, label = Species)) +
  # add scoring of the transects
  geom_point(aes(colour = habitat_type), size = 2) +
  # add the outer hull for each habitat type
  geom_encircle(aes(colour = habitat_type, fill = habitat_type), 
                s_shape = 1, alpha = 0.2, size = 1, expand = 0)+
  # add a little space on either side of the plot (for text)
  scale_x_continuous(expand = c(.01,.05))+
  # adjust plot layout
  theme(legend.position = 'bottom',
        panel.grid = element_blank())

# print plot
p1
```


```
## Warning: Removed 1 rows containing missing values (geom_text).
```

<img src="nmds_files/figure-html/unnamed-chunk-12-1.png" width="1008" />

Alternatively we can highlight the differnt sites using color.

Note that all that has changed here is several replacements from `colour = habitat_type` to `colour = Site` - same for `fill =`.


```r
p2 <- ggplot(data_nmds_results$spots,
       aes(x = NMDS1, y = NMDS2,
           group = Site)) + 
  # nMDS needs equal axis
  coord_equal() +
  # add lines for x = 0 and y = 0
  geom_hline(yintercept = 0, color = 'lightgray', alpha = .4) +
  geom_vline(xintercept = 0, color = 'lightgray', alpha = .4) +
  # add indication of species influence
  geom_text(inherit.aes = FALSE,
            data = data_nmds_results$species,
            aes(x = NMDS1, y = NMDS2, label = Species)) +
  # add scoring of the transects
  geom_point(aes(colour = Site), size = 2) +
  # add the outer hull for each site
  geom_encircle(aes(colour = Site, fill = Site),
                s_shape = 1, alpha = 0.2, size = 1, expand = 0)+
  # add a little space on either side of the plot (for text)
  scale_x_continuous(expand = c(.01,.05))+
  # adjust plot layout
  theme(legend.position = 'bottom',
        panel.grid = element_blank())

# print plot
p2
```


```
## Warning: Removed 1 rows containing missing values (geom_text).
```

<img src="nmds_files/figure-html/unnamed-chunk-14-1.png" width="1008" />

we can export these plots using `ggsave()`.


```r
ggsave(filename = 'nmds_by_haitat_type.pdf', plot = p1, width = 9, height = 9)
ggsave(filename = 'nmds_by_site.pdf', plot = p2, width = 9, height = 9)
```

## Permanova

We can gain som insight from the nMDS plots, but if we want to check if the clusters on the nMDS are actually significantly differnt, we need to run a Permanova.

The Permanova is also based on ecological differnces between the transects (just like nMDS), but here we need to compute the ourselves (nMDS does this step internally).


```r
# Permanova ------------------
# compute ecological distances between the transects
distances <- vegdist(data_matrix, 
                     distance = 'bray')

# export the hapitat types according to order of the distance matrix
distances_groups <- transects[match(labels(distances),
                                       transects$Transect_id),]$habitat_type
```

One point we unfortunately could not discuss in class (because of time 😿), is the fact that like many other statistical tests, Permanova also has some assumptions about our data.

So, before we start, we need to check if the assumptions for permanova are met by our data [@Anderson01]:

<center> *"The only assumption of the test is that [...] the observations are independent and that they have similar distributions"*</center>

This means that our data should be *homogeneous dispersed* - our clusters should be of  approximately equal size (similar to equal variances when running an ANOVA or a t-test).
The idea here is that if the clusters are of differnt sizes, the permanova might give a significant result because of the size diffence and not because of a difference in location of the clusters on the nMDS.

To check for homogeneous dispersion, we can either use a statistical test (permutation test) or visualize the dispersion.


```r
# For the test for homogeneous dispersion we nee to compute the
# distance of each transect to the center of its respective cluster
# (s. https://k-hench.github.io/nmds/nmds.html#9 )
distances_betadispersion <- betadisper(distances, 
                                       distances_groups)

# test homogenous dispersion
permutest(distances_betadispersion) 
```

```
## 
## Permutation test for homogeneity of multivariate dispersions
## Permutation: free
## Number of permutations: 999
## 
## Response: Distances
##           Df Sum Sq  Mean Sq      F N.Perm Pr(>F)  
## Groups     1 0.1401 0.140096 4.0667    999  0.045 *
## Residuals 44 1.5158 0.034449                       
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Here, we hope for a p-value larger than 0.05, since a value smaller than 0.05 would tell us that the dispersions of our clusters are indeed differnt (that the clusters are of different size).

We can also inspect the distances of all transects to the center of their respective cluster visually for the same inforamtion:


```r
tibble(distances = distances_betadispersion$distances,
       habitat_type = distances_betadispersion$group) %>%
  ggplot(aes(x = habitat_type, y = distances, fill = habitat_type))+
  geom_boxplot()
```

<img src="nmds_files/figure-html/unnamed-chunk-19-1.png" width="672" />

At the stage when I checked the data it looked like our data is very borderline - actually the clusters for mangrove and reef seem to be of different sizes.

Fortunately the author of permanova basically told people to chill in such situations as long as their experimantal design was balanced [@Anderson13]:

<center> *"In contrast, PERMANOVA and Pillai’s trace were largely unaffected by heterogeneity for balanced designs"*</center>

So, lets quickly doublecheck:


```r
transects %>%
  group_by(habitat_type) %>%
  count()
```

```
## # A tibble: 2 x 2
## # Groups:   habitat_type [2]
##   habitat_type     n
##   <chr>        <int>
## 1 Mangrove        22
## 2 Reef            24
```

To me 22/24 seems reasonably balanced, so we finally run the permanova: 


```r
data_permanova <- adonis(formula = distances ~ distances_groups,
                         permutations = 999, method = 'bray')

print(data_permanova) 
```

```
## 
## Call:
## adonis(formula = distances ~ distances_groups, permutations = 999,      method = "bray") 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## distances_groups  1    3.9133  3.9133  14.904 0.25302  0.001 ***
## Residuals        44   11.5528  0.2626         0.74698           
## Total            45   15.4660                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

At this point it is up to you to interpret the results - as far as stats go we're done.
