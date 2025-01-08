library(ggplot2)
library(tidyverse) 
library(dplyr)

# create data for world coordinates using map_data() function 
world_coordinates <- map_data("world") 

# read in your e coli data:
data = read_csv("P_aeruginosa_genome_ids.csv")
# Remove rows where isolation_country is 0
data <- data %>% 
  filter(Isolation_country != "0") 

# Count the number of genomes per country
country_counts <- data %>% 
  count(Isolation_country, name = "count")


world_map <- map_data("world")

# Rename countries for consistency if needed
# Example: country_counts$isolation_country <- recode(country_counts$isolation_country, "USA" = "United States")

# Merge data with world map
world_data <- world_map %>%
  left_join(country_counts, by = c("region" = "Isolation_country"))


ggplot() +
  geom_polygon(
    data = world_data,
    aes(x = long, y = lat, group = group, fill = count),
    color = "white", size = 0.2
  ) +
  scale_fill_gradient(low = "skyblue", high = "navyblue", na.value = "gray") +
  theme_minimal() +
  labs(
    title = "World Map of P. aeruginosa genomes",
    fill = "Count"
  )


world_data_reduced <- subset(world_data, count > 0)
