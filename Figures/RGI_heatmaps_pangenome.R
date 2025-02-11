# Load necessary libraries
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggforce)
library(viridis)
library(ggpubr)
library(extrafont)
font_import()
loadfonts(device = "pdf")  # Ensures fonts are available for PDF export

# Example dataset
E_coli_summary <- data.frame(
  Category = c("Core", "Softcore", "Shell", "Cloud"),
  Value = c(649, 1205, 4623, 132772),  # Values to be displayed
  Radius = c(1, 2, 3, 4)  # Arbitrary increasing radius for visualization
)



view(E_coli_pangenome)

P_aeruginosa_summary <- data.frame(
  Category = c("Core", "Softcore", "Shell", "Cloud"),
  Value = c(841, 2059, 4206, 163060),  # Values to be displayed
  Radius = c(1, 2, 3, 4)  # Arbitrary increasing radius for visualization
)






df$Species <- "P. aeruginosa"
EC$Species <- "E. coli"
combined_df <- rbind(df, EC)


print(P)


PA <- read.csv("Pseudomonas_core_rgi_genes_reduced_names.csv")
EC_core <- read.csv("E_coli_core_rgi_genes_reduced_names.csv")
combined_core_df <- rbind(PA, EC_core)


# Plot nested circles
E_coli_summary$Category <- factor(E_coli_summary$Category, 
                                        levels = c("Cloud","Shell","Softcore","Core"))  # Desired order
E_coli_pangenome <- ggplot() +
  geom_circle(data = E_coli_summary, aes(x0 = 0, y0 = 0, r = Radius, fill = Category), alpha = 0.4) + 
  geom_text(data = E_coli_summary, aes(x = 0, y = -Radius + 0.6, label = paste(Value)), 
            size = 4, color = "black") +  
  coord_fixed() +  
  theme_void() +  
  scale_fill_manual(values = c("Core" = "#6247aa", "Softcore" = "#9163cb", 
                               "Shell" = "#c19ee0", "Cloud" = "#dec9e9"),
                    breaks = c("Cloud","Shell","Softcore","Core")) +
  labs(title = "A.") +  # Main title
  annotate("text", x = -2, y = 5, label = expression(i. ~ italic("E. coli")), size = 4, fontface = "bold") # Small "i"

P_aeruginosa_summary$Category <- factor(P_aeruginosa_summary$Category, 
                                        levels = c("Cloud","Shell","Softcore","Core"))  # Desired order

P_aeruginosa_pangenome <- ggplot() +
  geom_circle(data = P_aeruginosa_summary, aes(x0 = 0, y0 = 0, r = Radius, fill = Category), alpha = 0.4) + 
  geom_text(data = P_aeruginosa_summary, aes(x = 0, y = -Radius + 0.6, label = paste(Value)), 
            size = 4, color = "black") +  
  coord_fixed() +  
  theme_void() +  
  scale_fill_manual(values = c("Core" = "#054a91", "Softcore" = "#3e7cb1", 
                               "Shell" = "#81a4cd", "Cloud" = "#dbe4ee"), 
                    breaks = c("Cloud","Shell","Softcore","Core")) + 
  labs(title = "") +  
  annotate("text", x = -2, y = 5, label = expression(ii. ~ italic("P. aeruginosa")), size = 4, fontface = "bold")

# plot all AMR heatmaps
combined_df$Species <- factor(combined_df$Species, 
                              levels = c("E. coli", "P. aeruginosa"), 
                              labels = c("i. E. coli", "ii. P. aeruginosa"))

All_AMR <- ggplot(combined_df, aes(x = Best_Hit_ARO, y = DrugClass, fill = Percentage)) +
  geom_tile() +
  scale_fill_viridis(breaks = c(0, 15, 30, 45, 60, 75, 90)) +
  facet_wrap(~Species, ncol = 1, scales = "free_x") +  # Labels are already formatted
  labs(title = "B.", y = "Drug Class", x = "Gene") +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1),
        legend.position = "right",
        strip.text = element_text(face = "italic", color = "black"), 
        strip.background = element_rect(fill = "grey80", color = "grey80")) 

# plot core amr heatmaps

combined_core_df$Species <- factor(combined_core_df$Species, 
                              levels = c("E. coli", "P. aeruginosa"), 
                              labels = c("i. E. coli", "ii. P. aeruginosa"))

core_AMR <- ggplot(combined_core_df, aes(x = Best_Hit_ARO, y = DrugClass, fill = Percentage)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", breaks = c(95, 96, 97, 98, 99)) +
  facet_wrap(~Species, ncol = 1, scales = "free_x", 
             labeller = as_labeller(c("P. aeruginosa" = expression(italic("P. aeruginosa")), 
                                      "E. coli" = expression(italic("E. coli"))))) +
  labs(title = "C.", y = "Drug Class", x = "Gene") +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(hjust = 1),
        legend.position = "right",
        strip.text = element_text(face = "italic", color = "black"), # Make text readable
        strip.background = element_rect(fill = "grey80", color = "grey80")) # Add grey backgroun



grid.arrange(
  E_coli_pangenome, P_aeruginosa_pangenome, core_AMR, All_AMR,
  layout_matrix = rbind(
    c(1, 3, 3),
    c(2, 3, 3),
    c(4, 4, 4),
    c(4, 4, 4)
  )
)

