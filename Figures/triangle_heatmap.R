# Load the packages
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)



# Assume your matrix is stored in a data frame called 'my_matrix'
my_matrix <- read.csv("matrix_heatmap_input.csv", header = TRUE, check.names = FALSE, row.names = 1)

# Convert to matrix and melt into long format
my_matrix <- as.matrix(my_matrix)
long_matrix <- melt(my_matrix, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels
long_matrix$Gene1 <- factor(long_matrix$Gene1, levels = rownames(my_matrix))
long_matrix$Gene2 <- factor(long_matrix$Gene2, levels = colnames(my_matrix))

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COG0789", "COG2602")

# Create the heatmap using ggplot2
ggplot(long_matrix, aes(Gene2, Gene1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "black",  
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Interaction") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "Gene Interaction Heatmap")



##### If the same have a lighter colour 

# Read the matrix and convert it to matrix format
my_matrix <- read.csv("matrix_heatmap_input.csv", header = TRUE, check.names = FALSE, row.names = 1)
my_matrix <- as.matrix(my_matrix)

# Melt the matrix into a long format for ggplot
long_matrix <- melt(my_matrix, varnames = c("Gene1", "Gene2"), na.rm = TRUE)
#write.csv(long_matrix_ecoli, '/Users/lucydillon/Library/CloudStorage/OneDrive-Personal/Documents/Postdoc/Pangenome_MDR_project/long_format.csv', row.names=FALSE)

# Factorize the gene labels to preserve order
long_matrix$Gene1 <- factor(long_matrix$Gene1, levels = rownames(my_matrix))
long_matrix$Gene2 <- factor(long_matrix$Gene2, levels = colnames(my_matrix))

# Create a new column to classify interaction type
long_matrix$interaction_type <- apply(long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(my_matrix[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
long_matrix <- long_matrix[!is.na(long_matrix$interaction_type), ]

triangle_data <- data.frame(
  Gene1 = c("COG243", "COG449", "COG123"),  # Example genes
  Gene2 = c("COG693", "COG859", "COG1414"),
  value = c(1, -1, 1)  # Use 1 for red, -1 for blue
)

triangle_matrix = read.csv("All_coinfinder_interactions_data_points.csv",header = TRUE, check.names = FALSE, row.names = 1)
triangle_matrix <- as.matrix(triangle_matrix)
triangle_data <- melt(triangle_matrix, varnames = c("Gene1", "Gene2"), na.rm = TRUE)


# Factorize the gene labels in the triangle data to match the main matrix
triangle_data$Gene1 <- factor(triangle_data$Gene1, levels = rownames(my_matrix))
triangle_data$Gene2 <- factor(triangle_data$Gene2, levels = colnames(my_matrix))
# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COG0789", "COG2602")

# Create the heatmap with ggplot2
all_interactions_plot <- ggplot(long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) +
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species")) +
        #geom_point(data = subset(triangle_data, value == "1"), aes(Gene2, Gene1), 
         #    shape = 24, size = 1, color = "red", fill = "red") +  # Red triangle for 1
  
      #geom_point(data = subset(triangle_data, value == "-1"), aes(Gene2, Gene1), 
             #shape = 25, size = 1, color = "blue", fill = "blue") +  # Blue triangle for -1
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "A. All gene interactions. ")

# Amikacin_cefepime

# Assuming your matrix is called 'my_matrix'
selected_rows <- c(1, 2, 3, 12, 14, 17, 28)
selected_cols <- c(1, 2, 3, 12, 14, 17, 28)

# Selecting the specific rows and columns
amikacin_cefepime <- my_matrix[selected_rows, selected_cols]

# Melt the matrix into a long format for ggplot
AC_long_matrix <- melt(amikacin_cefepime, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels to preserve order
AC_long_matrix$Gene1 <- factor(AC_long_matrix$Gene1, levels = rownames(amikacin_cefepime))
AC_long_matrix$Gene2 <- factor(AC_long_matrix$Gene2, levels = colnames(amikacin_cefepime))

# Create a new column to classify interaction type
AC_long_matrix$interaction_type <- apply(AC_long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(amikacin_cefepime[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
AC_long_matrix <- AC_long_matrix[!is.na(AC_long_matrix$interaction_type), ]

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COGNA789", "COG26NA2")

# Create the heatmap with ggplot2
amikacin_cefepime_plot <- ggplot(AC_long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(AC_long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  # Improve aesthetics
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(AC_long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(AC_long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(AC_long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(AC_long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) + 
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species"),
                    guide="none") +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "B. Amikacin cefepime model gene interactions")


# ciprofloxacin_ertapenem 
# Assuming your matrix is called 'my_matrix'
selection_CE <- c(18,21,13)

# Selecting the specific rows and columns
ciprofloxacin_ertapenem <- my_matrix[selection_CE, selection_CE]

# Melt the matrix into a long format for ggplot
CE_long_matrix <- melt(ciprofloxacin_ertapenem, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels to preserve order
CE_long_matrix$Gene1 <- factor(CE_long_matrix$Gene1, levels = rownames(ciprofloxacin_ertapenem))
CE_long_matrix$Gene2 <- factor(CE_long_matrix$Gene2, levels = colnames(ciprofloxacin_ertapenem))

# Create a new column to classify interaction type
CE_long_matrix$interaction_type <- apply(CE_long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(ciprofloxacin_ertapenem[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
CE_long_matrix <- CE_long_matrix[!is.na(CE_long_matrix$interaction_type), ]

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COGNA789", "COG26NA2")

# Create the heatmap with ggplot2
ggplot(CE_long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "red", 
                               "ecoli_minus1_pseudo_plus1" = "blue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(CE_long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  # Improve aesthetics
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(CE_long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(CE_long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(CE_long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(CE_long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) +
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species")) +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "ciprofloxacin_ertapenem ")

# gentamicin_levofloxacin 
selection_GL <- c(25,13,19)

# Selecting the specific rows and columns
gentamicin_levofloxacin <- my_matrix[selection_CE, selection_CE]

# Melt the matrix into a long format for ggplot
GL_long_matrix <- melt(gentamicin_levofloxacin, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels to preserve order
GL_long_matrix$Gene1 <- factor(GL_long_matrix$Gene1, levels = rownames(gentamicin_levofloxacin))
GL_long_matrix$Gene2 <- factor(GL_long_matrix$Gene2, levels = colnames(gentamicin_levofloxacin))

# Create a new column to classify interaction type
GL_long_matrix$interaction_type <- apply(GL_long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(gentamicin_levofloxacin[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
GL_long_matrix <- GL_long_matrix[!is.na(GL_long_matrix$interaction_type), ]

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COGNA789", "COG26NA2")

# Create the heatmap with ggplot2
ggplot(GL_long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "red", 
                               "ecoli_minus1_pseudo_plus1" = "blue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(GL_long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  # Improve aesthetics
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(GL_long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(GL_long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(GL_long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(GL_long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) + 
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species")) +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "gentamicin_levofloxacin ")

# ceftriaxone_ciprofloxacin 

selection_CC <- c(19, 20, 31)

# Selecting the specific rows and columns
ceftriaxone_ciprofloxacin <- my_matrix[selection_CC, selection_CC]

# Melt the matrix into a long format for ggplot
CC_long_matrix <- melt(ceftriaxone_ciprofloxacin, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels to preserve order
CC_long_matrix$Gene1 <- factor(CC_long_matrix$Gene1, levels = rownames(ceftriaxone_ciprofloxacin))
CC_long_matrix$Gene2 <- factor(CC_long_matrix$Gene2, levels = colnames(ceftriaxone_ciprofloxacin))

# Create a new column to classify interaction type
CC_long_matrix$interaction_type <- apply(CC_long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(ceftriaxone_ciprofloxacin[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
CC_long_matrix <- CC_long_matrix[!is.na(CC_long_matrix$interaction_type), ]

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COGNA789", "COG26NA2")

# Create the heatmap with ggplot2
ggplot(CC_long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "red", 
                               "ecoli_minus1_pseudo_plus1" = "blue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(CC_long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  # Improve aesthetics
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(CC_long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(CC_long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(CC_long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(CC_long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) + 
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species")) +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "ceftriaxone_ciprofloxacin")

# ampicillin_levofloxacin  

selection_AL <- c(6, 7, 15)

# Selecting the specific rows and columns
ampicillin_levofloxacin <- my_matrix[selection_AL, selection_AL]

# Melt the matrix into a long format for ggplot
AL_long_matrix <- melt(ampicillin_levofloxacin, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels to preserve order
AL_long_matrix$Gene1 <- factor(AL_long_matrix$Gene1, levels = rownames(ampicillin_levofloxacin))
AL_long_matrix$Gene2 <- factor(AL_long_matrix$Gene2, levels = colnames(ampicillin_levofloxacin))

# Create a new column to classify interaction type
AL_long_matrix$interaction_type <- apply(AL_long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(ampicillin_levofloxacin[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
AL_long_matrix <- AL_long_matrix[!is.na(AL_long_matrix$interaction_type), ]

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COG0789", "COG2602")

# Create the heatmap with ggplot2
ampicillin_levofloxacin_plot <- ggplot(AL_long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "red", 
                               "ecoli_minus1_pseudo_plus1" = "blue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(AL_long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  # Improve aesthetics
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(AL_long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(AL_long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(AL_long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(AL_long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) + 
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species"),
                    guide="none") +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "D. Ampicillin Levofloxacin \nmodel gene interactions")



# cefepime_levofloxacin 

selection_CL <- c(16, 23, 30)

# Selecting the specific rows and columns
cefepime_levofloxacin  <- my_matrix[selection_CL, selection_CL]

# Melt the matrix into a long format for ggplot
CL_long_matrix <- melt(cefepime_levofloxacin, varnames = c("Gene1", "Gene2"), na.rm = TRUE)

# Factorize the gene labels to preserve order
CL_long_matrix$Gene1 <- factor(CL_long_matrix$Gene1, levels = rownames(cefepime_levofloxacin))
CL_long_matrix$Gene2 <- factor(CL_long_matrix$Gene2, levels = colnames(cefepime_levofloxacin))

# Create a new column to classify interaction type
CL_long_matrix$interaction_type <- apply(CL_long_matrix, 1, function(row) {
  gene1 <- as.character(row["Gene1"])
  gene2 <- as.character(row["Gene2"])
  value <- as.numeric(row["value"])
  
  # Check for the reverse interaction (E. coli vs P. aeruginosa)
  reverse_value <- as.numeric(cefepime_levofloxacin[gene2, gene1])
  
  # Classify based on the interaction type
  if (is.na(value) || is.na(reverse_value)) {
    return(NA)  # Handle missing values
  } else if (value == 0 || reverse_value == 0) {
    return("none")  # For any 0 interaction
  } else if (value == -1 && reverse_value == -1) {
    return("both_minus1")  # Both interactions are -1
  } else if (value == 1 && reverse_value == 1) {
    return("both_plus1")  # Both interactions are 1
  } else if (value == 1 && reverse_value == -1) {
    return("ecoli_plus1_pseudo_minus1")  # E. coli is 1, P. aeruginosa is -1
  } else if (value == -1 && reverse_value == 1) {
    return("ecoli_minus1_pseudo_plus1")  # E. coli is -1, P. aeruginosa is 1
  }
})

# Remove rows where interaction_type is NA (for cases where either side is NA)
CL_long_matrix <- CL_long_matrix[!is.na(CL_long_matrix$interaction_type), ]

# List of specific gene labels to highlight
highlight_genes <- c("COG2271", "COG0789", "COG2602")

# Create the heatmap with ggplot2
cefepime_levofloxacin_plot <- ggplot(CL_long_matrix, aes(Gene2, Gene1, fill = interaction_type)) +
  geom_tile(color = "white") +
  
  # Manually define colors for interaction types
  scale_fill_manual(values = c("both_minus1" = "lightblue", 
                               "both_plus1" = "lightcoral", 
                               "ecoli_plus1_pseudo_minus1" = "red", 
                               "ecoli_minus1_pseudo_plus1" = "blue", 
                               "none" = "black"),
                    name = "Interaction") +
  geom_point(data = subset(CL_long_matrix, interaction_type == "both_plus1"), 
             aes(Gene2, Gene1), color = "black", shape = 19, size = 2) +
  
  # Improve aesthetics
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1,
                                   face = ifelse(levels(CL_long_matrix$Gene2) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(CL_long_matrix$Gene2) %in% highlight_genes, "#8338ec", "black")),
        axis.text.y = element_text(size = 12,
                                   face = ifelse(levels(CL_long_matrix$Gene1) %in% highlight_genes, "bold", "plain"),
                                   color = ifelse(levels(CL_long_matrix$Gene1) %in% highlight_genes, "#8338ec", "black"))) + 
  scale_fill_manual(values = c("both_plus1" = "lightcoral", 
                               "ecoli_minus1_pseudo_plus1" = "lightblue", 
                               "ecoli_plus1_pseudo_minus1" = "lightcoral"),name = "Coinfinder interaction",  # Legend title
                    labels = c("Associated in both species",
                               "Disassociated in one species",
                               "Associated in one Species"),
                    guide="none") +
  labs(x = expression( ~ italic("P. aeruginosa")), 
       y = expression( ~ italic("E. coli")), 
       title = "C. Cefepime Levofloxacin \nmodel gene interactions")


gap <- textGrob("") 

grid.arrange(
  all_interactions_plot, amikacin_cefepime_plot, cefepime_levofloxacin_plot,ampicillin_levofloxacin_plot, gap,
  layout_matrix = rbind(
    c(1, 1, 1),
    c(1, 1, 1),
    c(1, 1, 1),
    c(2, 2,2,2,2,5, 3,3,3,5),
    c(2, 2,2,2,2,5, 4,4,4,5)
  )
)

