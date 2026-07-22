# ==========================================
# 🧬 PIPELINE: CLINICAL VARIANT TAGGER PROJECT
# Author: Alisha Namakula
# ==========================================

library(dplyr)
library(ggplot2)

#Load Amino Acid Reference Sheet csv
aa_reference <- read.csv("amino_acids_reference_sheet.csv")

# 1. DATA CLEANING & REORDERING
clean_missense_data <- UMD_variants_EU[, c("WT AA_1", "Mutant AA_1", "Codon")]

clean_reference_data <- amino_acids_reference_sheet %>% 
  select(V1, V2, V3, V9, V10, V11, V12)

step_one_data <- left_join(
  clean_missense_data,
  clean_reference_data,
  by = c("WT AA_1" = "V2")
)

reorder_codon <- step_one_data %>% relocate(Codon, .before = `WT AA_1`)

step_two_data <- left_join(
  reorder_codon, 
  clean_reference_data,
  by = c("Mutant AA_1" = "V2")
)

reorder_mutant_final_data <- step_two_data %>% relocate(`Mutant AA_1`, .before = `V1.y`)

# 2. RENAMING COLUMNS & ENFORCING NUMERIC TYPES
Renamed_Final_Data <- reorder_mutant_final_data %>%
  rename(
    WT_Name       = V1.x, 
    WT_MolWt      = V3.x, 
    WT_Charge     = V9.x,    
    WT_Hydro_ph2  = V10.x, 
    WT_Hydro_ph7  = V11.x, 
    WT_Category   = V12.x,
    Mut_Name      = V1.y, 
    Mut_MolWt     = V3.y, 
    Mut_Charge    = V9.y,  
    Mut_Hydro_ph2 = V10.y, 
    Mut_Hydro_ph7 = V11.y, 
    Mut_Category  = V12.y
  ) %>%
  mutate(
    WT_Charge     = as.numeric(WT_Charge),
    Mut_Charge    = as.numeric(Mut_Charge),
    WT_Hydro_ph2  = as.numeric(WT_Hydro_ph2),
    Mut_Hydro_ph2 = as.numeric(Mut_Hydro_ph2),
    WT_MolWt      = as.numeric(WT_MolWt),
    Mut_MolWt     = as.numeric(Mut_MolWt),
    WT_Hydro_ph7  = as.numeric(WT_Hydro_ph7),
    Mut_Hydro_ph7 = as.numeric(Mut_Hydro_ph7)
  )

# 3. OTHER FEATURES & CLINICAL TAGGER
Switch_Final_Data <- Renamed_Final_Data %>%
  mutate(Mutation_Route = paste0(`WT AA_1`, "->", `Mutant AA_1`))

Tagged_Data <- Switch_Final_Data %>%
  mutate(
    Mol_Wt_Pct_change = ((Mut_MolWt - WT_MolWt) / WT_MolWt) * 100,
    Hydro_ph7_Pct_change = ((Mut_Hydro_ph7 - WT_Hydro_ph7) / WT_Hydro_ph7) * 100,
    Clinical_status = case_when(
      is.na(Mol_Wt_Pct_change) | is.na(Hydro_ph7_Pct_change) ~ "CHECK FOR FRAMESHIFT",
      WT_Charge != Mut_Charge ~ "PATHOGENIC",
      abs(Mol_Wt_Pct_change) > 30 | abs(Hydro_ph7_Pct_change) > 30  ~ "PATHOGENIC",
      TRUE ~ "BENIGN"
    )
  )

reordered_Final_Tagged_Data <- Tagged_Data %>% relocate(Mutation_Route, .before = `WT AA_1`)

# 4. PIPELINE FILTERS

#See all data without missing values (from frameshifts)

Filter_1_data <- Tagged_Data %>%
  filter(!is.na(Mol_Wt_Pct_change)) %>%
  filter(!is.na(Hydro_ph7_Pct_change))

#See only pathogenic data
Filter_2_data <- Tagged_Data %>%
  filter(Clinical_status == "PATHOGENIC")

-
Removed_pct_columns_from_data <- Tagged_Data %>%
  select(-c(Mol_Wt_Pct_change, Hydro_ph7_Pct_change))

Final_Viewing_Data <- Removed_pct_columns_from_data

# ==========================================
# 5. DATA VISUALIZATION SECTION
# ==========================================

# Graph 1: Bar Chart
ggplot(Tagged_Data, aes(x = Clinical_status, fill = Clinical_status)) +
  geom_bar(color = "white", alpha = 0.85, width = 0.6) +
  scale_fill_manual(values = c(
    "BENIGN" = "#4FBF67", 
    "PATHOGENIC" = "#FF5A5F", 
    "CHECK FOR FRAMESHIFT" = "#FFB400"
  )) +
  labs(
    title = "Distribution of Pathogenicity Across Patient Missense Mutations",
    x = "Clinical Tag",
    y = "Missense Mutations (Count)",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12))

# Graph 2a: Histogram (Hydrophobicity)
ggplot(Tagged_Data, aes(x = Hydro_ph7_Pct_change)) +
  geom_histogram(fill = "#FF007F", color = "white", alpha = 0.8) +
  labs(
    title = "Normal Distribtuion of Hydrophobicity Changes Across Missense Mutations",
    x = "Percentage Change in Hydrophobicity at pH7 (%)",
    y = "Missense Mutations (Count)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12))

# Graph 2b: Histogram (Molecular Weight)
ggplot(Filter_1_data, aes(x = Mol_Wt_Pct_change)) +
  geom_histogram(binwidth = 8, fill = "#05C3DD", color = "white", alpha = 0.8) +
  labs(
    title = "Most Prevalent Molecular Weight Change Across Missense Mutations",
    x = "Percentage Change in Molecular Weight (%)",
    y = "Missense Mutations (Count)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12))

# Graph 3: Scatter Plot (Molecular Weight vs. Hydrophobicity)
ggplot(Filter_1_data, aes(x = Mol_Wt_Pct_change, y = Hydro_ph7_Pct_change, color = Clinical_status)) +
  geom_point(alpha = 0.67, size = 2.5) +
  scale_color_manual(values = c(
    "BENIGN" = "#4FBF67", 
    "PATHOGENIC" = "#FF5A5F", 
    "CHECK FOR FRAMESHIFT" = "#FFB400"
  )) +
  labs(
    title = "Molecular Weight and Hydrophobicity Percentage Changes in Missense Mutations",
    subtitle = "Sub-Contributing Factors to Pathogenicity",
    x = "Molecular Weight Change (%)",
    y = "Hydrophobicity Change at pH7 (%)",
    color = "Clinical Status"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12))