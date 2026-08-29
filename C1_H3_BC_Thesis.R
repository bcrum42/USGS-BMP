# ============================================================
# SCRIPT: H3_NMDS_Biomass.R
# PURPOSE: NMDS ordination of macroinvertebrate biomass to
#          test whether assemblage composition differs more
#          between habitat types than among streams
# H3: Biomass-weighted taxonomic identity will be more
#     variable between depositional habitats and riffles
#     than among streams
# Last updated: 2026-04-15
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(stringr)
library(ggrepel)
library(indicspecies)
library(flextable)
library(officer)
library(cowplot)


# ============================================================
# SECTION 1: LOAD DATA
# ============================================================

all_biomass_bytaxon <- read.csv("outputs/rep_biomass_by_taxon_4.7.26.csv",
                                stringsAsFactors = FALSE)



dep_rich <- read.csv("outputs/taxa_lists/dep_taxa.csv")

rif.rich <- read.csv("outputs/taxa_lists/rif_taxa.csv")

#build lookup sheet
taxon_lookup <- bind_rows(dep_rich, rif.rich) %>%
  select(Genus, Order, Family) %>%   # adjust column names to match yours
  distinct()

# Define once, use in both table builds
order_levels <- c(
  "Ephemeroptera", "Odonata", "Plecoptera", "Hemiptera",
  "Coleoptera", "Megaloptera", "Trichoptera", "Lepidoptera",
  "Diptera"
)

# ============================================================
# SECTION 2: HTYPE CORRECTIONS
# Corrects misclassified habitat types identified during QA.
# dep_biomass and rif_biomass are split first to apply
# depositional-only corrections, then recombined into all_raw.
# ============================================================

dep_biomass <- all_biomass_bytaxon %>% filter(htype != "riffle")
rif_biomass <- all_biomass_bytaxon %>% filter(htype == "riffle")

dep_biomass <- dep_biomass %>%
  filter(!(Site == "SPOU" & Sample == 4)) %>%
  mutate(htype = case_when(
    Site == "BEAV" & Sample == 2 & htype == "rts" ~ "ub",
    Site == "BEAV" & Sample == 4 & htype == "rts" ~ "ub",
    Site == "CEDR" & Sample == 2 & htype == "rts" ~ "sng",
    Site == "EIDS" & Sample == 1 & htype == "rts" ~ "ub",
    Site == "EIDS" & Sample == 2 & htype == "sng" ~ "ub",
    Site == "LONG" & Sample == 3 & htype == "sng" ~ "da",
    Site == "MEAD" & Sample == 1 & htype == "rts" ~ "ub",
    Site == "MEAD" & Sample == 4 & htype == "sng" ~ "da",
    Site == "MLPG" & Sample == 2 & htype == "sng" ~ "da",
    Site == "MOSS" & Sample == 3 & htype == "sng" ~ "da",
    Site == "MOSS" & Sample == 5 & htype == "sng" ~ "da",
    Site == "MUDD" & Sample == 4 & htype == "sng" ~ "da",
    Site == "NPAS" & Sample == 4 & htype == "sng" ~ "da",
    Site == "NPAS" & Sample == 5 & htype == "sng" ~ "da",
    Site == "PINE" & Sample == 2 & htype == "sng" ~ "da",
    Site == "PINE" & Sample == 3 & htype == "sng" ~ "da",
    Site == "PINE" & Sample == 4 & htype == "sng" ~ "da",
    Site == "PINE" & Sample == 5 & htype == "sng" ~ "da",
    Site == "PISG" & Sample == 3 & htype == "rts" ~ "ub",
    Site == "PISG" & Sample == 5 & htype == "sng" ~ "ub",
    Site == "POAG" & Sample == 2 & htype == "sng" ~ "da",
    Site == "POAG" & Sample == 5 & htype == "sng" ~ "da",
    Site == "SPOU" & Sample == 2 & htype == "sng" ~ "da",
    Site == "THOR" & Sample == 4 & htype == "ub"  ~ "rts",
    Site == "TOMS" & Sample == 5 & htype == "sng" ~ "da",
    TRUE ~ htype
  ))

#remove snag as of 8.26.2026 
dep_biomass <- dep_biomass %>% 
  mutate(htype = if_else(
    htype == "sng", "da", htype
  ))

# Recombine and add binary htype2 grouping (depositional vs riffle)
depositional_types <- c("rts", "veg", "ub", "da")

all_raw <- bind_rows(
  dep_biomass %>% mutate(htype = tolower(htype)),
  rif_biomass  %>% mutate(htype = "riffle")
) %>%
  mutate(
    htype2     = ifelse(htype == "riffle", "riffle", "depositional"),
    TAXANAME1  = ifelse(TAXANAME1 == "physa", "physidae", TAXANAME1)
  )


# ============================================================
# SECTION 3: BUILD COMMUNITY MATRIX
# Rows = individual habitat replicates (Site × Sample × htype)
# Columns = taxa (TAXANAME1)
# Values = raw Biomass_m2
# ============================================================



# Remove rows with missing biomass or taxon identity
biomass_clean <- all_raw %>%
  filter(!is.na(Biomass_m2), !is.na(TAXANAME1)) %>%
  mutate(rep_id = paste(Site, Sample, htype, sep = "_"))

# Pivot to wide format (one row per replicate, one column per taxon)
biomass_wide <- biomass_clean %>%
  select(rep_id, Site, Sample, htype2, htype, TAXANAME1, Biomass_m2) %>%
  pivot_wider(
    names_from  = TAXANAME1,
    values_from = Biomass_m2,
    values_fn   = sum,
    values_fill = 0
  )

# Separate metadata from the community matrix
meta        <- biomass_wide %>% select(rep_id, Site, Sample, htype, htype2)
comm_matrix <- biomass_wide %>%
  select(-rep_id, -Site, -Sample, -htype, -htype2) %>%
  as.matrix()
rownames(comm_matrix) <- meta$rep_id

# Drop all-zero rows (replicates with no biomass) and columns (taxa never recorded)
# These cause division-by-zero errors in Wisconsin standardization
zero_rows <- rowSums(comm_matrix) == 0
if (any(zero_rows)) {
  cat("Removing", sum(zero_rows), "all-zero replicates:\n")
  print(rownames(comm_matrix)[zero_rows])
  comm_matrix <- comm_matrix[!zero_rows, ]
  meta        <- meta[!zero_rows, ]
}

zero_cols <- colSums(comm_matrix) == 0
if (any(zero_cols)) {
  cat("Removing", sum(zero_cols), "all-zero taxa\n")
  comm_matrix <- comm_matrix[, !zero_cols]
}

cat("Final matrix:", nrow(comm_matrix), "replicates ×", ncol(comm_matrix), "taxa\n")
meta %>% count(htype)
meta %>% count(htype2)


# ============================================================
# SECTION 4: TRANSFORM AND STANDARDIZE
# Step 1: log(x+1) — reduces the influence of highly abundant taxa
# Step 2: Wisconsin double standardization — scales each species
#         by its maximum, then each site by its total; standard
#         practice for biomass-based NMDS
# ============================================================

comm_scaled <- wisconsin(log1p(comm_matrix))


# ============================================================
# SECTION 5: BRAY-CURTIS DISTANCE AND NMDS ORDINATION
# Bray-Curtis is appropriate for abundance/biomass data.
# k = 2 dimensions; stress < 0.10 = good, < 0.20 = acceptable.
# ============================================================

bc_dist <- vegdist(comm_scaled, method = "bray")

plot(bc_dist)

set.seed(123)
nmds <- metaMDS(
  comm_scaled,
  distance = "bray",
  k        = 2,
  trymax   = 100,
  trace    = FALSE
)

cat("NMDS stress:", nmds$stress, "\n")

stressplot(nmds)

# Per-axis R² (squared correlation between ordination distances and BC distances)
axis1_dist <- dist(scores(nmds, display = "sites")[, 1])
axis2_dist <- dist(scores(nmds, display = "sites")[, 2])
r2_axis1   <- round(cor(as.vector(axis1_dist), as.vector(bc_dist))^2, 3)
r2_axis2   <- round(cor(as.vector(axis2_dist), as.vector(bc_dist))^2, 3)
r2_nonmetric <- round(1 - nmds$stress^2, 3)   # non-metric R² from stress

cat("NMDS1 R²:", r2_axis1, "\n")
cat("NMDS2 R²:", r2_axis2, "\n")


# ============================================================
# SECTION 6: EXTRACT SCORES AND BUILD PLOTTING OBJECTS
# ============================================================

# Site scores joined to metadata
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  mutate(rep_id = rownames(.)) %>%
  left_join(meta, by = "rep_id")

write.csv(nmds_scores, "biomass_nmds_scores.csv")


# Centroids: htype2 (depositional vs riffle) and htype (6 types)
centroids_htype2 <- nmds_scores %>%
  group_by(htype2) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

centroids_htype <- nmds_scores %>%
  group_by(htype, htype2) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

# Caption label used across all figures
stress_label <- paste0("Stress = ", round(nmds$stress, 3),
                       "  |  Non-metric R² = ", r2_nonmetric)

# Shared display labels
htype_labels <- c(
  "da"           = "Debris accumulation",
  "rts"          = "Submerged roots",
 # "sng"          = "Snag",
  "ub"           = "Undercut bank",
  "veg"          = "Vegetation",
  "depositional" = "Non-riffle",
  "riffle"       = "Riffle"
)

# Shared display labels
htype_labels_long <- c(
  "da"           = "Debris accumulation",
  "rts"          = "Submerged roots",
 # "sng"          = "Snag",
  "ub"           = "Undercut bank",
  "veg"          = "Vegetation",
  "depositional" = "Non-riffle",
  "riffle"       = "Riffle"
)

# Color palettes
colors_htype2 <- c("depositional" = "#8B5E3C", "riffle" = "#84C441")

colors_htype <- c(
  "ub"     = "#c2b311",
  "da"     = "#4575B4",
  "rts"    = "#F97D0B",
  #"sng"    = "#D7191C",
  "veg"    = "#7B2D8B",
  "riffle" = "#4DAC26"
)


# ============================================================
# SECTION 7: INDICATOR SPECIES ANALYSIS
#
# Uses the IndVal.g method (Dufrene & Legendre 1997) via
# multipatt() from the indicspecies package.
#
# IndVal.g identifies taxa that are both highly faithful to
# (specificity) and highly abundant in (sensitivity) a given
# habitat group. The combined statistic ranges 0–1; higher = 
# stronger indicator. Significance is assessed by permutation
# (999 permutations, α = 0.05).
#
# Two analyses are run:
#   A) htype2 — depositional vs riffle (binary grouping)
#   B) htype  — all 6 habitat types
# ============================================================

# Shared exclusion list — ubiquitous taxa that don't aid habitat discrimination
excluded_taxa <- c("chironomidae", "gastropoda")

# -----------------------------------------------------------
# 7A: non-riffle vs RIFFLE (htype2)
# -----------------------------------------------------------
indval_htype2 <- multipatt(
  comm_scaled,
  meta$htype2,
  func    = "IndVal.g",
  control = how(nperm = 999)
)

indval_results <- as.data.frame(indval_htype2$sign) %>%
  tibble::rownames_to_column("taxon") %>%
  filter(p.value <= 0.05) %>%
  mutate(indicates = case_when(
    s.depositional == 1 & s.riffle == 0 ~ "depositional",
    s.riffle == 1 & s.depositional == 0 ~ "riffle",
    TRUE ~ "both"
  )) %>%
  filter(indicates != "both") %>%                    # exclusive indicators only
  filter(!taxon %in% excluded_taxa)                  # same exclusions as 7B

cat("Depositional indicators:", sum(indval_results$indicates == "depositional"), "\n")
cat("Riffle indicators:",       sum(indval_results$indicates == "riffle"),       "\n")

top5_htype2 <- indval_results %>%
  group_by(indicates) %>%
  slice_max(order_by = stat, n = 5) %>%
  ungroup()

species_scores_htype2 <- as.data.frame(scores(nmds, display = "species")) %>%
  tibble::rownames_to_column("taxon") %>%
  inner_join(top5_htype2, by = "taxon") %>%
  mutate(
    taxon    = str_to_sentence(taxon),
    fontface = ifelse(str_detect(taxon, "idae$"), "plain", "italic")
  )

# -----------------------------------------------------------
# 7B: ALL 6 HABITAT TYPES (htype)
# -----------------------------------------------------------
indval_htype <- multipatt(
  comm_scaled,
  meta$htype,
  func    = "IndVal.g",
  control = how(nperm = 999)
)

indval_htype_results <- as.data.frame(indval_htype$sign) %>%
  tibble::rownames_to_column("taxon") %>%
  filter(p.value <= 0.05)

htype_cols <- c("s.da", "s.rts", "s.sng", "s.ub", "s.veg", "s.riffle")

exclusive_only <- indval_htype_results %>%
  pivot_longer(cols = all_of(htype_cols),
               names_to  = "htype_col",
               values_to = "associated") %>%
  filter(associated == 1) %>%
  mutate(htype = str_remove(htype_col, "s\\.")) %>%
  group_by(taxon) %>%
  filter(n() == 1) %>%                               # exclusive indicators only
  ungroup() %>%
  filter(!taxon %in% excluded_taxa) %>%              # same exclusions as 7A
  group_by(htype) %>%
  slice_max(order_by = stat, n = 2, with_ties = FALSE) %>%
  ungroup()

cat("Exclusive indicators per habitat:\n")
print(exclusive_only %>% count(htype))

species_scores_htype <- as.data.frame(scores(nmds, display = "species")) %>%
  tibble::rownames_to_column("taxon") %>%
  inner_join(exclusive_only %>% select(taxon, htype), by = "taxon") %>%
  mutate(
    taxon    = str_to_sentence(taxon),
    fontface = ifelse(str_detect(taxon, "idae$"), "plain", "italic")
  )

# -----------------------------------------------------------
# 7.1a  BUILD TABLE DATA — htype2 (depositional vs riffle)
# -----------------------------------------------------------

sign_df2 <- as.data.frame(indval_htype2$sign) %>%
  tibble::rownames_to_column("taxon")

A_B2 <- data.frame(
  taxon = sign_df2$taxon,
  A = mapply(function(tax, idx) indval_htype2$A[tax, idx],
             sign_df2$taxon, sign_df2$index),
  B = mapply(function(tax, idx) indval_htype2$B[tax, idx],
             sign_df2$taxon, sign_df2$index)
)

table_data2 <- as.data.frame(indval_htype2$sign) %>%
  tibble::rownames_to_column("taxon") %>%
  filter(p.value <= 0.05) %>%
  filter(!taxon %in% excluded_taxa) %>%
  left_join(A_B2, by = "taxon") %>%
  mutate(
    indicates = case_when(
      s.depositional == 1 & s.riffle == 0 ~ "depositional",
      s.riffle == 1 & s.depositional == 0 ~ "riffle",
      TRUE ~ "both"
    )
  ) %>%
  filter(indicates != "both") %>%
  mutate(
    htype_label = htype_labels[indicates],
    taxon_fmt   = str_to_sentence(taxon),
    p_fmt       = ifelse(p.value < 0.001, "< 0.001",
                         format(round(p.value, 3), nsmall = 3)),
    stars = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  
  left_join(taxon_lookup %>% mutate(genus = tolower(Genus)),
            by = c("taxon" = "genus")) %>%
  mutate(
    Order_fct = factor(Order, levels = order_levels),
    Order_rank = ifelse(is.na(Order_fct), Inf, as.numeric(Order_fct))
  ) %>%
  arrange(htype_label, Order_rank, desc(stat)) %>%
  select(htype_label, Order, Family, taxon_fmt, A, B, stat, p_fmt, stars) %>%
  select(htype_label, Order, Family, taxon_fmt, A, B, stat, p_fmt, stars)

names(table_data2) <- c("Habitat", "Order", "Family", "Taxon", "A", "B", "IndVal", "p-value", "Sig")

# -----------------------------------------------------------
# 7.2a  BUILD FLEXTABLE
# -----------------------------------------------------------

ft2 <- flextable(table_data2) %>%
  
  merge_v(j = "Habitat") %>%
  valign(j = "Habitat", valign = "top") %>%
  set_header_labels(Sig = "") %>%
  
  colformat_double(j = c("A", "B", "IndVal"), digits = 2) %>%
  
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  bold(j = "Habitat", part = "body") %>%
  
  italic(
    i    = which(!str_detect(table_data2$Taxon, "idae$")),
    j    = "Taxon",
    part = "body"
  ) %>%
  
  align(part = "header", align = "center") %>%
  align(part = "body",   align = "center") %>%
  align(j = c("Habitat", "Order", "Family", "Taxon"), part = "body", align = "left") %>%
  # Add after width(j = "Taxon"):
  width(j = "Order",  width = 1.3) %>%
  width(j = "Family", width = 1.3) %>%
  
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header",
               border = fp_border(color = "black", width = 0.75)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.5)) %>%
  
  {
    hab_groups2 <- unique(table_data2$Habitat)
    ft_tmp <- .
    for (i in seq_along(hab_groups2)) {
      if (i %% 2 == 0) {
        rows <- which(table_data2$Habitat == hab_groups2[i])
        ft_tmp <- bg(ft_tmp, i = rows, bg = "#F5F5F5", part = "body")
      }
    }
    ft_tmp
  } %>%
  
  width(j = "Habitat", width = 2.2) %>%
  width(j = "Taxon",   width = 1.6) %>%
  width(j = c("A", "B"), width = 0.55) %>%
  width(j = "IndVal",  width = 0.75) %>%
  width(j = "p-value", width = 0.75) %>%
  width(j = "Sig",     width = 0.4)

# -----------------------------------------------------------
# 7.3a  EXPORT TO .docx AND .png
# -----------------------------------------------------------


#





doc2 <- read_docx() %>%
  body_add_par(
    paste0(
      "Table X. Exclusive indicator taxa for non-riffle (depositional) and riffle ",
      "habitat types identified by IndVal.g analysis (indicspecies; ",
      "De Cáceres & Legendre 2009). Only taxa exclusively associated with a ",
      "single habitat type are shown. ",
      "A = specificity (probability site belongs to group given taxon presence); ",
      "B = fidelity (probability of taxon presence given membership in group); ",
      "IndVal = \u221A(A \u00D7 B). Significance assessed by 999 permutations. ",
      "Chironomidae and Gastropoda excluded (ubiquitous; non-discriminating)."
    ),
    style = "Normal"
  ) %>%
  body_add_flextable(ft2)

print(doc2, target = "outputs/indval_htype2_table.docx")

print(ft2)

save_as_image(ft2,
              path   = "outputs/indval_htype2_table.png",
              zoom   = 3,
              expand = 10)

cat("Table saved to outputs/indval_htype2_table.docx\n")


# -----------------------------------------------------------
# 7.1b  BUILD TABLE DATA
# Pull stat, p.value back from the full results
# -----------------------------------------------------------

# Define labels outside the pipe
htype_labels <- c(
  da     = "Debris accumulation",
  rts    = "Submerged roots",
 # sng    = "Snag",
  ub     = "Undercut bank",
  veg    = "Submerged vegetation",
  riffle = "Riffle",
  depositional= "Non-riffle"
)

# Extract A and B using the index from $sign to select the right column
sign_df <- as.data.frame(indval_htype$sign) %>%
  tibble::rownames_to_column("taxon")

A_B <- data.frame(
  taxon = sign_df$taxon,
  A = mapply(function(tax, idx) indval_htype$A[tax, idx],
             sign_df$taxon, sign_df$index),
  B = mapply(function(tax, idx) indval_htype$B[tax, idx],
             sign_df$taxon, sign_df$index)
)

# Then join into table_data
table_data <- as.data.frame(indval_htype$sign) %>%
  tibble::rownames_to_column("taxon") %>%
  filter(p.value <= 0.05) %>%
  left_join(A_B, by = "taxon") %>%
  mutate(
    htype_label = case_when(
      s.da     == 1 ~ htype_labels["da"],
      s.rts    == 1 ~ htype_labels["rts"],
    #  s.sng    == 1 ~ htype_labels["sng"],
      s.ub     == 1 ~ htype_labels["ub"],
      s.veg    == 1 ~ htype_labels["veg"],
      s.riffle == 1 ~ htype_labels["riffle"],
      TRUE          ~ "Multiple"
    ),
    taxon_fmt = str_to_sentence(taxon),
    p_fmt = ifelse(p.value < 0.001, "< 0.001", format(round(p.value, 3), nsmall = 3)),
    stars = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  left_join(taxon_lookup %>% mutate(genus = tolower(Genus)),
            by = c("taxon" = "genus")) %>%
  mutate(
    Order_fct = factor(Order, levels = order_levels),
    Order_rank = ifelse(is.na(Order_fct), Inf, as.numeric(Order_fct))
  ) %>%
  arrange(htype_label, Order_rank, desc(stat)) %>%
  select(htype_label, Order, Family, taxon_fmt, A, B, stat, p_fmt, stars) %>%
  select(htype_label, Order, Family, taxon_fmt, A, B, stat, p_fmt, stars)

names(table_data) <- c("Habitat", "Order", "Family", "Taxon", "A", "B", "IndVal", "p-value", "Sig")

# -----------------------------------------------------------
# 7.2b BUILD FLEXTABLE
# -----------------------------------------------------------

ft <- flextable(table_data) %>%
  
  merge_v(j = "Habitat") %>%
  valign(j = "Habitat", valign = "top") %>%
  set_header_labels(Sig = "") %>%  
  
  # --- Numeric formatting ---
  colformat_double(j = c("A", "B", "IndVal"), digits = 2) %>%       # removed A and B
  
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  bold(j = "Habitat", part = "body") %>%
  
  italic(
    i    = which(!str_detect(table_data$Taxon, "idae$")),
    j    = "Taxon",
    part = "body"
  ) %>%
  
  align(part = "header", align = "center") %>%
  align(part = "body",   align = "center") %>%
  align(j = c("Habitat", "Order", "Family", "Taxon"), part = "body", align = "left") %>%
  
  width(j = "Order",  width = 1.3) %>%
  width(j = "Family", width = 1.3) %>%
  
  
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header",
               border = fp_border(color = "black", width = 0.75)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.5)) %>%
  
  {
    hab_groups <- unique(table_data$Habitat)
    ft_tmp <- .
    for (i in seq_along(hab_groups)) {
      if (i %% 2 == 0) {
        rows <- which(table_data$Habitat == hab_groups[i])
        ft_tmp <- bg(ft_tmp, i = rows, bg = "#F5F5F5", part = "body")
      }
    }
    ft_tmp
  } %>%
  
  # --- Column widths ---
  width(j = "Habitat", width = 2.2) %>%
  width(j = "Taxon",   width = 1.6) %>%
  width(j = "IndVal",  width = 0.75) %>%            # removed A and B
  width(j = "p-value", width = 0.75) %>%
  width(j = c("A", "B"), width = 0.55) %>%
  width(j = "IndVal",    width = 0.75) %>%
  width(j = "Sig",        width = 0.4)

# -----------------------------------------------------------
# 7.3b EXPORT TO .docx AND .png
# -----------------------------------------------------------

###################################################






#########################################################



doc <- read_docx() %>%
  body_add_par(
    paste0(
      "Table X. Exclusive indicator taxa per habitat type identified by ",
      "IndVal.g analysis (indicspecies; De Cáceres & Legendre 2009). ",
      "Only taxa exclusively associated with a single habitat type are shown. ",
      "A = specificity (probability site belongs to group given taxon presence); ",
      "B = fidelity (probability of taxon presence given membership in group); ",
      "IndVal = \u221A(A \u00D7 B). Significance assessed by 999 permutations. ",
      "Chironomidae and Gastropoda excluded (ubiquitous; non-discriminating)."
    ),
    style = "Normal"
  ) %>%
  body_add_flextable(ft)

print(doc, target = "outputs/indval_htype_table.docx")

print(ft)

save_as_image(ft,
              path   = "outputs/indval_htype_table.png",
              zoom   = 3,
              expand = 10)

cat("Table saved to outputs/indval_htype_table.docx\n")



# ============================================================
# SECTION 8: FIGURES
# Shared theme applied to all panels via theme_bw(base_size=20)
# ============================================================

shared_theme <- theme_bw(base_size = 12) +
  theme(
    # legend.position       = "bottom",
    # Option A: black major only, drop minor entirely
    panel.grid.major = element_line(colour = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    # legend.text           = element_text(size = 18),
    # legend.title          = element_text(size = 18),
    # axis.text             = element_text(size = 20),
    # axis.title            = element_text(size = 20),
    # plot.caption          = element_text(size = 18, hjust = 0),
    plot.caption.position = "plot"
  )



# -----------------------------------------------------------
# FIGURE 1: Depositional vs riffle — centroids + 95% ellipses
# -----------------------------------------------------------

fig2_dep_rif <- ggplot(nmds_scores,
                       aes(x = NMDS1, y = NMDS2, color = htype2 )) +
  stat_ellipse(aes(group = htype2), type = "t", level = 0.95, linewidth = 1,
               linetype=2) +
  geom_point(alpha = 0.7, size = 5, ) +
  geom_point(data = centroids_htype2,
             aes(x = NMDS1, y = NMDS2, fill = htype2),
             color="black", stroke=1.5,
             size = 8, shape = 23, show.legend = FALSE ) +
  scale_color_manual(values = colors_htype2, labels = htype_labels) +
  scale_fill_manual(values  = colors_htype2, labels = htype_labels) +
  labs(
    x       = paste0("NMDS1 (R² = ", r2_axis1, ")"),
    y       = paste0("NMDS2 (R² = ", r2_axis2, ")"),
    color   = "Habitat Type",
    fill    = "Habitat Type",
    caption = paste0(stress_label, "     |     Diamonds = centroids | Ellipses = 95% CI"),
  ) +
  shared_theme+
  theme(
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=10),
    legend.text = element_text(size=10),
    panel.background = element_rect(fill= "transparent"),
    plot.background = element_rect(fill="transparent"),
    legend.background = element_rect(fill="transparent"),
    legend.box.background = element_rect(fill="transparent"),
  )

fig2_dep_rif
ggsave("outputs/H3_figs/H3_fig2_dep_rif_centroids.png", fig2_dep_rif,
       width = 16, height = 8, dpi = 300, bg = "transparent")

# -----------------------------------------------------------
# FIGURE 2: All 6 habitat types — centroids + 95% ellipses
# -----------------------------------------------------------

fig5_6htype <- ggplot(nmds_scores,
                      aes(x = NMDS1, y = NMDS2, color = htype, fill= "transparent")) +
  stat_ellipse(aes(group = htype, color = htype), type = "t", level = 0.95,
               linewidth = 0.8, linetype=2)+
  geom_point(alpha =  0.7, size = 5) +
  geom_point(data = centroids_htype,
             aes(x = NMDS1, y = NMDS2, fill = htype),
             color="black", stroke=1.5,
             size = 8, shape = 23, show.legend = FALSE ) +
  scale_color_manual(values = colors_htype, labels = htype_labels) +
  scale_fill_manual(values  = colors_htype, labels = htype_labels) +
  labs(
    x       = paste0("NMDS1 (R² = ", r2_axis1, ")"),
    y       = paste0("NMDS2 (R² = ", r2_axis2, ")"),
    color   = "Habitat Type",
    fill    = "Habitat Type",
    caption = paste0(stress_label, " | Diamonds = centroids | Ellipses = 95% CI")
  ) +
  shared_theme

fig5_6htype
ggsave("outputs/H3_figs/H3_fig5_6htype_centroids.png", fig5_6htype,
       width = 16, height = 8, dpi = 300, bg = "transparent")

#===================================================================#
#cowplot for the two plots I want for the paper 8.26.2026
#===================================================================#

# one shared caption for the whole figure
shared_caption <- paste0(stress_label,
                         "     |     Diamonds = centroids | Ellipses = 95% CI")

# make both panels visually consistent + strip their individual captions
transparent_bits <- theme(
  axis.text.x           = element_text(size = 12),
  axis.text.y           = element_text(size = 12),
  legend.text           = element_text(size = 12),
  legend.title          = element_text(size = 11),
  plot.title            = element_text(size = 12, face = "bold", hjust = 0.5),
  panel.background      = element_rect(fill = "transparent", color = NA),
  plot.background       = element_rect(fill = "transparent", color = NA),
  legend.background     = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA)
)

p_left <- fig2_dep_rif +
  labs(caption = NULL, title = "Riffle and other habitats") +
  transparent_bits

p_right <- fig5_6htype +
  labs(caption = NULL, title = "All habitats") +
  transparent_bits

# side-by-side panels
panels <- plot_grid(
  p_left, p_right,
  labels     = c("A", "B"),
  label_size = 13,
  ncol       = 1,
  align      = "v",
  axis       = "lr"
)

# caption as its own thin strip underneath
caption_row <- ggdraw() +
  draw_label(shared_caption, size = 12, x = 0.5, hjust = 0.5)

fig6_combined <- plot_grid(
  panels, caption_row,
  ncol        = 1,
  rel_heights = c(1, 0.05)
) +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

fig6_combined

ggsave("outputs/H3_figs/H3_fig6_combined_dep_rif_6htype.png", fig6_combined,
       width = 10, height =10, dpi = 300, bg = "transparent")



# ============================================================
# SECTION 9: PERMANOVA
# Tests whether assemblage composition differs significantly
# by habitat type and/or stream identity.
# Models are run sequentially (by = "terms") so variance is
# partitioned in order: Site first, then htype.
# Betadisper checks the homogeneity-of-dispersion assumption;
# a significant result means groups differ in spread, not just
# location — interpret PERMANOVA results with caution if so.
# ============================================================

# Helper: convert adonis2 output to a tidy summary data frame
tidy_adonis <- function(model_label, adonis_obj) {
  as.data.frame(adonis_obj) %>%
    mutate(variable = rownames(.)) %>%
    filter(variable != "Total") %>%
    transmute(
      Model    = model_label,
      Variable = variable,
      SumOfSqs = round(SumOfSqs, 2),
      R2       = round(R2, 3),
      F        = round(`F`, 3),
      p.value  = `Pr(>F)`
    )
}

set.seed(123)

# -----------------------------------------------------------
# PART A: htype2 — depositional vs riffle
# -----------------------------------------------------------

cat("\n--- PART A: htype2 (depositional vs riffle) ---\n")

perm_a1 <- adonis2(bc_dist ~ htype2,            data = meta, permutations = 999, by = "terms")
perm_a2 <- adonis2(bc_dist ~ Site,              data = meta, permutations = 999, by = "terms")
perm_a3 <- adonis2(bc_dist ~ Site + htype2,     data = meta, permutations = 999, by = "terms")
perm_a4 <- adonis2(bc_dist ~ Site + htype2 + Site:htype2,
                                                 data = meta, permutations = 999, by = "terms")




# Betadisper: test whether depositional and riffle groups have equal dispersion
bd_htype2      <- betadisper(bc_dist, meta$htype2)
bd_htype2_perm <- permutest(bd_htype2, permutations = 999)
cat("\nBetadisper — htype2:\n"); print(bd_htype2_perm)

bd_htype2_perm

table_a <- bind_rows(
  tidy_adonis("BC ~ htype2",                      perm_a1),
  tidy_adonis("BC ~ Site",                        perm_a2),
  tidy_adonis("BC ~ Site + htype2",               perm_a3),
  tidy_adonis("BC ~ Site + htype2 + Site:htype2", perm_a4)
)
print(table_a, row.names = FALSE)
write.csv(table_a, "outputs/H3_PERMANOVA_htype2_summary.csv", row.names = FALSE)

# -----------------------------------------------------------
# PART B: htype — all 6 habitat types
# -----------------------------------------------------------

cat("\n--- PART B: htype (6 habitat types) ---\n")

perm_b1 <- adonis2(bc_dist ~ htype,            data = meta, permutations = 999, by = "terms")
perm_b2 <- adonis2(bc_dist ~ Site,             data = meta, permutations = 999, by = "terms")
perm_b3 <- adonis2(bc_dist ~ Site + htype,     data = meta, permutations = 999, by = "terms")
perm_b4 <- adonis2(bc_dist ~ Site + htype + Site:htype,
                                                data = meta, permutations = 999, by = "terms")

# Betadisper: test whether the 6 habitat groups have equal dispersion
bd_htype      <- betadisper(bc_dist, meta$htype)
bd_htype_perm <- permutest(bd_htype, permutations = 999)
cat("\nBetadisper — htype:\n"); print(bd_htype_perm)

table_b <- bind_rows(
  tidy_adonis("BC ~ htype",                    perm_b1),
  tidy_adonis("BC ~ Site",                     perm_b2),
  tidy_adonis("BC ~ Site + htype",             perm_b3),
  tidy_adonis("BC ~ Site + htype + Site:htype", perm_b4)
)
print(table_b, row.names = FALSE)
write.csv(table_b, "outputs/H3_PERMANOVA_htype_summary.csv", row.names = FALSE)

# -----------------------------------------------------------
# Combined table for reporting
# -----------------------------------------------------------

bind_rows(
  table_a %>% mutate(Grouping = "htype2 (dep vs rif)"),
  table_b %>% mutate(Grouping = "htype (6 types)")
) %>%
  select(Grouping, Model, Variable, SumOfSqs, R2, p.value) %>%
  write.csv("outputs/H3_PERMANOVA_combined.csv", row.names = FALSE)

cat("\nAll outputs saved to outputs/\n")

# ============================================================
# SECTION 10: PERMANOVA + DISPERSION SUMMARY TABLES
# Combines the adonis2 model-building sequence with the
# betadisper/permutest dispersion test into two publication-
# ready flextables:
#   Table A: htype2 (non-riffle vs riffle)
#   Table B: htype  (all 6 mesohabitat types)
#
# Run this AFTER Section 9 (PERMANOVA) in H3_NMDS_Biomass.R —
# it reuses table_a, table_b, bd_htype2_perm, bd_htype_perm,
# meta, and htype_labels from that section.
# ============================================================

# -----------------------------------------------------------
# 10.0 HELPERS
# -----------------------------------------------------------

# Pulls the "Groups" row out of a permutest.betadisper object
# and formats it to match the tidy_adonis() column structure,
# so it can be row-bound onto the PERMANOVA model table.
tidy_betadisper <- function(model_label, bd_perm_obj) {
  tab <- as.data.frame(bd_perm_obj$tab)
  tab %>%
    tibble::rownames_to_column("variable") %>%
    filter(variable == "Groups") %>%
    transmute(
      Model    = model_label,
      Variable = "Dispersion (PERMDISP)",
      SumOfSqs = round(`Sum Sq`, 2),
      R2       = NA_real_,
      F        = round(`F`, 3),
      p.value  = `Pr(>F)`
    )
}

# Adds a formatted p-value column and significance stars,
# consistent with the indicator-species tables (ft, ft2)
add_stars <- function(df) {
  df %>%
    mutate(
      p_fmt = ifelse(p.value < 0.001, "< 0.001", format(round(p.value, 3), nsmall = 3)),
      stars = case_when(
        p.value <= 0.001 ~ "***",
        p.value <= 0.01  ~ "**",
        p.value <= 0.05  ~ "*",
        TRUE              ~ ""
      )
    )
}

# Shared flextable builder — takes the combined data frame and
# a vector of row indices where a new Model block starts (for
# shading), and returns a formatted flextable.
build_permanova_ft <- function(df) {
  
  ft <- flextable(df) %>%
    merge_v(j = "Model") %>%
    valign(j = "Model", valign = "top") %>%
    set_header_labels(
      Model    = "Model",
      Variable = "Term",
      SumOfSqs = "Sum Sq",
      R2       = "R\u00b2",
      F        = "F",
      p_fmt    = "p",
      stars    = ""
    ) %>%
    colformat_double(j = c("SumOfSqs", "F"), digits = 2) %>%
    colformat_double(j = "R2", digits = 3, na_str = "\u2014") %>%
    
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    bold(j = "Model", part = "body") %>%
    italic(i = which(df$Variable == "Dispersion (PERMDISP)"),
           j = "Variable", part = "body") %>%
    
    align(part = "header", align = "center") %>%
    align(part = "body",   align = "center") %>%
    align(j = c("Model", "Variable"), part = "body", align = "left") %>%
    
    border_remove() %>%
    hline_top(part = "header",
              border = fp_border(color = "black", width = 1.5)) %>%
    hline_bottom(part = "header",
                 border = fp_border(color = "black", width = 0.75)) %>%
    hline_bottom(part = "body",
                 border = fp_border(color = "black", width = 1.5)) %>%
    
    # thin rule under each model block, and above the dispersion row
    hline(i = which(df$Variable == "Dispersion (PERMDISP)") - 1,
          border = fp_border(color = "grey60", width = 0.5), part = "body") %>%
    
    width(j = "Model",    width = 1.9) %>%
    width(j = "Variable", width = 1.7) %>%
    width(j = "SumOfSqs", width = 0.8) %>%
    width(j = "R2",       width = 0.7) %>%
    width(j = "F",        width = 0.7) %>%
    width(j = "p_fmt",    width = 0.7) %>%
    width(j = "stars",    width = 0.35)
  
  # alternating shading by Model block
  model_groups <- unique(df$Model)
  for (i in seq_along(model_groups)) {
    if (i %% 2 == 0) {
      rows <- which(df$Model == model_groups[i])
      ft <- bg(ft, i = rows, bg = "#F5F5F5", part = "body")
    }
  }
  ft
}

# -----------------------------------------------------------
# 10A. TABLE A — htype2 (non-riffle vs riffle)
# -----------------------------------------------------------

table_a_full <- bind_rows(
  table_a,
  tidy_betadisper("Dispersion test", bd_htype2_perm)
) %>%
  add_stars() %>%
  select(Model, Variable, SumOfSqs, R2, F, p_fmt, stars)

ft_a <- build_permanova_ft(table_a_full)

doc_a <- read_docx() %>%
  body_add_par(
    paste0(
      "Table X. PERMANOVA model-building sequence and multivariate dispersion ",
      "test (PERMDISP) for macroinvertebrate biomass composition between non-riffle ",
      "and riffle habitats. Models were fit sequentially on Bray-Curtis dissimilarities ",
      "via adonis2 (vegan; by = \"terms\", 999 permutations), partitioning variance in the ",
      "order terms are entered. Homogeneity of multivariate dispersion was assessed via ",
      "betadisper/permutest (999 permutations); group sample sizes were balanced (n = 105 per group)."
    ),
    style = "Normal"
  ) %>%
  body_add_flextable(ft_a)

print(doc_a, target = "outputs/H3_PERMANOVA_dispersion_htype2_table.docx")
print(ft_a)
save_as_image(ft_a,
              path   = "outputs/H3_PERMANOVA_dispersion_htype2_table.png",
              zoom   = 3,
              expand = 10)

cat("Table A saved to outputs/H3_PERMANOVA_dispersion_htype2_table.docx\n")

# -----------------------------------------------------------
# 10B. TABLE B — htype (all 6 mesohabitat types)
# -----------------------------------------------------------

table_b_full <- bind_rows(
  table_b,
  tidy_betadisper("Dispersion test", bd_htype_perm)
) %>%
  add_stars() %>%
  select(Model, Variable, SumOfSqs, R2, F, p_fmt, stars)

ft_b <- build_permanova_ft(table_b_full)

doc_b <- read_docx() %>%
  body_add_par(
    paste0(
      "Table X. PERMANOVA model-building sequence and multivariate dispersion ",
      "test (PERMDISP) for macroinvertebrate biomass composition among six mesohabitat ",
      "types (debris accumulation, submerged roots, snag, undercut bank, vegetation, riffle). ",
      "Models were fit sequentially on Bray-Curtis dissimilarities via adonis2 (vegan; ",
      "by = \"terms\", 999 permutations), partitioning variance in the order terms are entered. ",
      "Homogeneity of multivariate dispersion was assessed via betadisper/permutest ",
      "(999 permutations)."
    ),
    style = "Normal"
  ) %>%
  body_add_flextable(ft_b)

print(doc_b, target = "outputs/H3_PERMANOVA_dispersion_htype_table.docx")
print(ft_b)
save_as_image(ft_b,
              path   = "outputs/H3_PERMANOVA_dispersion_htype_table.png",
              zoom   = 3,
              expand = 10)

cat("Table B saved to outputs/H3_PERMANOVA_dispersion_htype_table.docx\n")


# ============================================================
# SECTION 11: PAIRWISE PERMANOVA + SNG (n=5) SENSITIVITY CHECK
#
# Two things happen here:
#   11A. Pairwise PERMANOVA across all 6 htype levels, Holm-
#        corrected, with pairs involving `sng` flagged as
#        low-power (n = 5) rather than silently trusted.
#   11B. Sensitivity re-run of the omnibus PERMANOVA + betadisper
#        with `sng` excluded entirely, to check whether the
#        htype effect and its magnitude hold without it.
#
# Run after Section 9. Reuses bc_dist, meta, comm_scaled,
# htype_labels from earlier sections.
# ============================================================

library(purrr)

# -----------------------------------------------------------
# 11A. PAIRWISE PERMANOVA (all 6 htype levels)
# -----------------------------------------------------------

pairwise_permanova <- function(dist_obj, group_var, perms = 999,
                               low_power_n = 10) {
  
  group_var <- as.character(group_var)
  dist_mat  <- as.matrix(dist_obj)
  levels_g  <- sort(unique(group_var))
  pairs     <- combn(levels_g, 2, simplify = FALSE)
  
  n_by_group <- table(group_var)
  
  results <- map_dfr(pairs, function(p) {
    idx       <- group_var %in% p
    sub_dist  <- as.dist(dist_mat[idx, idx])
    sub_group <- factor(group_var[idx], levels = p)
    
    mod <- adonis2(sub_dist ~ sub_group, permutations = perms)
    
    tibble(
      Group1 = p[1],
      Group2 = p[2],
      n1     = n_by_group[[p[1]]],
      n2     = n_by_group[[p[2]]],
      SumOfSqs = round(mod$SumOfSqs[1], 2),
      R2       = round(mod$R2[1], 3),
      F        = round(mod$`F`[1], 3),
      p.raw    = mod$`Pr(>F)`[1]
    )
  })
  
  results %>%
    mutate(
      p.adj      = round(p.adjust(p.raw, method = "holm"), 3),
      low_power  = (n1 < low_power_n | n2 < low_power_n),
      stars      = case_when(
        p.adj <= 0.001 ~ "***",
        p.adj <= 0.01  ~ "**",
        p.adj <= 0.05  ~ "*",
        TRUE           ~ ""
      ),
      flag = ifelse(low_power, "low n \u2014 interpret cautiously", "")
    ) %>%
    arrange(desc(low_power), p.adj)
}

set.seed(123)
pairwise_htype <- pairwise_permanova(bc_dist, meta$htype, perms = 999, low_power_n = 10)

cat("\n--- Pairwise PERMANOVA: htype (Holm-corrected) ---\n")
print(pairwise_htype, n = Inf)

write.csv(pairwise_htype, "outputs/H3_PERMANOVA_pairwise_htype.csv", row.names = FALSE)

# ---- Flextable version for the manuscript / supplement ----

pairwise_htype_ft_data <- pairwise_htype %>%
  mutate(
    Habitat1 = htype_labels[Group1],
    Habitat2 = htype_labels[Group2],
    n1_n2    = paste0(n1, " vs ", n2),
    p_fmt    = ifelse(p.adj < 0.001, "< 0.001", format(round(p.adj, 3), nsmall = 3))
  ) %>%
  select(Habitat1, Habitat2, n1_n2, SumOfSqs, R2, F, p_fmt, stars, flag)

names(pairwise_htype_ft_data) <- c("Habitat 1", "Habitat 2", "n (1 vs 2)",
                                   "Sum Sq", "R\u00b2", "F", "p (Holm)", "", "Flag")

ft_pairwise <- flextable(pairwise_htype_ft_data) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 9, part = "body") %>%
  fontsize(size = 9, part = "header") %>%
  bold(part = "header") %>%
  italic(i = which(pairwise_htype_ft_data$Flag != ""), j = "Flag", part = "body") %>%
  color(i = which(pairwise_htype_ft_data$Flag != ""), color = "#8B3A3A", part = "body") %>%
  align(part = "header", align = "center") %>%
  align(part = "body",   align = "center") %>%
  align(j = c("Habitat 1", "Habitat 2", "Flag"), part = "body", align = "left") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 0.75)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
  bg(i = which(pairwise_htype_ft_data$Flag != ""), bg = "#FBEAEA", part = "body") %>%
  width(j = c("Habitat 1", "Habitat 2"), width = 1.4) %>%
  width(j = "Flag", width = 1.6)

doc_pairwise <- read_docx() %>%
  body_add_par(
    paste0(
      "Table X. Pairwise PERMANOVA comparisons among all six mesohabitat types, based on ",
      "Bray-Curtis dissimilarities of biomass composition (adonis2, 999 permutations per pair; ",
      "p-values Holm-corrected for 15 comparisons). Comparisons involving the snag (sng) habitat ",
      "type (n = 5) are flagged and highlighted; given the small sample size for this group, ",
      "these pairwise results should be interpreted with caution and are not treated as strong ",
      "evidence of compositional difference or equivalence."
    ),
    style = "Normal"
  ) %>%
  body_add_flextable(ft_pairwise)

print(doc_pairwise, target = "outputs/H3_PERMANOVA_pairwise_htype_table.docx")
save_as_image(ft_pairwise,
              path   = "outputs/H3_PERMANOVA_pairwise_htype_table.png",
              zoom   = 3,
              expand = 10)

cat("Pairwise table saved to outputs/H3_PERMANOVA_pairwise_htype_table.docx\n")

# -----------------------------------------------------------
# 11B. SENSITIVITY CHECK — omnibus PERMANOVA + betadisper
#      with sng (n=5) EXCLUDED
# -----------------------------------------------------------

keep_idx      <- meta$htype != "sng"
meta_nosng    <- meta[keep_idx, ] %>% mutate(htype = droplevels(factor(htype)))
comm_nosng    <- comm_scaled[keep_idx, ]
bc_dist_nosng <- vegdist(comm_nosng, method = "bray")

cat("\n--- Sensitivity check: htype with sng excluded ---\n")
meta_nosng %>% count(htype) %>% print()

set.seed(123)
perm_nosng_1 <- adonis2(bc_dist_nosng ~ htype,        data = meta_nosng, permutations = 999, by = "terms")
perm_nosng_2 <- adonis2(bc_dist_nosng ~ Site,          data = meta_nosng, permutations = 999, by = "terms")
perm_nosng_3 <- adonis2(bc_dist_nosng ~ Site + htype,  data = meta_nosng, permutations = 999, by = "terms")

bd_nosng      <- betadisper(bc_dist_nosng, meta_nosng$htype)
bd_nosng_perm <- permutest(bd_nosng, permutations = 999)

cat("\nPERMANOVA (htype, sng excluded):\n"); print(perm_nosng_1)
cat("\nBetadisper (htype, sng excluded):\n"); print(bd_nosng_perm)

table_nosng <- bind_rows(
  tidy_adonis("BC ~ htype (no sng)",        perm_nosng_1),
  tidy_adonis("BC ~ Site (no sng)",         perm_nosng_2),
  tidy_adonis("BC ~ Site + htype (no sng)", perm_nosng_3),
  tidy_betadisper("Dispersion test (no sng)", bd_nosng_perm)
) %>%
  add_stars() %>%
  select(Model, Variable, SumOfSqs, R2, F, p_fmt, stars)

write.csv(table_nosng, "outputs/H3_PERMANOVA_sensitivity_no_sng.csv", row.names = FALSE)

cat("\n--- Comparison: full model vs. sng-excluded sensitivity check ---\n")
cat("Full model    (BC ~ htype): R2 =", round(perm_b1$R2[1], 3),
    " F =", round(perm_b1$`F`[1], 3), " p =", perm_b1$`Pr(>F)`[1], "\n")
cat("No-sng model  (BC ~ htype): R2 =", round(perm_nosng_1$R2[1], 3),
    " F =", round(perm_nosng_1$`F`[1], 3), " p =", perm_nosng_1$`Pr(>F)`[1], "\n")

cat("\nFull model    dispersion F =", round(bd_htype_perm$tab$F[1], 2),
    " p =", bd_htype_perm$tab$`Pr(>F)`[1], "\n")
cat("No-sng model  dispersion F =", round(bd_nosng_perm$tab$F[1], 2),
    " p =", bd_nosng_perm$tab$`Pr(>F)`[1], "\n")




