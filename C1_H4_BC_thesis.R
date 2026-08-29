# ============================================================
# SCRIPT: C1_H4_BC_thesis_traitbiomass_v2.R
#
# TRAIT-GROUP BIOMASS  (direct summation)
#
#   trait_biomass_i = sum(trait_j × Biomass_m2_j)   for all taxa j
#
#   where trait_j is a binary (0/1) affinity for trait i.
#   Result: absolute biomass of taxa possessing each trait (mg DM m^-2).
#   Taxa coded 1 for multiple traits contribute their biomass to each
#   relevant trait category independently (reflects biological reality
#   of overlapping functional roles under a binary coding scheme;
#   see Poff et al. 2006; Cummins et al. 2022).
#
#   PART A: Riffle vs. Depositional, site-level averaged. Stats: Wilcoxon.
#   PART B: All 6 habitat types, replicate-level. Stats: KW + Dunn (BH).
#   PART C: Combined multi-panel figures (cowplot).
#   PART D: KW + Dunn results tables (.csv + flextable .docx).
#
# INPUTS:
#   - outputs/rep_biomass_by_taxon_4.7.26.csv
#   - outputs/traits_final.csv
# OUTPUTS:
#   - outputs/H4.1_trait_biomass_collapsed.csv
#   - outputs/H4.1_trait_biomass_by_habtype.csv
#   - outputs/H4_4.15/H4.1_traitbio_FFG_combined.png
#   - outputs/H4_4.15/H4.1_traitbio_Habit_combined.png
#   - outputs/H4_4.15/H4.1_traitbio_KW_Dunn_results.csv
#   - outputs/H4_4.15/H4.1_traitbio_KW_Dunn_table.docx
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstatix)
library(purrr)
library(scales)

if (!requireNamespace("multcompView", quietly = TRUE)) install.packages("multcompView")
if (!requireNamespace("cowplot",      quietly = TRUE)) install.packages("cowplot")
library(multcompView)
library(cowplot)

dir.create("outputs/H4_4.15", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SECTION 1: LOAD DATA
# ============================================================

all_biomass  <- read.csv("outputs/rep_biomass_by_taxon_4.7.26.csv",
                         stringsAsFactors = FALSE)
traits_final <- read.csv("outputs/traits_final.csv",
                         stringsAsFactors = FALSE)

# ============================================================
# SECTION 2: TAXON ROLLUP MAPS  (Chessie BIBI, Smith et al. 2017)
# ============================================================
rollup_map <- c(
  "oligochaeta"   = "clitellata",
  "hirundinea"    = "clitellata",
  "annelida"      = "clitellata",
  "corbicula"     = "bivalvia",
  "pisidiidae"    = "bivalvia",
  "sphaeriidae"   = "bivalvia",
  "psididiidae"   = "bivalvia",
  "ancylidae"     = "gastropoda",
  "hydrobiidae"   = "gastropoda",
  "lymnaeidae"    = "gastropoda",
  "physa"         = "gastropoda",
  "physidae"      = "gastropoda",
  "planorbidae"   = "gastropoda",
  "pleuroceridae" = "gastropoda",
  "viviparidae"   = "gastropoda",
  "bithyniidae"   = "gastropoda",
  "turbellaria"   = "turbellaria"
)
family_rollup <- c(
  "aluetes"        = "dryopidae",
  "argiogomphus"   = "gomphidae",
  "ceratopogoniae" = "ceratopogonidae",
  "hydraticus"     = "dytiscidae",
  "prasocuris"     = "chrysomelidae",
  "sphenophorus"   = "curculionidae",
  "tricyphonoa"    = "limoniidae",
  "galenorish"     = "corixidae",
  "isaweon"        = "baetidae",
  "sacodes"        = "scirtidae",
  "sarabandus"     = "scirtidae"
)

# Genus proxy assignments (mirrors Trait_Generation.R Section 1G).
# NOTE: these are NOT applied in the mutate() rollup below, and that's
# intentional. Trait_Generation.R's lookup_taxon_traits() writes the
# proxy-derived trait row back out under the ORIGINAL taxon name
# (TAXANAME1 = "pedicidae"), borrowing only Dicranota's/Psychoda's FFG
# and habit values — not the genus label itself. So traits_final already
# has a row keyed "pedicidae" and a row keyed "psychodini", and the
# left_join below will match on those names with no renaming needed.
# Listed here for documentation/QC only.
genus_proxy_map <- c(
  "pedicidae"  = "dicranota",
  "psychodini" = "psychoda"
)

exclude_always <- c(
  "hydracarina", "copepoda", "daphnia",
  "ostracod", "saldidae", "collembola"
)

# pedicidae and psychodini removed — they now have real trait data via
# genus_proxy_map (see Trait_Generation.R Section 1G) and should flow
# through to trait-based metrics like any other matched taxon.
document_only <- c(
  "diptera", "ephemeroptera", "plecoptera"
)

all_biomass <- all_biomass %>%
  filter(!TAXANAME1 %in% exclude_always) %>%
  mutate(TAXANAME1 = case_when(
    TAXANAME1 %in% names(rollup_map)    ~ rollup_map[TAXANAME1],
    TAXANAME1 %in% names(family_rollup) ~ family_rollup[TAXANAME1],
    TRUE ~ TAXANAME1
  ))

cat("--- Rollup QC ---\n")
cat("Unique taxa after rollup:", n_distinct(all_biomass$TAXANAME1), "\n")
cat("Genus proxy taxa retained under original name (traits borrowed,",
    "not rolled up):", paste(names(genus_proxy_map), collapse = ", "), "\n\n")

# Sanity check: confirm genus proxy taxa present in the biomass data
# actually made it into traits_final under their own name (catches the
# case where Trait_Generation.R wasn't re-run after this change).
# NOTE: requires traits_final (from outputs/traits_final.csv) to already
# be loaded in this session — move this block after that read_csv() call
# if it isn't loaded yet at this point in your script.
missing_proxy_traits <- names(genus_proxy_map)[
  names(genus_proxy_map) %in% all_biomass$TAXANAME1 &
    !names(genus_proxy_map) %in% traits_final$TAXANAME1
]
if (length(missing_proxy_traits) > 0) {
  warning(
    "Genus proxy taxa present in biomass but missing from traits_final — ",
    "re-run Trait_Generation.R: ",
    paste(missing_proxy_traits, collapse = ", ")
  )
}
# ============================================================
# SECTION 3: COLOR PALETTES
# ============================================================



hab_colors <- c(
    "ub"     = "#c2b311",
    "da"     = "#4575B4",
    "rts"    = "#F97D0B",
    #"sng"    = "#D7191C",
    "veg"    = "#7B2D8B",
    "riffle" = "#4DAC26"
  )

collapsed_colors <- c("depositional" = "#8B5E3C", "riffle" = "#84C441")
collapsed_names  <- c("depositional" = "Depositional", "riffle" = "Riffle")

hab_full_names <- c(
  "ub"     = "UB",
  "da"     = "DA",
  "rts"    = "RTS",
 # "sng"    = "Snag",
  "veg"    = "VEG",
  "riffle" = "Riffle"
)

htyp2_names <- c(
  "depositional" = "Non-riffle",
  "riffle" = "Riffle"
)

# ============================================================
# SECTION 4: SPLIT + RENAMES + ALL_RAW
# ============================================================


dep_biomass <- all_biomass %>% filter(htype != "riffle")
rif_biomass <- all_biomass %>% filter(htype == "riffle")

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

dep_biomass <- dep_biomass %>% 
  mutate(htype = if_else(htype == "sng","da", htype))

unique(dep_biomass$htype)

dep_raw_labeled <- dep_biomass %>% mutate(htype = tolower(htype))
rif_labeled     <- rif_biomass %>% mutate(htype = "riffle")
all_raw         <- bind_rows(dep_raw_labeled, rif_labeled)

# ============================================================
# SECTION 5: TRAIT MATRIX + TRAIT VECTORS
# ============================================================

trait_matrix <- traits_final %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -9999))) %>%
  dplyr::select(TAXANAME1,
                starts_with("Feed_"),
                starts_with("Habit_")) %>%
  dplyr::select(-any_of(c("Feed_Breadth", "Habit_Breadth"))) %>%
  distinct(TAXANAME1, .keep_all = TRUE) %>%
  mutate(across(-TAXANAME1, ~ replace_na(.x, 0)))

trait_cols <- names(trait_matrix)[-1]

feeding_groups <- c(
  "Collector-filterer", "Collector-gatherer",
  "Scraper/Grazer", "Predator", "Shredder", "Herbivore/Peircer"
)

habits <- c(
  "Burrower", "Climber", "Clinger", "Crawler",
  "Skater", "Sprawler", "Swimmer"
)

label_trait <- function(trait) dplyr::case_when(
  trait == "Feed_CF"          ~ "Collector-filterer",
  trait == "Feed_CG"          ~ "Collector-gatherer",
  trait == "Feed_SC"          ~ "Scraper/Grazer",
  trait == "Feed_PH"          ~ "Piercer/Herbivore",
  trait == "Feed_PA"          ~ "Parasite",
  trait == "Feed_PR"          ~ "Predator",
  trait == "Feed_SH"          ~ "Shredder",
  trait == "Habit_Burrower"   ~ "Burrower",
  trait == "Habit_Climber"    ~ "Climber",
  trait == "Habit_Clinger"    ~ "Clinger",
  trait == "Habit_Planktonic" ~ "Planktonic",
  trait == "Habit_Skater"     ~ "",
  trait == "Habit_Sprawler"   ~ "Sprawler",
  trait == "Habit_Swimmer"    ~ "Swimmer",
  TRUE ~ trait
)

# ============================================================
# SECTION 6: SITE-LEVEL BIOMASS (collapsed: riffle vs depositional)
# Biomass summed per taxon across reps within a site, divided by
# n reps → mean per-m² value. Unit of replication: Site.
# ============================================================

raw_collapsed <- all_raw %>%
  mutate(
    htype_collapsed = ifelse(htype == "riffle", "riffle", "depositional"),
    htype_collapsed = factor(htype_collapsed,
                             levels = c("depositional", "riffle"))
  )

n_reps_site <- raw_collapsed %>%
  group_by(Site, htype_collapsed) %>%
  summarise(n_reps = n_distinct(Sample), .groups = "drop")

raw_collapsed_site <- raw_collapsed %>%
  group_by(Site, htype_collapsed, TAXANAME1) %>%
  summarise(Biomass_m2 = sum(Biomass_m2, na.rm = TRUE), .groups = "drop") %>%
  left_join(n_reps_site, by = c("Site", "htype_collapsed")) %>%
  mutate(Biomass_m2 = Biomass_m2 / n_reps) %>%
  dplyr::select(-n_reps)

# ============================================================
# SECTION 7: DIRECT TRAIT-GROUP BIOMASS CALCULATION
#
# For each replicate (or site), trait_biomass for trait i is:
#   sum(binary_trait_i × Biomass_m2) across all taxa
#
# Taxa coded 1 for multiple traits contribute their biomass to
# each relevant category independently, consistent with the
# binary coding approach of Poff et al. (2006) and the FFG
# biomass summation framework of Cummins et al. (2022) and
# Yegon et al. (2021).
# ============================================================

trait_bio_taxa <- raw_collapsed_site %>%
  filter(!TAXANAME1 %in% document_only) %>%
  left_join(trait_matrix, by = "TAXANAME1") %>%
  group_by(TAXANAME1)
  summarise(
    across(
      all_of(trait_cols),
      ~ sum(.x * Biomass_m2, na.rm = TRUE),
      .names = "traitbio_{.col}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(cols      = starts_with("traitbio_"),
               names_to  = "Trait",
               values_to = "trait_biomass") %>%
  mutate(
    Trait      = gsub("^traitbio_", "", Trait),
    TraitGroup = case_when(
      startsWith(Trait, "Feed_")  ~ "Feeding Group",
      startsWith(Trait, "Habit_") ~ "Habit"
    ),
    TraitLabel = label_trait(Trait)
  ) %>%
  filter(!is.na(TraitGroup)) 


  write.csv(trait_bio_taxa, "outputs/taxa_trait_bio_7.27.2")
  
# ---- Part A: site-level (collapsed: riffle vs. depositional) ----
traitbio_collapsed <- raw_collapsed_site %>%
  filter(!TAXANAME1 %in% document_only) %>%
  left_join(trait_matrix, by = "TAXANAME1") %>%
  group_by(Site, htype_collapsed) %>%
  summarise(
    across(
      all_of(trait_cols),
      ~ sum(.x * Biomass_m2, na.rm = TRUE),
      .names = "traitbio_{.col}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(cols      = starts_with("traitbio_"),
               names_to  = "Trait",
               values_to = "trait_biomass") %>%
  mutate(
    Trait      = gsub("^traitbio_", "", Trait),
    TraitGroup = case_when(
      startsWith(Trait, "Feed_")  ~ "Feeding Group",
      startsWith(Trait, "Habit_") ~ "Habit"
    ),
    TraitLabel = label_trait(Trait)
  ) %>%
  filter(!is.na(TraitGroup)) %>% 
  filter(!TraitLabel== "Crawler")






  

write.csv(
  traitbio_collapsed %>%
    dplyr::select(Site, htype_collapsed, TraitGroup, TraitLabel, trait_biomass),
  "outputs/H4.1_trait_biomass_collapsed.csv",
  row.names = FALSE
)


trait_raw_debug<- traitbio_collapsed %>% 
  group_by(Site, htype_collapsed,TraitGroup) %>% 
  summarise(group_biomass = sum(trait_biomass, na.rm=TRUE), .groups='drop') %>% 
  left_join(
    raw_collapsed_site %>% 
     group_by(Site, htype_collapsed) %>% 
  summarise(raw_total_biomass = sum(Biomass_m2, na.rm =TRUE), .groups='drop'), 
    by = c("Site", "htype_collapsed"))
  

write.csv(trait_raw_debug, "outputs/collapsedbiomassdebug.csv")

# ---- Part B: replicate-level (all 6 habitat types) ----
traitbio_habtype <- all_raw %>%
  filter(!TAXANAME1 %in% document_only) %>%
  left_join(trait_matrix, by = "TAXANAME1") %>%
  group_by(Site, Sample, htype) %>%
  summarise(
    across(
      all_of(trait_cols),
      ~ sum(.x * Biomass_m2, na.rm = TRUE),
      .names = "traitbio_{.col}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(cols      = starts_with("traitbio_"),
               names_to  = "Trait",
               values_to = "trait_biomass") %>%
  mutate(
    Trait      = gsub("^traitbio_", "", Trait),
    TraitGroup = case_when(
      startsWith(Trait, "Feed_")  ~ "Feeding Group",
      startsWith(Trait, "Habit_") ~ "Habit"
    ),
    TraitLabel = label_trait(Trait),
    htype      = factor(htype,
                        levels = c("ub", "da", "rts", "sng", "veg", "riffle"))
  ) %>%
  filter(!is.na(TraitGroup))



write.csv(
  traitbio_habtype %>%
    dplyr::select(Site, Sample, htype, TraitGroup, TraitLabel, trait_biomass),
  "outputs/H4.1_trait_biomass_by_habtype.csv",
  row.names = FALSE
)

# ============================================================
# PART A: RIFFLE vs. DEPOSITIONAL — SITE-LEVEL
# Stats: Wilcoxon rank-sum. Letters: a / b.
# ============================================================

plot_collapsed <- function(trait_label, trait_group) {

  df <- traitbio_collapsed %>%
    filter(TraitLabel == trait_label, !is.na(trait_biomass))
  if (nrow(df) == 0 || n_distinct(df$htype_collapsed) < 2) {
    message("Skipping (insufficient data): ", trait_label); return(invisible(NULL))
  }

  n_lab <- df %>% group_by(htype_collapsed) %>% summarise(n = n(), .groups = "drop")

  wx     <- wilcox.test(trait_biomass ~ htype_collapsed, data = df, exact = FALSE)
  wx_lab <- sprintf("Wilcoxon rank-sum: W = %.0f, p = %.3f",
                    wx$statistic, wx$p.value)

  means2 <- df %>% group_by(htype_collapsed) %>%
    summarise(m = mean(trait_biomass, na.rm = TRUE), .groups = "drop") %>%
    arrange(m)
  grp_letters <- if (wx$p.value < 0.05)
    setNames(c("a", "b"), as.character(means2$htype_collapsed)) else
    setNames(c("a", "a"), as.character(means2$htype_collapsed))

  letter_pos <- df %>% group_by(htype_collapsed) %>%
    summarise(top = max(trait_biomass, na.rm = TRUE), .groups = "drop") %>%
    mutate(letter = grp_letters[as.character(htype_collapsed)],
           y      = top * 1.6)

  p <- ggplot(df, aes(htype_collapsed, trait_biomass, fill = htype_collapsed)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, width = 0.5, alpha = 0.85) +
    geom_point(stat = "summary", fun = mean,
               shape = 23, size = 3, fill = "white", color = "black") +
    scale_fill_manual(values = collapsed_colors) +
    scale_y_continuous(trans = "log1p", labels = scales::comma)+
    scale_x_discrete(labels = function(x) {
      n <- n_lab$n[match(x, as.character(n_lab$htype_collapsed))]
      paste0(collapsed_names[x], "\nn = ", n)
    }) +
    labs(title    = paste(trait_group, "\u2014", trait_label),
         subtitle = wx_lab,
         x        = "Habitat type",
         y        = expression("Log(x+1) Trait-group biomass  (mg DM m"^-2*")")) +
    theme_bw(base_size = 20) +
    theme(
      legend.position  = "none",
      title            = element_text(size = 16),
      axis.title       = element_text(size = 20),
      axis.text.x      = element_text(size = 20),
      axis.text.y      = element_text(size = 25),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    geom_text(data = letter_pos, aes(htype_collapsed, y, label = letter),
              inherit.aes = FALSE, size = 10, fontface = "bold", vjust = 0)

  fname <- tolower(gsub("[^a-zA-Z0-9]", "_", trait_label))
  ggsave(paste0("outputs/H4_4.15/H4.1_plot_traitbio_", fname, "_collapsed.png"),
         p, width = 10, height = 6, dpi = 300, bg = "transparent")
  invisible(p)
}

cat("--- Part A: collapsed trait-group biomass plots ---\n")
for (tl in feeding_groups) plot_collapsed(tl, "Feeding Group")
for (tl in habits)         plot_collapsed(tl, "Habit")

# ============================================================
# PART B: ALL 6 HABITAT TYPES — REPLICATE-LEVEL
# Stats: Kruskal-Wallis + Dunn (BH). Compact letter displays.
# ============================================================

# ---- shared internals ----
.traitbio_prep <- function(trait_label, fixed_levels = NULL) {
  df <- traitbio_habtype %>%
    filter(TraitLabel == trait_label, !is.na(trait_biomass)) %>%
    # --- collapse replicates: one mean per Site x htype ---
    group_by(Site, htype, TraitLabel) %>%
    summarise(trait_biomass = mean(trait_biomass, na.rm = TRUE), .groups = "drop")
  
  if (nrow(df) == 0 || n_distinct(df$htype) < 2) return(NULL)
  
  ord <- if (is.null(fixed_levels)) {
    df %>% group_by(htype) %>%
      summarise(m = median(trait_biomass, na.rm = TRUE), .groups = "drop") %>%  # mean -> median
      arrange(m) %>% pull(htype) %>% as.character()
  } else as.character(fixed_levels)
  
  df    <- df %>% mutate(htype = factor(htype, levels = ord))
  n_lab <- df %>% group_by(htype) %>% summarise(n = n(), .groups = "drop") %>%
    mutate(htype = factor(htype, levels = ord))
  
  kw <- kruskal.test(trait_biomass ~ htype, data = df)
  dunn_full <- df %>% dunn_test(trait_biomass ~ htype, p.adjust.method = "BH")
  
  top_all <- max(df$trait_biomass, na.rm = TRUE)
  
  if (kw$p.value < 0.05) {
    pv <- dunn_full$p.adj
    names(pv) <- paste(dunn_full$group1, dunn_full$group2, sep = "-")
    grp_letters <- multcompLetters(pv, threshold = 0.05)$Letters
    
    letter_pos <- data.frame(htype  = names(grp_letters),
                             letter = unname(grp_letters),
                             stringsAsFactors = FALSE) %>%
      mutate(y = top_all * 1.15,
             htype = factor(htype, levels = levels(df$htype)))
  } else {
    # non-significant overall KW: no letters to display
    letter_pos <- data.frame(htype = character(0), letter = character(0),
                             y = numeric(0))
  }
  
  list(df = df, n_lab = n_lab, kw = kw, letter_pos = letter_pos)
}

# ---- panel plot (for cowplot assembly) ----
plot_traitbio_panel <- function(trait_label, trait_group, fixed_levels = NULL) {
  prep <- .traitbio_prep(trait_label, fixed_levels)
  if (is.null(prep)) return(NULL)
  df <- prep$df; n_lab <- prep$n_lab; kw <- prep$kw; letter_pos <- prep$letter_pos
  kw_lab <- sprintf("KW: \u03C7\u00B2 = %.2f, df = %d, p = %.3f",
                    kw$statistic, kw$parameter, kw$p.value)
  
  ggplot(df, aes(htype, trait_biomass, fill = htype)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.2, width = 0.6, alpha = 0.85) +
    scale_y_continuous(
      trans  = "log1p",
      labels = scales::comma,
      n.breaks = 4)+
    scale_fill_manual(values = hab_colors) +
    scale_x_discrete(labels = function(x) {
      n <- n_lab$n[match(x, as.character(n_lab$htype))]
      paste0(hab_full_names[x], "\nn = ", n)
    }) +
    labs(title    = trait_label,
        ######## #subtitle = kw_lab,####### add back for the manuscript plots
         x = NULL, y = NULL) +
    geom_text(data = letter_pos, aes(htype, y, label = letter),
              inherit.aes = FALSE, size = 8, fontface = "bold",
              vjust = 0, color = "black") +
    theme_classic(base_size = 12) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(size = 20, face = "bold", color = "black"),
      plot.subtitle      = element_text(size = 20, color = "black"),
      axis.text.x        = element_text(size =22, angle =0, hjust = 1,
                                        vjust = 1, lineheight = 0.9,
                                        color = "black"),
      axis.text.y        = element_text(size = 22, color = "black"),
      axis.title.x       = element_text(size = 22, color = "black"),
      axis.title.y       = element_text(size = 20, color = "black"),
      # L-shaped frame: keep left + bottom axis lines only
      axis.line          = element_line(colour = "black"),
      axis.line.x.top    = element_blank(),
      axis.line.y.right  = element_blank(),
      panel.border       = element_blank(),
      panel.background   = element_rect(fill = "transparent", colour = NA),
      plot.background    = element_rect(fill = "transparent", colour = NA),
      panel.grid.major   = element_line(colour = "grey85"),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(t = 8, r = 10, b = 8, l = 8)
    )
}

# ---- shared axis labels ----
shared_y <- ggdraw() +
  draw_label(expression(" Log(x+1) Trait-group biomass (mg DM m"^-2*")"),
             fontface = "bold", size = 25, angle = 90)
shared_x <- ggdraw() + draw_label("Habitat type", fontface = "bold", size =20)

make_combined <- function(plot_list, n_cols, title_text = NULL) {
  grid   <- plot_grid(plotlist = plot_list, ncol = n_cols,
                      align = "hv", axis = "tblr")
  
  # Two "Habitat type" labels, one per column, aligned to the grid (not the y-label gutter)
  shared_x_left  <- ggdraw() + draw_label("Habitat type", fontface = "bold", size = 20)
  shared_x_right <- ggdraw() + draw_label("Habitat type", fontface = "bold", size = 20)
  shared_x_row   <- plot_grid(shared_x_left, shared_x_right, ncol = 2)
  
  grid_x <- plot_grid(grid, shared_x_row, ncol = 1, rel_heights = c(1, 0.03))
  body   <- plot_grid(shared_y, grid_x, ncol = 2, rel_widths = c(0.04, 1))
  
  if (is.null(title_text)) return(body)
  title  <- ggdraw() + draw_label(title_text, fontface = "bold", size = 16)
  plot_grid(title, body, ncol = 1, rel_heights = c(0.04, 1))
}

cat("--- Part B/C: per-habitat trait-group biomass panels ---\n")
ffg_plots   <- purrr::compact(lapply(feeding_groups, plot_traitbio_panel,
                                     trait_group = "Feeding Group"))
habit_plots <- purrr::compact(lapply(habits,         plot_traitbio_panel,
                                     trait_group = "Habit"))



ffg_combined <- make_combined(ffg_plots, n_cols = 2)
ggsave("outputs/H4.1_traitbio_FFG_combined.png",
       ffg_combined, width = 14, height = 20, dpi = 600, bg = "transparent")

ffg_combined

habit_combined <- make_combined(habit_plots, n_cols = 2)
ggsave("outputs/H4.1_traitbio_Habit_combined.png",
       habit_combined, width = 14, height = 21, dpi = 600, bg = "transparent")
habit_combined

# ============================================================
# PART D: KRUSKAL-WALLIS + DUNN RESULTS TABLE
# ============================================================


library(flextable)
library(officer)

sig_stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))

fmt_p <- function(p) ifelse(is.na(p), "",
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

kw_dunn_rows <- function(trait_label, trait_group) {
  df <- traitbio_habtype %>%
    filter(TraitLabel == trait_label, !is.na(trait_biomass)) %>%
    group_by(Site, htype, TraitLabel) %>%
    summarise(trait_biomass = mean(trait_biomass, na.rm = TRUE), .groups = "drop")
  if (nrow(df) == 0 || n_distinct(df$htype) < 2) return(NULL)
  
  kw   <- kruskal.test(trait_biomass ~ htype, data = df)
  dunn <- df %>% dunn_test(trait_biomass ~ htype, p.adjust.method = "BH")
  
  kw_row <- data.frame(
    TraitGroup = trait_group,
    Trait      = trait_label,
    Comparison = "Overall (Kruskal-Wallis)",
    Statistic  = unname(kw$statistic),
    df         = unname(kw$parameter),
    p          = kw$p.value,
    p_adj      = NA_real_,
    Sig        = sig_stars(kw$p.value),
    row_type   = "kw",
    stringsAsFactors = FALSE
  )
  
  dunn_rows <- data.frame(
    TraitGroup = trait_group,
    Trait      = trait_label,
    Comparison = paste(hab_full_names[dunn$group1], "vs",
                       hab_full_names[dunn$group2]),
    Statistic  = dunn$statistic,
    df         = NA_integer_,
    p          = dunn$p,
    p_adj      = dunn$p.adj,
    Sig        = sig_stars(dunn$p.adj),
    row_type   = "dunn",
    stringsAsFactors = FALSE
  )
  
  rbind(kw_row, dunn_rows)
}

ffg_stats   <- bind_rows(lapply(feeding_groups, kw_dunn_rows,
                                trait_group = "Feeding Group"))
habit_stats <- bind_rows(lapply(habits,         kw_dunn_rows,
                                trait_group = "Habit"))
all_stats   <- bind_rows(ffg_stats, habit_stats)

write.csv(all_stats,
          "outputs/H4_6_24/H4.1_traitbio_KW_Dunn_results.csv",
          row.names = FALSE)

make_kw_dunn_flextable <- function(stats_df, sig_only = TRUE) {

  tbl <- stats_df %>%
    filter(row_type == "kw" |
           (!sig_only) |
           (row_type == "dunn" & !is.na(p_adj) & p_adj < 0.05))

  display <- tbl %>%
    mutate(
      Stat = sprintf("%.2f", Statistic),
      dff  = ifelse(is.na(df), "", as.character(df)),
      pval = ifelse(row_type == "kw", fmt_p(p), fmt_p(p_adj))
    ) %>%
    dplyr::select(Trait, Comparison, Stat, dff, pval, Sig)

  kw_rows   <- which(tbl$row_type == "kw")
  dunn_rows <- which(tbl$row_type == "dunn")

  flextable(display) %>%
    set_header_labels(Trait = "Trait", Comparison = "Comparison",
                      Stat = "Z / \u03C7\u00B2", dff = "df",
                      pval = "p", Sig = "Sig") %>%
    merge_v(j = "Trait") %>%
    valign(j = "Trait", valign = "top") %>%
    bold(i = kw_rows, part = "body") %>%
    italic(i = dunn_rows, j = "Comparison", part = "body") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    bold(part = "header") %>%
    align(part = "all", align = "center") %>%
    align(j = c("Trait", "Comparison"), align = "left", part = "body") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(color = "black", width = 0.75)) %>%
    hline_bottom(part = "body",   border = fp_border(color = "black", width = 1.5)) %>%
    width(j = "Trait",      width = 1.7) %>%
    width(j = "Comparison", width = 1.7) %>%
    width(j = "Stat",       width = 0.8) %>%
    width(j = "dff",        width = 0.4) %>%
    width(j = "pval",       width = 0.8) %>%
    width(j = "Sig",        width = 0.5)
}

ft_ffg   <- make_kw_dunn_flextable(ffg_stats)
ft_habit <- make_kw_dunn_flextable(habit_stats)

cap_ffg <- paste0(
  "Table X. Kruskal-Wallis and Dunn post-hoc (BH-adjusted) results for ",
  "trait-group biomass (feeding groups) across the six habitat types. ",
  "Trait-group biomass was calculated as the summed dry mass of all taxa ",
  "assigned to each binary trait category within a replicate (mg DM m\u207B\u00B2). ",
  "Bold rows: overall KW test (\u03C7\u00B2, df, p). ",
  "Italic rows: significant pairwise Dunn comparisons (Z, BH-adjusted p)."
)

doc <- read_docx() %>%
  body_add_par(cap_ffg, style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft_ffg) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table X. Habits (same structure as above).", style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft_habit)

print(doc, target = "outputs/H4_6_24/H4.1_traitbio_KW_Dunn_table.docx")

try({
  save_as_image(ft_ffg,   "outputs/H4_6_24/H4.1_traitbio_KW_Dunn_FFG.png",
                zoom = 3, expand = 10)
  save_as_image(ft_habit, "outputs/H4_6_24/H4.1_traitbio_KW_Dunn_Habit.png",
                zoom = 3, expand = 10)
}, silent = TRUE)

cat("All outputs written to outputs/H4_6_24/\n")

# ============================================================
# PART D (addition): KW + DUNN TABLE — TOTAL BIOMASS ONLY
# ============================================================

bmas_by_hab <- all_raw %>% 
  group_by(htype, Sample, Site) %>% 
  rep_bio <- (summarise(sum(Biomass_m2) .groups ="drop"))) %>% 
  group_by(htype, Site) %>% 
  hab_biomass <- (summarise(mean(rep_bio) ,groups ="drop")))

bmas_by_hab


totalbio_stats <- kw_dunn_rows("Total Biomass", "Total Biomass")

write.csv(totalbio_stats,
          "outputs/H4_6_24/H4.1_totalbiomass_KW_Dunn_results.csv",
          row.naes = FALSE)

# sig_only = FALSE so you see ALL pairwise habitat comparisons,
# not just the significant ones (unlike the FFG/Habit tables)
ft_totalbio <- make_kw_dunn_flextable(totalbio_stats, sig_only = FALSE)

cap_totalbio <- paste0(
  "Table X. Kruskal-Wallis and Dunn post-hoc (BH-adjusted) results for total ",
  "macroinvertebrate community biomass across the six habitat types. Total ",
  "biomass was calculated as the summed dry mass of all taxa within a replicate ",
  "(mg DM m\u207B\u00B2), excluding taxa retained at coarse taxonomic resolution for ",
  "documentation only. Bold row: overall KW test (\u03C7\u00B2, df, p). Italic rows: ",
  "pairwise Dunn comparisons (Z, BH-adjusted p)."
)

doc_totalbio <- read_docx() %>%
  body_add_par(cap_totalbio, style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft_totalbio)

print(doc_totalbio, target = "outputs/H4_6_24/H4.1_totalbiomass_KW_Dunn_table.docx")

try(save_as_image(ft_totalbio, "outputs/H4_6_24/H4.1_totalbiomass_KW_Dunn.png",
                  zoom = 3, expand = 10), silent = TRUE)

ft_totalbio


# ---- end of script ----



#results stats

summary_traits <- traitbio_habtype %>% 
  group_by(TraitGroup, TraitLabel, htype) %>% 
  summarise(
    mean=mean(trait_biomass, na.rm=TRUE),
    median=median(trait_biomass, na.rm=TRUE),
    n = n(),
    .groups ="drop"
  ))
summary_traits

write.csv(summary_traits,"outputs/summary_traits7.2.26.csv")


# ============================================================
# SECTION 7B: TOTAL COMMUNITY BIOMASS (all taxa, unweighted by trait)
# Added as an extra "trait" so it can reuse .traitbio_prep() /
# plot_traitbio_panel() unchanged — flows into the FFG panel grid.
# ============================================================
totalbio_htype2 <- raw_collapsed %>% 
  filter(!TAXANAME1 %in% document_only) %>% 
  group_by(Site, Sample, htype_collapsed) %>% 
  summarise(total_biomass = sum(Biomass_m2, na.rm= TRUE), .groups= "drop") %>% 
  group_by(Site, htype_collapsed) %>% 
  summarise(mean_biomass= mean(total_biomass, na.rm= TRUE), .groups= "drop")


#produce median values
median_rif <-totalbio_htype2 %>% 
  group_by(htype_collapsed) %>% 
  summarise(median= median(mean_biomass, na.rm = TRUE), .groups= "drop")
median_rif

totaldens_htype2 <- raw_collapsed %>% 
  filter(!TAXANAME1 %in% document_only) %>% 
  group_by(Site, Sample, htype_collapsed) %>% 
  summarise(total_density = sum(Density_m2, na.rm= TRUE), .groups= "drop") %>% 
  group_by(Site, htype_collapsed) %>% 
  summarise(mean_density= mean(total_density, na.rm= TRUE), .groups= "drop")

median_rif_dens <-totaldens_htype2 %>% 
  group_by(htype_collapsed) %>% 
  summarise(median= median(mean_density, na.rm = TRUE), .groups= "drop")
median_rif_dens

#density KW test 
KW_dens <- kruskal.test(mean_density ~ htype_collapsed, data=totaldens_htype2)
KW_dens

totalbio_habtype <- all_raw %>% 
  filter(!TAXANAME1 %in% document_only) %>%
  group_by(Site, Sample, htype) %>%
  summarise(trait_biomass = sum(Biomass_m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Trait      = "Total_Biomass",
    TraitGroup = "Feeding Group",      # bucket it with the FFG panels
    TraitLabel = "Total Biomass",
    htype      = factor(htype, levels = c("ub", "da", "rts", "veg", "riffle"))
  )

traitbio_habtype <- bind_rows(traitbio_habtype, totalbio_habtype)

write.csv(
  totalbio_habtype %>%
    dplyr::select(Site, Sample, htype, TraitGroup, TraitLabel, trait_biomass),
  "outputs/H4.1_total_biomass_by_habtype.csv",
  row.names = FALSE
)

feeding_groups_plus <- c(feeding_groups, "Total Biomass")

ffg_plots   <- purrr::compact(lapply(feeding_groups_plus, plot_traitbio_panel,
                                     trait_group = "Feeding Group"))


ffg_combined <- make_combined(ffg_plots, n_cols = 2)
ggsave("outputs/H4_6_24/H4.1_traitbio_FFG_combined.png",
       ffg_combined, width = 14, height = 17, dpi = 600, bg = "transparent")

ffg_combined

habits_plus <- c(habits, "Total Biomass")

habit_plots <- purrr::compact(lapply(habits_plus, plot_traitbio_panel,
                                     trait_group = "Habit"))
habit_combined <- make_combined(habit_plots, n_col=2)

habit_combined

ggsave("outputs/H4_6_24/H4.1_traitbio_habit_combined.png",
              habit_combined, width = 14, height = 17, dpi = 600, bg = "transparent")



########################################################
#biomass boxplot 

kw_rif_nf <- kruskal.test(mean_biomass~ htype_collapsed, data= totalbio_htype2)
kw_rif_nf

kw_label <- sprintf("Kruskal-Wallis, χ² = %.2f, df = %d, p = %.3f",
                    kw_rif_nf$statistic, kw_rif_nf$parameter, kw_rif_nf$p.value)



htype2_box <- ggplot(totalbio_htype2, aes(htype_collapsed, mean_biomass, fill=htype_collapsed)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.2, width = 0.6, alpha = 0.85) +
  scale_y_continuous(trans = "log1p", labels = scales::comma,
                     expand = expansion(mult = c(0.05, 0.12)),
                     ) +
  scale_fill_manual(values = collapsed_colors)+
  scale_x_discrete(labels = htyp2_names) +
  labs(y= "Log+1 mean site biomass (mg/m2)",
       subtitle = kw_label)+
  theme_bw(base_size = 16) +
  theme_classic(base_size = 16) +
  theme(
    axis.title.y = element_text(size=14),
    axis.title.x = element_blank(),
    legend.position  = "none",
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 14),
    axis.text.x      = element_text(size = 14, angle = 0, hjust = 0.5),
    # L-shaped frame: keep left + bottom axis lines only
    axis.line        = element_line(colour = "black"),
    axis.line.x.top  = element_blank(),
    axis.line.y.right = element_blank(),
    panel.border     = element_blank(),          # remove the full box from theme_bw
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank())
htype2_box  

ggsave("outputs/H4_6_24/log_biomass_boxplot.png", 
       htype2_box, width = 6, height = 6, dpi = 600, bg = "transparent")


# ============================================================
# MEDIAN TRAIT-GROUP BIOMASS — FEEDING GROUPS, ALL 6 HABITAT TYPES
# ============================================================

ffg_median_habtype <- traitbio_habtype %>%
  filter(TraitGroup == "Feeding Group", !is.na(trait_biomass)) %>%
  group_by(TraitLabel, htype) %>%
  summarise(
    median_biomass = median(trait_biomass, na.rm = TRUE),
    mean_biomass   = mean(trait_biomass, na.rm = TRUE),
    n              = n(),
    .groups = "drop"
  ) %>%
  arrange(TraitLabel, htype)

print(ffg_median_habtype, n=1000)

# Wide version — one row per FFG, one column per habitat, for a quick
# read (median mg DM m^-2 per feeding group x habitat)
ffg_median_wide <- ffg_median_habtype %>%
  dplyr::select(TraitLabel, htype, median_biomass) %>%
  pivot_wider(names_from = htype, values_from = median_biomass)

ffg_median_wide

write.csv(ffg_median_habtype,
          "outputs/H4_6_24/H4.1_FFG_median_biomass_by_habtype.csv",
          row.names = FALSE)
# ============================================================
# COMPOSITE MEDIAN TRAIT-GROUP BIOMASS TABLE — FFG + HABIT COMBINED
# Excludes Piercer/Herbivore and Parasite. Reflects the single-trait
# revision: each taxon's biomass now contributes to exactly one FFG
# and one Habit category, so these medians are no longer inflated by
# the old multi-flag double-counting.
# ============================================================

hab_levels <- c("ub", "da", "rts", "veg", "riffle")

composite_median <- traitbio_habtype %>%
  filter(
    !TraitLabel %in% c("Total Biomass", "Piercer/Herbivore", "Parasite", "Habit_Crawler"),
    !is.na(trait_biomass)
  ) %>%
  mutate(htype = factor(htype, levels = hab_levels)) %>%
  filter(!is.na(htype)) %>%                      # drops snag + any stray codes
  group_by(TraitGroup, TraitLabel, htype, .drop = FALSE) %>%
  summarise(
    median_biomass = median(trait_biomass, na.rm = TRUE),
    n              = n(),
    .groups = "drop"
  ) %>%
  mutate(cell = ifelse(n == 0, "\u2014", sprintf("%.2f (n=%d)", median_biomass, n)))

composite_wide <- composite_median %>%
  dplyr::select(TraitGroup, TraitLabel, htype, cell) %>%
  pivot_wider(names_from = htype, values_from = cell, names_expand = TRUE) %>%
  dplyr::select(TraitGroup, TraitLabel, all_of(hab_levels)) %>%
  arrange(TraitGroup, TraitLabel)

make_composite_median_flextable <- function(df) {
  hab_short <- c(ub = "Undercut\nBank", da = "Debris\nAccum.",
                 rts = "Submerged\nRoots",
                 veg = "Vegetation", riffle = "Riffle")
  
  trait_rows <- which(df$TraitGroup == "Feeding Group")
  
  ft <- flextable(df) %>%
    set_header_labels(
      TraitGroup = "Trait Group", TraitLabel = "Trait",
      ub = hab_short[["ub"]], da = hab_short[["da"]], rts = hab_short[["rts"]],
      veg = hab_short[["veg"]], riffle = hab_short[["riffle"]]
    ) %>%
    merge_v(j = "TraitGroup") %>%
    valign(j = "TraitGroup", valign = "top") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>%
    align(part = "all", align = "center") %>%
    align(j = c("TraitGroup", "TraitLabel"), align = "left", part = "body") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(color = "black", width = 0.75)) %>%
    hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
    width(j = "TraitGroup", width = 0.9) %>%
    width(j = "TraitLabel", width = 1.2) %>%
    width(j = hab_levels, width = 0.95)
  
  if (length(trait_rows) > 0 && max(trait_rows) < nrow(df)) {
    ft <- hline(ft, i = max(trait_rows),
                border = fp_border(color = "black", width = 0.75))
  }
  ft
}

ft_composite_median <- make_composite_median_flextable(composite_wide)
ft_composite_median

cap_composite <- paste0(
  "Table X. Median trait-group biomass (mg DM m\u207B\u00B2; n = number of replicates in ",
  "parentheses) for functional feeding groups and habits across the six mesohabitat ",
  "types. Piercer/Herbivore and Parasite feeding groups are excluded due to low ",
  "representation. Each taxon was assigned exactly one feeding group and one habit ",
  "trait (single-trait assignment; see Trait_Generation.R), so biomass is not ",
  "double-counted across trait categories."
)

doc_composite <- read_docx() %>%
  body_add_par(cap_composite, style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft_composite_median)

print(doc_composite, target = "outputs/H4_6_24/H4.1_composite_median_biomass_table.docx")

try(save_as_image(ft_composite_median,
                  "outputs/H4_6_24/H4.1_composite_median_biomass_table.png",
                  zoom = 3, expand = 10), silent = TRUE)

write.csv(composite_wide,
          "outputs/H4_6_24/H4.1_composite_median_biomass_table.csv",
          row.names = FALSE)

