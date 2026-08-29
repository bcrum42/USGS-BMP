  #============================================================
  # SCRIPT: H2_Biomass_NMDS.R
  # PURPOSE: Calculate replicate-level biomass for dep hab and
  #          riffle macroinvertebrates for downstream analysis
  # H2: Macroinvertebrate richness in depositional habitats and
  #     riffles greatest in streams with most dep hab diversity
  # Last updated: 2026-03-30
  # Biomass code written using code created from Sergio Sabat-Bonilla
  # code written with the help of Claude AI 
  # ============================================================
  
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
  
  ###########################set wd####################################
  setwd("G:/Shared drives/Aquatic Entomology Lab/Projects/USGS BMP/Year 1_SHEN_AG/Brice's Corner/BC thesis materials/BC Thesis")
  
  ##########################load data ###########################
  #########
  
  
  
  rif.data1 <- read.csv("G:/Shared drives/Aquatic Entomology Lab/Projects/USGS BMP/Year 1_SHEN_AG/project data/Year 1 Shen Recorded Macros/R_Scripts_&_Outputs/Chapter 2 - Typology 1 analysis/Taxa_Abundance_Length.csv")
  
  dep.data  <- read.csv("outputs/typ1_hab_specific_data_3.10.26.csv")
  
  taxa_info <- read.csv("Taxa_Length-Mass_Information_BC.csv")
  
  drop_table <- read.csv("drop_table_3.2.26.csv", check.names = FALSE)
  
  
  ############################################
  #Riffle data cleaning 
  #################################################################
  
  # pivot drop table to long and identify reps to remove
  drop_long <- drop_table %>%
    pivot_longer(
      cols      = `1`:`5`,
      names_to  = "Sample",
      values_to = "drop"
    ) %>%
    filter(drop == "x") %>%
    mutate(Sample = as.integer(Sample)) %>%
    select(Site, Sample)
  
  # drop unusable replicates and excluded sites
  rif.data1$Sample <- as.integer(rif.data1$Sample)
  
  rif.data2 <- rif.data1 %>%
    anti_join(drop_long, by = c("Site", "Sample")) %>%
    filter(!Site %in% c("WEST", "MIDD", "MLBK", "MLRH", "BULL")) %>%
    mutate(
      loc = "riffle",
      across(c(Order, Family, Genus), trimws)
    )
  
  # fix misspelled taxa
  rif.data2 <- rif.data2 %>%
    mutate(
      Genus = case_when(
        Genus == "Hygenius brevistylus" ~ "Hagenius brevistylus",
        Genus == "Percomaini/Psychoda"  ~ "Psychoda",
        Genus == "Trianoides"           ~ "Triaenodes",
        TRUE                            ~ Genus
      )
    ) %>%
    filter(Family != "Collembola")

#fix issues with corydalus
 rif.data2 <- rif.data2 %>%
    mutate(LIFE_STAGE = case_when(
      Genus =="Corydalus"~"L",
      TRUE ~ LIFE_STAGE))
  
  # unique taxa reference table
  rif_taxa <- rif.data2 %>%
    distinct(Order, Family, Genus) %>%
    arrange(Order)
###############################
#write clean riffle habitat code 

write.csv(rif.data2, "clean_riffle_abundance_3.30.26.csv")
###############################
  

########################################
  #dep data cleaning
########################################
  
  dep.data1 <- dep.data %>%
    # rename CASH to EIDS to match riffle site codes
    mutate(
      Site = replace(Site, Site == "CASH", "EIDS"),
      across(c(Order, Family, Genus), trimws))

  #FILTER SITE THAT IS NOT PRESENT IN RIFFLE DATA 
  dep.data1 <- dep.data1 %>% 
    filter(!(Site == "SPOU" & Sample == 4))

  
  # fix misspelled taxa and metadata errors
  dep.data1 <- dep.data1 %>%
    mutate(
      # genus fixes
      Genus = case_when(
        Genus == "Dinetus"               ~ "Dineutus",
        Genus == "Ceraptogoninae*"        ~ "Ceratopogoniae",
        Genus == "Stenacron_Maccaffertium"~ "Stenacron/Maccaffertium",
        Genus == "Sperchopis"             ~ "Sperchopsis",
        Genus == "lype diversa"           ~ "lype",
        TRUE                              ~ Genus
      ),
      # family fixes
      Family = case_when(
        Family == "Leptohyidae"      ~ "Leptohyphidae",
        Family == "Lecutridae"       ~ "Leuctridae",
        Family == "Helicophsychidae" ~ "Helicopsychidae",
        Family == "Chironominae"     ~ "Chironomidae",
        Family == "Velidae"          ~ "Veliidae",
        Family == "Vellidae"         ~ "Veliidae",
        Order == "Trichoptera" & 
          Genus == "Neureclipsis"    ~ "Polycentropodidae",
        TRUE                         ~ Family
      ),
      # order fixes
      Order = case_when(
        Order == "Coleoptea"      ~ "Coleoptera",
        Order == "Simuliidae"     ~ "Diptera",
        Order == "Empididae"      ~ "Diptera",
        Order == "Hydroptilidae"  ~ "Trichoptera",
        Order == "Philopotamidae" ~ "Trichoptera",
        Order == "Rhyacophilidae" ~ "Trichoptera",
        Order == "Tabanidae"      ~ "Diptera",
        Order == "trichoptera"    ~ "Trichoptera",
        Order == "Bivalva"        ~ "Bivalvia",
        Family == "Veliidae" & 
          Genus == "Steinovelia"  ~ "Hemiptera",
        TRUE                      ~ Order
      )
    ) 
  
  # ============================================================
  # SECTION 4: TAXON NAME AND LIFE STAGE STANDARDIZATION
  # ============================================================
  
  # --- 4A: create TAXANAME and TAXANAME1 columns ---
  dep.data1 <- dep.data1 %>%
    mutate(TAXANAME = Genus)
  
  dep.data1.1 <- dep.data1 %>%
    mutate(
      TAXANAME1 = tolower(trimws(TAXANAME)),
      # strip parenthetical life stage info from name e.g. "Hydraena(A)"
      TAXANAME1 = gsub("\\s*\\([^)]*\\)", "", TAXANAME1)
    )
  
  # --- 4B: taxon name corrections (typos, fallback names) ---
  dep.data1.1 <- dep.data1.1 %>%
    mutate(TAXANAME1 = case_when(
      # typos / transposed letters
      TAXANAME1 == "simulium/stegopterna" ~ "simulium",
      TAXANAME1 == "vivipariadae"      ~ "viviparidae",
      TAXANAME1 == "psuedolimnophila"  ~ "pseudolimnophila",
      TAXANAME1 == "peltodtyes"        ~ "peltodytes",
      TAXANAME1 == "salididae"         ~ "saldidae",
      TAXANAME1 == "sperchopis"        ~ "sperchopsis",
      TAXANAME1 == "crytpolabis"       ~ "cryptolabis",
      TAXANAME1 == "brachydentridae"   ~ "brachycentridae",
      TAXANAME1 == "stenonoma"         ~ "stenonema",
      TAXANAME1 == "orchotrichia"      ~ "ochrotrichia",
      TAXANAME1 == "anchytarsis"       ~ "anchytarsus",
      TAXANAME1 == "argyra"            ~ "argya",
      TAXANAME1 == "dinetus"           ~ "dineutus",
      TAXANAME1 == "psididiidae"       ~ "pisidiidae",
      TAXANAME1 == "pisididiidae"      ~ "pisidiidae",
      TAXANAME1 == "psididiidae"       ~ "pisidiidae",
      TAXANAME1 == "bithyniidae"       ~ "planorbidae",
      TAXANAME1 == "ceratopogoninae"   ~ "ceratopogonidae",
      TAXANAME1 == "lype diversa"      ~ "lype",
      # fallback to family/order level equations
      TAXANAME1 == "libellula"                ~ "libellulidae",
      TAXANAME1 == "nemoura"                  ~ "nemouridae",
      TAXANAME1 == "stenacron/maccaffertium"  ~ "stenacron_maccaffertium",
      TAXANAME1 == "dannella"                 ~ "ephemerellidae",
      TAXANAME1 == "hexagenia"                ~ "ephemeridae",
      TAXANAME1 == "nectopsyche"              ~ "leptoceridae",
      TAXANAME1 == "pelocoris"                ~ "pelocoris",
      TRUE                                    ~ TAXANAME1
    ))
  
  # --- 4C: order corrections that affect life stage ---
  dep.data1.1 <- dep.data1.1 %>%
    mutate(Order = case_when(
      TAXANAME1 == "hesperocorixa" ~ "Hemiptera",
      TAXANAME1 == "hydrobius"     ~ "Coleoptera",
      TRUE                         ~ Order
    ))
  
  # --- 4D: life stage assignment ---
  
  # taxa excluded from biomass (no equations / non-aquatic)
  non_insect_taxa <- c(
    "lirceus", "pleuroceridae", "oligochaeta",
    "physidae", "viviparidae", "planorbidae",
    "corbicula", "lymnaeidae", "bithyniidae",
    "ancylidae", "sphaeriidae", "psididiidae",
    "crangonyx", "copepoda", "ostracod", "daphnia",
    "hydracarina", "hirundinea", "annelida",
    "hydrobiidae", "collembola", "pisidiidae"
  )
  
  junk_taxa <- c("terrestrial", "terrestrial?")
  
  hemimetabolous_orders <- c(
    "Odonata", "Plecoptera", "Ephemeroptera", "Hemiptera"
  )
  
  # taxa only collected as adults
  adult_only_taxa <- c(
    "aluetes", "gymnochthebius", "listronotus",
    "plateumaris", "prasocuris", "sphenophorus"
  )
  
  # hemiptera where only adult equation exists in taxa_info
  hemiptera_adult_equation <- c(
    "gerris", "halobates", "rhagovelia",
    "saldidae", "pelocoris", "chiloxanthus"
  )
  
  dep.data1.1 <- dep.data1.1 %>%
    mutate(LIFE_STAGE = case_when(
      is.na(TAXANAME1)                         ~ NA_character_, # true NAs
      TAXANAME1 %in% junk_taxa                 ~ NA_character_, # junk
      TAXANAME1 %in% non_insect_taxa           ~ NA_character_, # non-insects
      grepl("\\(l\\)", tolower(Genus))         ~ "L",           # explicit larva
      grepl("\\(a\\)", tolower(Genus))         ~ "A",           # explicit adult
      TAXANAME1 %in% adult_only_taxa           ~ "A",           # adult-only taxa
      TAXANAME1 %in% hemiptera_adult_equation  ~ "A",           # hemiptera adult eq
      TAXANAME1 == "hydrobius" & 
        Order == "Coleoptera"                  ~ "L",           # hydrobius larva
      Order %in% hemimetabolous_orders         ~ "N",           # hemimetabolous
      TRUE                                     ~ "L"            # holometabolous
    ))
  
############################
#Download fixed depositonal habitat datasheets with fixed code 

write.csv(dep.data1.1, "cleaned_typ1_dep_hab_abundance3.30.2026.csv")
  
############################
  
  
  # drop zero-count rows
  size_class_columns <- grep("^X[0-9]+$", names(dep.data1.1), value = TRUE)
  
  dep.data1.1 <- dep.data1.1 %>%
    filter(rowSums(across(all_of(size_class_columns)), na.rm = TRUE) > 0)
  
######################################################################
  # SECTION 5: PREPARE TAXA_INFO LOOKUP TABLE
  # ============================================================
  
  # create TAXANAME1 column (first word, lowercase)
  taxa_info <- taxa_info %>%
    mutate(
      TAXANAME1 = tolower(trimws(word(TAXANAME, 1))),
      TAXANAME  = tolower(trimws(TAXANAME))
    )
  
  # deduplicate — keep first equation per taxon/life stage combo
  # see taxa_info_unmatched_entries.csv for full duplicate list
  taxa_info_dedup <- taxa_info %>%
    group_by(TAXANAME1, LIFE_STAGE) %>%
    slice(1) %>%
    ungroup()
  
  taxa_info_dedup_rif <- taxa_info %>%
    group_by(TAXANAME, LIFE_STAGE) %>%
    slice(1) %>%
    ungroup()
  
  ################################################################
  # SECTION 6: QC — CHECK FOR UNMATCHED TAXA
  ################################################################
  
  dep.data2 <- left_join(dep.data1.1, taxa_info_dedup,
                         by = c("TAXANAME1", "LIFE_STAGE"))
  
  no_match <- dep.data2 %>%
    filter(is.na(Slope)) %>%
    select(TAXANAME1, LIFE_STAGE, Order, Genus) %>%
    distinct() %>%
    arrange(TAXANAME1)
  
  # export QC files
  write.csv(no_match,     "outputs/no_match_taxa_context.csv", row.names = FALSE)
  write.csv(taxa_info,    "outputs/taxa_info_lookup.csv",      row.names = FALSE)
  
  cat("Unmatched taxa:", nrow(no_match), "\n")
  print(no_match)
# Ostracods and daphnia are ommitted and have not length mass regressions
  

  ###############################################################
  # SECTION 7: DEPOSITIONAL HABITAT BIOMASS CALCULATIONS
  # Goal: replicate-level taxon biomass for NMDS
  # No correction factor — raw abundance used directly
  # Unit of replication: Site x Sample x TAXANAME1 x htype
  ###############################################################
  
  dep_biomass_rep <- dep.data1.1 %>%
    pivot_longer(
      cols      = all_of(size_class_columns),
      names_to  = "Size_Class",
      values_to = "Abundance"
    ) %>%
    mutate(Length = as.numeric(sub("X", "", Size_Class))) %>%#make sure that size class columns are numeric
    filter(Abundance > 0) %>%#remove any abundances that are zero 
    group_by(Site, TAXANAME1, LIFE_STAGE, Sample, Length, htype) %>%#group samples by site, name, life stage, sample, length, and habitat 
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%# sumarize the abundance of each combination calcuated above for each length
    left_join(taxa_info_dedup, by = c("TAXANAME1", "LIFE_STAGE")) %>%#join the length mass information to the bug data 
    mutate(
      DM                  = round(Intercept * (Length ^ Slope), 4),#calculate the dry mass of each taxa and round to 4 decimal places 
      Biomass_x_Abundance = round(DM * Abundance, 4),#multiply dry mass of each length by the abundance of taxa in each length class and round 
      Density_per_m2      = round(Abundance * 4, 4),#multiply the densities by 4, since we collected using a .25 m2 sampling area
      Q_area_corrected    = round(Biomass_x_Abundance * 4, 4)#multiply the biomass by 4 since, we collected using a 0.25 m2 sampling area 
    )
  
  dep_biomass_bytaxon <- dep_biomass_rep %>%
    group_by(Site, Sample, TAXANAME1, htype) %>%
    summarise(
      Biomass_m2 = sum(Q_area_corrected, na.rm = TRUE),
      Density_m2 = sum(Density_per_m2,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(htype = as.character(htype))
  
  # verify
  cat("--- Dep hab biomass summary ---\n")
  cat("Unique sites:", n_distinct(dep_biomass_bytaxon$Site), "\n")
  cat("Unique taxa:", n_distinct(dep_biomass_bytaxon$TAXANAME1), "\n")
  cat("Total rows:", nrow(dep_biomass_bytaxon), "\n")
  dep_biomass_bytaxon %>% count(htype)

  ##############################################################
  # SECTION 8: RIFFLE BIOMASS CALCULATIONS
  # Goal: replicate-level taxon biomass for NMDS
  # No correction factor — matches dep hab approach
  # Unit of replication: Site x Sample x TAXANAME
  ##############################################################
  
  size_class_columns_rif <- grep("^X[0-9]+$", names(rif.data2), value = TRUE)
  
####################################
  # standardize riffle TAXANAME to match taxa_info TAXANAME1
  rif.data2 <- rif.data2 %>%
    mutate(TAXANAME1 = tolower(trimws(word(TAXANAME, 1))))
  
  # then join on TAXANAME1 instead of TAXANAME
  rif_biomass_rep <- rif.data2 %>%
    mutate(across(all_of(size_class_columns_rif),
                  ~ replace(., is.na(.), 0))) %>%
    pivot_longer(
      cols      = all_of(size_class_columns_rif),
      names_to  = "Size_Class",
      values_to = "Abundance"
    ) %>%
    mutate(Length = as.numeric(sub("X", "", Size_Class))) %>%
    filter(Abundance > 0) %>%
    group_by(Site, TAXANAME, TAXANAME1, LIFE_STAGE, Sample, Length) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    left_join(taxa_info_dedup, by = c("TAXANAME1", "LIFE_STAGE")) %>%  # join on TAXANAME1
    mutate(
      DM                  = round(Intercept * (Length ^ Slope), 4),
      Biomass_x_Abundance = round(DM * Abundance, 4),
      Density_per_m2      = round(Abundance * 4, 4),
      Q_area_corrected    = round(Biomass_x_Abundance * 4, 4)
    )
  
  # check how many NAs remain after fix
  rif_biomass_rep %>%
    filter(is.na(Slope)) %>%
    count(TAXANAME1, LIFE_STAGE, sort = TRUE)

  rif_biomass_bytaxon <- rif_biomass_rep %>%
    group_by(Site, Sample, TAXANAME1, LIFE_STAGE) %>%
    summarise(
      Biomass_m2 = sum(Q_area_corrected, na.rm = TRUE),
      Density_m2 = sum(Density_per_m2,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      htype     = "riffle"
    )
  
  # verify
  cat("--- Riffle biomass summary ---\n")
  cat("Unique sites:", n_distinct(rif_biomass_bytaxon$Site), "\n")
  cat("Unique taxa:", n_distinct(rif_biomass_bytaxon$TAXANAME1), "\n")
  cat("Total rows:", nrow(rif_biomass_bytaxon), "\n")
  
  ##############################################################
  # SECTION 9: COMBINE AND EXPORT
  #############################################################
  
  all_biomass_bytaxon <- bind_rows(
    rif_biomass_bytaxon %>%
      select(Site, Sample, TAXANAME1, Biomass_m2, Density_m2, htype),
    dep_biomass_bytaxon %>%
      select(Site, Sample, TAXANAME1, Biomass_m2, Density_m2, htype)
  )
  
  # verify
  cat("--- Combined biomass summary ---\n")
  cat("Total rows:", nrow(all_biomass_bytaxon), "\n")
  all_biomass_bytaxon %>% count(htype, sort = TRUE)
  
  # export
  write.csv(all_biomass_bytaxon,
            "outputs/rep_biomass_by_taxon_4.7.26.csv",
            row.names = FALSE)
head(all_biomass_bytaxon)
  
  
  
  #############################
  
  # Check how simulium appears in the RAW riffle data
  rif.data2 %>%
    filter(grepl("simulium|stegopterna", tolower(Genus))) %>%
    distinct(Order, Family, Genus, TAXANAME1) %>%
    arrange(Genus)
  
  # Check how simulium appears in the RAW depositional data
  dep.data1.1 %>%
    filter(grepl("simulium|stegopterna", tolower(TAXANAME1))) %>%
    distinct(Order, Family, Genus, TAXANAME1, htype) %>%
    arrange(TAXANAME1)
  
  # Check what ended up in the final combined biomass file
  all_biomass_bytaxon %>%
    filter(grepl("simulium|stegopterna", tolower(TAXANAME1))) %>%
    count(TAXANAME1, htype, sort = TRUE)
  
  ##############################################################
  # SECTION 9B: VERIFY SITE x SAMPLE ALIGNMENT BETWEEN HABITATS
  ##############################################################
  
  # Get unique Site x Sample combos for each habitat type
  rif_reps <- rif_biomass_bytaxon %>%
    distinct(Site, Sample) %>%
    mutate(in_riffle = TRUE)
  
  dep_reps <- dep_biomass_bytaxon %>%
    distinct(Site, Sample) %>%
    mutate(in_dep = TRUE)
  
  # Full join to find mismatches
  rep_comparison <- full_join(rif_reps, dep_reps, by = c("Site", "Sample")) %>%
    mutate(
      in_riffle = replace_na(in_riffle, FALSE),
      in_dep    = replace_na(in_dep,    FALSE),
      status    = case_when(
        in_riffle & in_dep   ~ "matched",
        in_riffle & !in_dep  ~ "riffle only",
        !in_riffle & in_dep  ~ "dep only"
      )
    ) %>%
    arrange(status, Site, Sample)
  
  # Summary
  cat("--- Site x Sample alignment ---\n")
  rep_comparison %>% count(status) %>% print()
  
  # Print any mismatches for inspection
  mismatches <- rep_comparison %>% filter(status != "matched")
  if (nrow(mismatches) > 0) {
    cat("\nMismatched replicates:\n")
    print(mismatches)
  } else {
    cat("\nAll replicates matched across riffle and depositional habitats.\n")
  }
  
  # Also check sites specifically
  rif_sites <- rif_biomass_bytaxon %>% distinct(Site)
  dep_sites <- dep_biomass_bytaxon %>% distinct(Site)
  
  cat("\n--- Sites in riffle but not dep ---\n")
  print(anti_join(rif_sites, dep_sites, by = "Site"))
  
  cat("\n--- Sites in dep but not riffle ---\n")
  print(anti_join(dep_sites, rif_sites, by = "Site"))
  
  
  
  ######################################
  
  ######################################
  
  # ============================================================
  # PART A: CWM — RIFFLE vs. DEPOSITIONAL, SITE-LEVEL AVERAGED
  # Raw biomass summed per taxon across reps within each site,
  # divided by n reps → mean per-m² value per site.
  # Unit of replication: Site. Stats: Wilcoxon rank-sum.
  # ============================================================
  
  raw_collapsed <- all_biomass_bytaxon %>%
    mutate(
      htype_collapsed = ifelse(htype == "riffle", "riffle", "depositional"),
      htype_collapsed = factor(htype_collapsed,
                               levels = c("depositional", "riffle"))
    )
  
write.csv(raw_collapsed, "outputs/biomass_density_4.28.2026.csv")
  
  n_reps_site_raw <- raw_collapsed %>%
    group_by(Site, htype_collapsed) %>%
    summarise(n_reps = n_distinct(Sample), .groups = "drop")
  
  
  cat("--- Part A: reps per site used for site-level averaging ---\n")
  n_reps_site_raw %>%
    group_by(htype_collapsed, n_reps) %>%
    summarise(n_sites = n(), .groups = "drop") %>%
    print()
  
  raw_collapsed_site <- raw_collapsed %>%
    group_by(Site, htype_collapsed, TAXANAME1) %>%
    summarise(Biomass_m2 = sum(Biomass_m2, na.rm = TRUE), .groups = "drop") %>%
    left_join(n_reps_site_raw, by = c("Site", "htype_collapsed")) %>%
    mutate(Biomass_m2 = Biomass_m2 / n_reps) %>%
    dplyr::select(-n_reps)
  
  write.csv(raw_collapsed_site, "outputs/taxa_biomass_mg.m2_site_4.28.2026.csv")
  
  Site_raw <-raw_collapsed_site %>% 
    group_by(Site, htype_collapsed) %>% 
    summarise(Biomass_m2 = sum(Biomass_m2, na.rm = TRUE))
  
  write.csv(Site_raw, "outputs/Biomass_mg.m2_site_4.28.2026.csv")
  
  
  # --- Carry both Biomass and Density through site-level averaging ---
  raw_collapsed_site <- raw_collapsed %>%
    group_by(Site, htype_collapsed, TAXANAME1) %>%
    summarise(
      Biomass_m2  = sum(Biomass_m2,  na.rm = TRUE),
      Density_m2  = sum(Density_m2,  na.rm = TRUE),   # add density here
      .groups = "drop"
    ) %>%
    left_join(n_reps_site_raw, by = c("Site", "htype_collapsed")) %>%
    mutate(
      Biomass_m2 = Biomass_m2 / n_reps,
      Density_m2 = Density_m2 / n_reps
    ) %>%
    dplyr::select(-n_reps)
  
  write.csv(raw_collapsed_site, "outputs/taxa_biomass_mg.m2_site_4.28.2026.csv", row.names = FALSE)
  
  
  # --- Ranked taxa per site x habitat type ---
  taxa_ranked <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    arrange(desc(Biomass_m2), .by_group = TRUE) %>%          # most abundant first
    mutate(
      rank_biomass = row_number(),
      pct_biomass  = Biomass_m2 / sum(Biomass_m2, na.rm = TRUE) * 100
    ) %>%
    ungroup() %>%
    arrange(Site, htype_collapsed, rank_biomass)
  
  write.csv(taxa_ranked, "outputs/taxa_ranked_biomass_density_site_4.28.2026.csv", row.names = FALSE)
  
  # --- Ranked taxa per site x habitat type ---
  taxa_ranked_density <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    arrange(desc(Density_m2), .by_group = TRUE) %>%          # most abundant first
    mutate(
      rank_density = row_number(),
      pct_density  = Density_m2 / sum(Density_m2, na.rm = TRUE) * 100
    ) %>%
    ungroup() %>%
    arrange(Site, htype_collapsed, rank_density)
  
  write.csv(taxa_ranked_density, "outputs/taxa_ranked_density_site_4.28.2026.csv", row.names = FALSE)
  
  
  # --- Ranked taxa per site x habitat type by 5  ---
  taxa_ranked_by5 <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    arrange(desc(Biomass_m2), .by_group = TRUE) %>%          # most abundant first
    mutate(
      rank_biomass = row_number(),
      pct_biomass  = Biomass_m2 / sum(Biomass_m2, na.rm = TRUE) * 100
    ) %>%
    filter(rank_biomass <= 5) %>% 
    ungroup() %>%
    arrange(Site, htype_collapsed, rank_biomass)
  
  write.csv(taxa_ranked_by5, "outputs/taxa_ranked_biomass_density_top5_site_4.28.2026.csv", row.names = FALSE)
  
  taxa_ranked_density
  
  
  ##############################################################
  # TOP 3 TAXA BY BIOMASS AND DENSITY PER SITE x HABITAT
  # Input: raw_collapsed_site (already rep-averaged mg/m2 and no./m2)
  ##############################################################
  
  # --- Top 3 by biomass ---
  top3_biomass <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    mutate(
      pct_biomass  = Biomass_m2 / sum(Biomass_m2, na.rm = TRUE) * 100,
      rank_biomass = dense_rank(desc(Biomass_m2))
    ) %>%
    filter(rank_biomass <= 3) %>%
    ungroup() %>%
    dplyr::select(Site, htype_collapsed, rank_biomass, TAXANAME1,
                  Biomass_m2, pct_biomass) %>%
    arrange(Site, htype_collapsed, rank_biomass)
  
  # --- Top 3 by density ---
  top3_density <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    mutate(
      pct_density  = Density_m2 / sum(Density_m2, na.rm = TRUE) * 100,
      rank_density = dense_rank(desc(Density_m2))
    ) %>%
    filter(rank_density <= 3) %>%
    ungroup() %>%
    dplyr::select(Site, htype_collapsed, rank_density, TAXANAME1,
                  Density_m2, pct_density) %>%
    arrange(Site, htype_collapsed, rank_density)
  
  # --- Side-by-side table: rank 1-3 for each metric ---
  top3_combined <- full_join(
    top3_biomass  %>% rename(rank = rank_biomass,
                             taxon_biomass = TAXANAME1),
    top3_density  %>% rename(rank = rank_density,
                             taxon_density = TAXANAME1),
    by = c("Site", "htype_collapsed", "rank")
  ) %>%
    dplyr::select(Site, htype_collapsed, rank,
                  taxon_biomass, Biomass_m2, pct_biomass,
                  taxon_density, Density_m2, pct_density) %>%
    arrange(Site, htype_collapsed, rank)
  
  write.csv(top3_biomass,  "outputs/top3_taxa_biomass_site_htype.csv", row.names = FALSE)
  write.csv(top3_density,  "outputs/top3_taxa_density_site_htype.csv", row.names = FALSE)
  write.csv(top3_combined, "outputs/top3_taxa_biomass_density_site_htype.csv", row.names = FALSE)
  
  # quick look
  top3_combined %>% print(n = 30)
  
  
  
  # ============================================================
  # TOP 3 TAXA TABLES WITH SCIENTIFIC NAMES
  # Produces two long-format tables (one row per taxon):
  #   top3_biomass_tbl  - top 3 taxa by mean biomass (mg/m2)
  #   top3_density_tbl  - top 3 taxa by mean density (no./m2)
  # Taxon names rendered "Family (Genus)"
  # Requires: raw_collapsed_site, dep.data1.1, rif.data2
  # ============================================================
  
  # ---- 1. Build TAXANAME1 -> Order / Family / Genus lookup ----
  taxon_lookup_raw <- bind_rows(
    dep.data1.1 %>% dplyr::select(TAXANAME1, Order, Family, Genus),
    rif.data2   %>% dplyr::select(TAXANAME1, Order, Family, Genus)
  ) %>%
    mutate(
      # strip parenthetical life-stage tags e.g. "Hydraena(A)"
      across(c(Order, Family, Genus),
             ~ trimws(gsub("\\s*\\([^)]*\\)", "", as.character(.x)))),
      across(c(Order, Family, Genus), ~ na_if(.x, ""))
    ) %>%
    filter(!is.na(TAXANAME1))
  
  # most frequent value, used to resolve conflicting metadata
  mode_chr <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }
  
  taxon_lookup <- taxon_lookup_raw %>%
    group_by(TAXANAME1) %>%
    summarise(
      Order     = mode_chr(Order),
      Family    = mode_chr(Family),
      Genus     = mode_chr(Genus),
      n_genera  = n_distinct(Genus[!is.na(Genus)]),  # >1 = pooled taxon, check it
      .groups   = "drop"
    )
  
  # ---- 2. Build display names ----
  proper_case <- function(x) {
    out <- vapply(x, function(s) {
      if (is.na(s)) return(NA_character_)
      s     <- gsub("_", "/", s)
      parts <- strsplit(s, "/")[[1]]
      parts <- paste0(toupper(substr(parts, 1, 1)), tolower(substr(parts, 2, nchar(parts))))
      paste(parts, collapse = "/")
    }, character(1), USE.NAMES = FALSE)
    out
  }
  
  norm_name <- function(x) gsub("[ _]+", "/", tolower(trimws(as.character(x))))
  
  taxon_names <- taxon_lookup %>%
    mutate(
      # is TAXANAME1 itself the genus?
      is_genus   = !is.na(Genus) & norm_name(Genus) == norm_name(TAXANAME1),
      Rank_level = case_when(
        is_genus                        ~ "Genus",
        grepl("idae$|inae$", TAXANAME1) ~ "Family",
        TRUE                            ~ "Higher"
      ),
      Family_disp = proper_case(ifelse(Rank_level == "Family", TAXANAME1, Family)),
      Genus_disp  = ifelse(Rank_level == "Genus", proper_case(TAXANAME1), NA_character_),
      Order_disp  = proper_case(Order),
      Taxon = case_when(
        Rank_level == "Genus"  ~ Genus_disp,          # lowest = genus
        Rank_level == "Family" ~ Family_disp,         # lowest = family
        TRUE                   ~ proper_case(TAXANAME1)  # order/class/higher
      )
    ) %>%
    dplyr::select(TAXANAME1, Taxon, Order = Order_disp, Family = Family_disp,
                  Genus = Genus_disp, Rank_level, n_genera)
  
  # QC: taxa with no family resolved, or family-level names pooling >1 genus
  cat("--- Taxon name QC ---\n")
  taxon_names %>%
    filter(is.na(Family) | (Rank_level != "Genus" & n_genera > 1)) %>%
    arrange(Rank_level, TAXANAME1) %>%
    print(n = Inf)
  
  # ---- 3. TABLE 1: top 3 by biomass ----
  top3_biomass_tbl <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    mutate(
      Pct  = Biomass_m2 / sum(Biomass_m2, na.rm = TRUE) * 100,
      Rank = dense_rank(desc(Biomass_m2))
    ) %>%
    filter(Rank <= 3) %>%
    ungroup() %>%
    left_join(taxon_names, by = "TAXANAME1") %>%
    transmute(
      Site,
      Habitat       = stringr::str_to_title(as.character(htype_collapsed)),
      Rank,
      Taxon,
      Order, Family, Genus, Rank_level,
      Biomass_mg_m2 = round(Biomass_m2, 2),
      Pct_biomass   = round(Pct, 1)
    ) %>%
    arrange(Site, Habitat, Rank)
  
  # ---- 4. TABLE 2: top 3 by density ----
  top3_density_tbl <- raw_collapsed_site %>%
    group_by(Site, htype_collapsed) %>%
    mutate(
      Pct  = Density_m2 / sum(Density_m2, na.rm = TRUE) * 100,
      Rank = dense_rank(desc(Density_m2))
    ) %>%
    filter(Rank <= 3) %>%
    ungroup() %>%
    left_join(taxon_names, by = "TAXANAME1") %>%
    transmute(
      Site,
      Habitat        = stringr::str_to_title(as.character(htype_collapsed)),
      Rank,
      Taxon,
      Order, Family, Genus, Rank_level,
      Density_no_m2  = round(Density_m2, 1),
      Pct_density    = round(Pct, 1)
    ) %>%
    arrange(Site, Habitat, Rank)
  
  # ---- 5. Export ----
  write.csv(top3_biomass_tbl,
            "outputs/Table_top3_taxa_biomass_site.csv", row.names = FALSE)
  write.csv(top3_density_tbl,
            "outputs/Table_top3_taxa_density_site.csv", row.names = FALSE)
  
  aacat("\n--- Top 3 by biomass ---\n"); print(top3_biomass_tbl, n = 30)
  cat("\n--- Top 3 by density ---\n"); print(top3_density_tbl, n = 30)
  
  # Any taxon that failed to get a name
  cat("\nRows with unresolved Taxon:",
      sum(is.na(top3_biomass_tbl$Taxon)) + sum(is.na(top3_density_tbl$Taxon)), "\n")
  # install.packages(c("flextable", "officer"))
  library(flextable)
  library(officer)
  
  set_flextable_defaults(
    font.family   = "Times New Roman",
    font.size     = 9,
    padding       = 3,
    border.color  = "grey30",
    table.layout  = "autofit"
  )
  
  # ---- Builder ------------------------------------------------
  make_top3_ft <- function(dat, value_col, pct_col,
                           value_main, value_unit, caption_txt) {
    d <- dat %>%
      mutate(Value = .data[[value_col]],
             Pct   = .data[[pct_col]]) %>%
      arrange(Site, Habitat, Rank)
    
    flextable(d, col_keys = c("Site","Habitat","Rank","Taxon","Value","Pct")) %>%
      compose(part = "header", j = "Value",
              value = as_paragraph(value_main, " (", value_unit, as_sup("−2"), ")")) %>%
      set_header_labels(Site = "Site", Habitat = "Habitat", Rank = "Rank",
                        Taxon = "Taxon", Pct = "% of total") %>%
      italic(i = ~ Rank_level == "Genus", j = "Taxon") %>%   # genera only
      merge_v(j = c("Site", "Habitat")) %>%
      valign(j = c("Site", "Habitat"), valign = "top", part = "body") %>%
      align(j = c("Rank","Value","Pct"), align = "right", part = "all") %>%
      align(j = c("Site","Habitat","Taxon"), align = "left", part = "all") %>%
      colformat_double(j = "Value", digits = 2, big.mark = ",") %>%
      colformat_double(j = "Pct",   digits = 1) %>%
      theme_booktabs() %>%
      hline(i = ~ !is.na(Site), border = fp_border(color = "grey80", width = 0.5)) %>%
      bold(part = "header") %>%
      set_caption(caption_txt) %>%
      autofit()
  }
  # ---- Table 1: biomass ---------------------------------------
  ft_biomass <- make_top3_ft(
    top3_biomass_tbl,
    value_col   = "Biomass_mg_m2",
    pct_col     = "Pct_biomass",
    value_main  = "Mean biomass",
    value_unit  = "mg m",
    caption_txt = paste("Table 1. Three most dominant macroinvertebrate taxa by",
                        "mean dry biomass in depositional and riffle habitats at",
                        "each study site. Values are replicate means; percentages",
                        "are of total assemblage biomass within a site and habitat.")
  )
  
  # ---- Table 2: density ---------------------------------------
  ft_density <- make_top3_ft(
    top3_density_tbl,
    value_col   = "Density_no_m2",
    pct_col     = "Pct_density",
    value_main  = "Mean density",
    value_unit  = "ind. m",
    caption_txt = paste("Table 2. Three most dominant macroinvertebrate taxa by",
                        "mean density in depositional and riffle habitats at each",
                        "study site. Values are replicate means; percentages are",
                        "of total assemblage density within a site and habitat.")
  )
  
  # ---- View ---------------------------------------------------
  ft_biomass
  ft_density
  
  # ---- Export to the figure repository ------------------------
  dir.create("outputs/", showWarnings = FALSE, recursive = TRUE)
  
  # Word (best for pasting into a thesis chapter, keeps italics)
  save_as_docx(
    `Table 1. Top 3 taxa by biomass` = ft_biomass,
    `Table 2. Top 3 taxa by density` = ft_density,
    path = "outputs/Top3_taxa_tables.docx",
    pr_section = prop_section(
      page_size = page_size(orient = "portrait"),
      type = "continuous"
    )
  )
  
  
  # PNG (needs webshot2 + chromote installed)
  # save_as_image(ft_biomass, path = "outputs/figures/Table1_top3_biomass.png")
  # save_as_image(ft_density, path = "outputs/figures/Table2_top3_density.png")
  
  # PowerPoint, if your figure repo takes editable vector objects
  # save_as_pptx(`Top3 biomass` = ft_biomass,
  #              `Top3 density` = ft_density,
  #              path = "outputs/figures/Top3_taxa_tables.pptx")
  
  