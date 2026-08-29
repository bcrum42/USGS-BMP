# =============================================================================
# Correlation plot: biomass NMDS axes vs environmental variables
# Files are read from the working directory.
# =============================================================================

library(corrplot)

# ---- read data --------------------------------------------------------------

nmds <- read.csv("biomass_nmds_scores.csv")
hab  <- read.csv("G:/Shared drives/Aquatic Entomology Lab/Projects/USGS BMP/USGS Teams Data Backup/5.29.2026 unpublished and updated data/Habitat/Typology 1/Typ1_Summary_Metrics.csv")
lu   <- read.csv("dephab_usgsLU_2.25.2026.csv")
wq   <- read.csv("Shenandoah_2021_WQ_with_site_IDs.csv")

hab$Stream_Slope


# ---- water quality: average the two 2021 visits per site --------------------
# The WQ columns contain USGS qualifier codes ("nd", "x", "--"), so they read in
# as text. as.numeric turns those codes into NA and leaves real numbers alone.

wq_site <- data.frame(Site = wq$Short_local_site_ID)
wq_site$SpCond     <- as.numeric(wq$Specific_conductance_us_per_cm_at_25C)
wq_site$DO         <- as.numeric(wq$Dissolved_oxygen_mg_per_l)
wq_site$pH         <- as.numeric(wq$pH_standard_units)
wq_site$Turbidity  <- as.numeric(wq$Turbidity_FNU)
wq_site$Alkalinity <- as.numeric(wq$Alkalinity_mg_per_l_of_CaCO3_wf)
wq_site$TotalN     <- as.numeric(wq$Total_N_mg_per_l_wu)
wq_site$TotalP     <- as.numeric(wq$Phosphorus_mg_per_l_of_P_wu)
#wq_site$DOC        <- as.numeric(wq$DOC_mg_per_l)
#wq_site$Chloride   <- as.numeric(wq$Chloride_mg_per_l_wf_NWQL)
#wq_site$Sulfate    <- as.numeric(wq$Sulfate_mg_per_l_wf)
#wq_site$Nitrate    <- as.numeric(wq$Nitrate_mg_per_l_of_N_wf_NWIS)

wq_site <- aggregate(. ~ Site, data = wq_site, FUN = mean, na.rm = TRUE,
                     na.action = na.pass)

range(wq_site$)
# ---- habitat and land cover (already one row per site) ----------------------

hab_site <- data.frame(Site         = hab$Internal_Site_ID,
                       Embeddedness = hab$Embed_Mean,
                       D50          = hab$Peb_D50,
                       stream_slope = hab$Stream_Slope)

lu_site <- data.frame(Site       = lu$SITE_ID,
                      ShannonHab = lu$H,
                      ws_ag      = lu$ws_ag,
                      ws_forest  = lu$ws_forest,
                      ws_urban   = lu$ws_urban,
                      rip_ag     = lu$rip_ag,
                      rip_forest = lu$RIP_Forest_trans,
                      rip_urban  = lu$rip_urban)


# ---- combine ----------------------------------------------------------------

env <- merge(wq_site, hab_site, by = "Site")
env <- merge(env, lu_site, by = "Site")

# Averaging the two visits leaves floating-point dust: a pH of 7.7 can come out
# as 7.69999999999999929, which Spearman then ranks as untied. Round it away.
env[-1] <- round(env[-1], 10)

dat <- merge(nmds[c("Site", "NMDS1", "NMDS2")], env, by = "Site")

vars <- c("NMDS1", "NMDS2",
          "SpCond", "DO", "pH", "Turbidity", "Alkalinity", "TotalN", "TotalP",
          "Embeddedness", "D50", "ShannonHab", "stream_slope",
          "ws_ag", "ws_forest", "ws_urban",
          "rip_ag", "rip_forest", "rip_urban")


# ---- correlations -----------------------------------------------------------

M <- cor(dat[vars], method = "spearman", use = "pairwise.complete.obs")
P <- cor.mtest(dat[vars], method = "spearman", exact = FALSE)$p
dimnames(P) <- dimnames(M)
p
M

p <- cor.mtest(dat[vars], conf.level = 0.95)$p   #
# ---- plots ------------------------------------------------------------------

# everything vs everything
png("corrplot_full.png", width = 10, height = 10, units = "in", res = 300)
plot <- corrplot(M, p.mat = p, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.85,
         addCoef.col = "black", number.cex = 0.62, number.digits = 2,
         sig.level = 0.05, insig = "blank",
         mar = c(0, 0, 2, 0), title = "...")
plot
# just the two NMDS axes against the predictors
preds <- setdiff(vars, c("NMDS1", "NMDS2"))
png("corrplot_nmds.png", width = 11, height = 3.5, units = "in", res = 300)
# insig = "blank" leaves non-significant cells empty; that keeps the printed
# coefficients from colliding with significance stars in the same cell.
corrplot(M[c("NMDS1", "NMDS2"), preds], p.mat = P[c("NMDS1", "NMDS2"), preds],
         method = "circle", tl.col = "black", tl.srt = 90, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.65, cl.pos = "r",
         sig.level = 0.05, insig = "blank",
         mar = c(0, 0, 2, 0),
         title = "Biomass NMDS axes vs environment (Spearman, p < 0.05)")
dev.off()



# ---- save the numbers -------------------------------------------------------

write.csv(round(M, 3), "correlation_r.csv")
write.csv(signif(P, 3), "correlation_p.csv")
