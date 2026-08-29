#originator: Brice Crum 
#Project: USGS BMP 
#Thesis chapter: Chapter 1 
#Updated: added ws_forest model (m4) and AIC-based model comparison table

#Hypothesis: H1: Depositional habitat diversity will be greatest in 
#stream reaches with the most percent forest land cover in the riparian
#zone because of more organic matter inputs and stream bank stability 
#that promotes reach-level retention of habitats 

#set working directory to Aquatic Entomology Lab Drive 
setwd( "G:/Shared drives/Aquatic Entomology Lab/Projects/USGS BMP/Year 1_SHEN_AG
       /Brice's Corner/BC thesis materials/BC Thesis")

#read.DPLYR 
library(dplyr)
library(stringr)
library(ggplot2)
library(vegan)

#read in raw depositional habitat data 
hab <- read.csv("Inputs/DEP_HAB.csv")

#create a column with the lumped vegetation data 
hab$veg <- hab$EV + hab$SV + hab$MB + hab$Ag

#remove the vegetation categories that we are not using for the thesis 
hab <- subset(hab, select = -c(EV, SV, MB, Ag)) 

#remove Snag 

hab <- hab%>% 
  mutate(DA_LW = DA+LW) %>% 
  select(-DA, -LW) 


hab


hab_vars <- hab[, c("UB","DA_LW","RT","veg")]

hab_vars <- hab_vars %>% 
  rename(
    "Undercut bank"      = UB,
    "Debris accumulation" = DA_LW,
    "Emergent roots"      = RT,
    "Instream vegetation" = veg
  )
hab

#------------------------------------#
#generate a df with the means and sd 

hab_summary <- data.frame(
  variable = names(hab)[-1], 
  mean = sapply(hab[,-1], mean, na.rm = TRUE),
  sd   = sapply(hab[,-1], sd,   na.rm = TRUE)
)

#-----------------------------------#

#remove the site column for generation of Shannon's diversity 
hdf <- hab[, -c(1)]

#generate Shannon's diversity 
H <- diversity(hdf, "shannon")
print(H)

#combine H with data frame
hab$H <- H


#add the number of habitats to the plot
hab$number_hab <- rowSums(hab_vars >0, na.rm=TRUE)
hab$number_hab



#--------------------------land use variable management-----------------------#

ws_forest <- c(
  "forest","other_tree_canopy",
  "riverine_wetlands_tree_canopy","riverine_wetlands_forest",
  "terrene_wetlands_tree_canopy","terrene_wetlands_forest",
  "tidal_wetlands_tree_canopy","tidal_wetlands_forest"
)
ws_agri <- c(
  "cropland_barren","cropland_herbaceous",
  "pasture_hay_barren","pasture_hay_herbaceous","pasture_hay_scrub_shrub",
  "orchard_vineyard_herbaceous","orchard_vineyard_scrub_shrub",
  "harvested_forest_barren","harvested_forest_herbaceous"
)
ws_urban <- c(
  "roads","structures","other_impervious",
  "tree_canopy_over_roads","tree_canopy_over_structures",
  "tree_canopy_over_other_impervious","tree_canopy_over_turf_grass",
  "turf_grass","transitional_barren","extractive_barren","solar_field_herbaceous"
)
ws_water  <- c("lakes_and_reservoirs","riverine_ponds","terrene_ponds","lotic_water_fresh","bare_shore")
ws_wet_nf <- c(
  "riverine_wetlands_barren","riverine_wetlands_herbaceous","riverine_wetlands_scrub_shrub",
  "terrene_wetlands_barren","terrene_wetlands_herbaceous","terrene_wetlands_scrub_shrub",
  "tidal_wetlands_herbaceous","tidal_wetlands_scrub_shrub"
)
ws_succ <- c(
  "suspended_succession_barren","suspended_succession_herbaceous","suspended_succession_scrub_shrub",
  "natural_succession_barren","natural_succession_herbaceous","natural_succession_scrub_shrub"
)
rip_imperv <- c("rip_impervious_roads","rip_impervious_structures","rip_impervious_other",
                "rip_tree_canopy_over_impervious")
rip_turf   <- c("rip_tree_canopy_over_turf_grass","rip_turf_grass","rip_pervious_developed_other")
rip_wet_nf <- c("rip_wetlands_riverine_non_forested","rip_wetlands_terrene_non_forested",
                "rip_wetlands_tidal_non_forested")
rip_agri   <- c("rip_cropland","rip_pasture_hay","rip_harvested_forest")

#-----------------------------------------------------------------------------#

#open watershed land use from the Microsoft Teams channel 2/23/2026
LU <- read.csv("Inputs/WS_4_landuse.csv") # data dated 4/16/2025 on the drive

#subtract riparian water cover from forest cover — values overlap so subtract
LU$RIP_Forest_trans <- LU$RIP_FOREST - LU$RIP_WATER

#create watershed agriculture variable by summing all agricultural subcategories
LU$ws_ag <- rowSums(LU[, toupper(ws_agri)], na.rm = TRUE)

#create riparian agriculture variable
LU$rip_ag <- rowSums(LU[, toupper(rip_agri)], na.rm = TRUE)

#create watershed forest variable
LU$ws_forest <- rowSums(LU[, toupper(ws_forest)], na.rm = TRUE)

LU$rip_urban <- rowSums(LU[, toupper(rip_imperv)], na.rm= TRUE)

LU$ws_urban <- rowSums(LU[, toupper(ws_urban)], na.rm=TRUE)

#remove fish watersheds from the LU object 
LU <- LU %>% 
  filter(!grepl("_FISH$", SITE_ID))

#remove "_MAIN" from the LU site code to allow join with dep hab data
LU <- LU %>% 
  mutate(SITE_ID = str_remove(SITE_ID, "_MAIN"))

#create duplicate column of SITE_ID as SITE for joining
LU$SITE <- LU$SITE_ID

#inner join depositional habitat data and land use data by SITE
df <- LU %>% 
  inner_join(hab, LU, by = "SITE")

#write CSV for future use 
write.csv(df, file = "dephab_usgsLU_2.25.2026.csv")

threshold <- 40

high <- dplyr::filter(df, RIP_Forest_trans >= threshold)
low  <- dplyr::filter(df, RIP_Forest_trans<  threshold)

highrng <- range(high$UB, high$RT, high$DA_LW, high$veg)
highrng

lowrng <- range(low$UB, low$RT, low$DA_LW, low$veg)
lowrng

#-----------exploratory histograms-----------#

hist(df$RIP_Forest_trans)  # riparian forest — bimodal, peak at 10-30% and 60-70%
hist(df$ws_ag)             # watershed agriculture
hist(df$rip_ag)            # riparian agriculture
hist(df$ws_forest)         # watershed forest cover
hist(df$H)                 # Shannon diversity — approximately normal


#----------exploratory regressions---------#

#agricultural land
ag_lm <- lm(LU$rip_ag~LU$ws_ag)
 
 summary(ag_lm)
 
#forest land 
fs_lm <- lm(LU$RIP_Forest_trans~LU$ws_forest)

summary(fs_lm)

df1 <- subset(LU, select = c(RIP_Forest_trans, ws_forest, ws_ag, rip_ag, ws_urban, rip_urban))
              
cor(df1)
#============================================================#
#        MODEL BUILDING — all single-predictor LMs           #
#============================================================#

# m1: riparian forest cover (corrected for water overlap)
m1 <- lm(H ~ RIP_Forest_trans, df)
summary(m1) # Adj R2 = 0.04, not significant

# m2: watershed-scale agriculture
m2 <- lm(H ~ ws_ag, df)
summary(m2) # Adj R2 = 0.045, not significant

# m3: riparian-scale agriculture
m3 <- lm(H ~ rip_ag, df)
summary(m3) # Adj R2 = 0.045, not significant

# m4: watershed-scale forest cover  *** NEW ***
m4 <- lm(H ~ ws_forest, df)
summary(m4)

#-----------check m4 assumptions-----------#

fit_vals_m4  <- m4$fitted.values
stud_resid_m4 <- rstudent(m4)

# 1. Constant variance
plot(fit_vals_m4, stud_resid_m4, pch = 16,
     xlab = "Fitted values", ylab = "Studentized residuals",
     main = "m4: Residuals vs Fitted (ws_forest)")
abline(h = 0, col = "blue")

# 2. Normality of residuals
qqnorm(stud_resid_m4, main = "m4: Normal Q-Q (ws_forest)")
abline(0, 1, col = "red")
shapiro.test(stud_resid_m4)

# 3. Linearity — scatter with regression line
plot(df$ws_forest, df$H,
     xlab = "% Watershed forest cover",
     ylab = "Shannon depositional habitat diversity",
     main = "m4: ws_forest ~ H")
abline(m4, col = "red")

#============================================================#
#           AIC-BASED MODEL COMPARISON TABLE                 #
#============================================================#

model_comparison <- data.frame(
  Model     = c("m1", "m2", "m3", "m4"),
  Predictor = c("Riparian forest (corrected)", 
                "Watershed agriculture", 
                "Riparian agriculture",
                "Watershed forest"),
  Scale     = c("Riparian", "Watershed", "Riparian", "Watershed"),
  AIC       = c(AIC(m1), AIC(m2), AIC(m3), AIC(m4)),
  AdjR2     = c(summary(m1)$adj.r.squared,
                summary(m2)$adj.r.squared,
                summary(m3)$adj.r.squared,
                summary(m4)$adj.r.squared),
  P_value   = c(summary(m1)$coefficients[2, 4],
                summary(m2)$coefficients[2, 4],
                summary(m3)$coefficients[2, 4],
                summary(m4)$coefficients[2, 4])
)

model_comparison <- model_comparison[order(model_comparison$AIC), ]
print(model_comparison)

# Write model comparison table to CSV
write.csv(model_comparison, file = "outputs/H1_model_comparison.csv", row.names = FALSE)

#============================================================#
#      GGPLOT FOR BEST MODEL (update "best_model" below)     #
#============================================================#

# Inspect model_comparison and set best_model to whichever has lowest AIC
# For now m4 is included — swap in whichever wins

best_model  <- m1
best_pred   <- df$RIP_Forest_trans
best_xlabel <- "Riparian Forest Cover (%)"

bm_r2      <- round(summary(best_model)$r.squared, 3)
bm_p       <- round(summary(best_model)$coefficients[2, 4], 3)
bm_int     <- round(coef(best_model)[1], 2)
bm_slope   <- round(coef(best_model)[2], 3)
bm_p_label <- ifelse(bm_p < 0.001, "p < 0.001", paste0("p = ", bm_p))

bm_label <- paste0(
  "y = ", bm_int, " + ", bm_slope, "x\n",
  "R² = ", bm_r2, ",  ", bm_p_label
)

best_plot <- ggplot(data = df, aes(x = RIP_Forest_trans, y = H, color = number_hab)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    x       = best_xlabel,
    y       = "Habitat diversity (H')",
    color   = "Number of habitats",
    caption = bm_label
  ) +
  theme(
    legend.position  = "bottom",
    axis.title.x     = element_text(size = 25),
    axis.title.y     = element_text(size = 25),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    legend.text      = element_text(size = 15),
    legend.title     = element_text(size = 20),
    plot.caption     = element_text(size = 20, hjust = 0)
  )

best_plot


ggsave("outputs/h1_best_model_plot.png", plot = best_plot,
       width = 16, height = 8, dpi = 600, bg = "transparent")

#============================================================#
#              HABITAT ACCUMULATION CURVE                    #
#============================================================#

hab_accum <- specaccum(
  comm         = hdf,
  method       = "random",
  permutations = 999
)

accum_df <- data.frame(
  Sites    = hab_accum$sites,
  Richness = hab_accum$richness,
  SD       = hab_accum$sd
)

hab_accum_plot <- ggplot(accum_df, aes(x = Sites, y = Richness)) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD),
              fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_hline(yintercept = ncol(hdf),
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  annotate("text",
           x     = -Inf, y = -Inf,
           label = paste0("Max habitats = ", ncol(hdf)),
           hjust = -0.1, vjust = -0.5,
           size  = 4.5, color = "grey40") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                     limits = c(0, ncol(hdf) + 0.5)) +
  labs(
    title    = "Habitat Accumulation Curve",
    subtitle = "Randomized site accumulation (999 permutations) | Ribbon = ±1 SD",
    x        = "Number of sites",
    y        = "Habitats detected"
  ) +
  theme_classic() +
  theme(
    axis.title.x  = element_text(size = 20),
    axis.title.y  = element_text(size = 20),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 12),
    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40")
  )

hab_accum_plot

ggsave("outputs/h1_habitat_accumulation.png", plot = hab_accum_plot,
       width = 16, height = 8, dpi = 600, bg = "transparent")



library(flextable)
library(officer)
library(dplyr)

# 1. DATA
df <- data.frame(
  Model     = c("m1", "m2", "m3", "m4"),
  Predictor = c("Riparian Forest (%)",
                "Watershed Forest (%)",
                "Watershed Agriculture (%)",
                "Riparian Agriculture (%)"),
  Int_Est   = c(0.503, 0.477, 0.935, 0.878),
  Int_SE    = c(0.111, 0.144, 0.186, 0.154),
  B_Est     = c(0.005,  0.005, -0.005, -0.004),
  B_SE      = c(0.003,  0.003,  0.004,  0.003),
  t         = c(1.756,  1.463, -1.501, -1.468),
  p         = c(0.090,  0.155,  0.145,  0.153),
  R2        = c(0.099,  0.071,  0.074,  0.071),
  Adj_R2    = c(0.067,  0.038,  0.041,  0.038),
  AIC       = c(19.87,  20.80,  20.69,  20.79),
  stringsAsFactors = FALSE
)

# 2. RENAME COLUMNS
new_names <- c(
  "Model",
  "Predictor",
  "Intercept",
  "SE (Intercept)",
  "\u03b2",
  "SE (\u03b2)",
  "t",
  "p",
  "R\u00b2",
  "Adj. R\u00b2",
  "AIC"
)

df_display <- setNames(df, new_names)

# 3. BUILD FLEXTABLE
ft <- flextable(df_display) %>%
  
  # colwidths: 1 + 1 + 2 + 3 + 4 = 11
  add_header_row(
    values    = c("", "", "Intercept", "Predictor Coefficient", "Model Fit"),
    colwidths = c(1,  1,   2,           3,                        4)
  ) %>%
  
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  
  align(part = "header", align = "center") %>%
  align(part = "body",   align = "center") %>%
  align(j = 1:2, part = "body", align = "left") %>%
  
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = 1, part = "header",
        border = fp_border(color = "black", width = 0.75)) %>%
  hline_bottom(part = "header",
               border = fp_border(color = "black", width = 0.75)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.5)) %>%
  
  bg(i = seq(2, nrow(df_display), by = 2),
     bg = "#F5F5F5", part = "body") %>%
  
  width(j = 1,     width = 0.50) %>%
  width(j = 2,     width = 1.80) %>%
  width(j = 3:4,   width = 0.75) %>%
  width(j = 5:7,   width = 0.65) %>%
  width(j = 8,     width = 0.60) %>%
  width(j = 9:10,  width = 0.65) %>%
  width(j = 11,    width = 0.60)

# 4. EXPORT
doc <- read_docx() %>%
  body_add_par(
    paste0("Table S1. Simple linear regression models predicting depositional habitat ",
           "diversity (H\u2032) from watershed and riparian land cover predictors ",
           "(n = 30 sites, df = 28 for all models). \u03b2 = unstandardized regression ",
           "coefficient. Riparian Forest was arcsine square-root transformed prior to analysis."),
    style = "Normal"
  ) %>%
  body_add_flextable(ft)

print(doc, target = "h1model_summary_table.docx")
cat("Done!\n")

print(ft)

