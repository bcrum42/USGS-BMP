# ============================================================
# TOP 3 TAXA - FLEXTABLE OUTPUT WITH ITALICIZED GENERA
# Requires: top3_biomass_tbl, top3_density_tbl
# Genus names inside "Family (Genus)" are italicized; family
# and higher-level names stay roman.
# ============================================================

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
    mutate(
      Value = .data[[value_col]],
      Pct   = .data[[pct_col]],
      # split the display name so only the genus gets italics
      fam_txt = case_when(
        !is.na(Family)                ~ Family,
        is.na(Family) & !is.na(Genus) ~ "",     # genus only -> all italic
        TRUE                          ~ Taxon   # higher-level -> roman
      ),
      gen_txt = ifelse(is.na(Genus), "", Genus),
      open    = ifelse(!is.na(Genus) & !is.na(Family), " (", ""),
      close   = ifelse(!is.na(Genus) & !is.na(Family), ")",  "")
    ) %>%
    arrange(Site, Habitat, Rank)

  ft <- flextable(
    d,
    col_keys = c("Site", "Habitat", "Rank", "Taxon", "Value", "Pct")
  ) %>%
    # italic genus inside the taxon cell
    compose(
      j = "Taxon",
      value = as_paragraph(as_chunk(fam_txt), as_chunk(open),
                           as_i(gen_txt),     as_chunk(close))
    ) %>%
    # headers
    compose(part = "header", j = "Value",
            value = as_paragraph(value_main, " (", value_unit,
                                 as_sup("−2"), ")")) %>%
    set_header_labels(
      Site    = "Site",
      Habitat = "Habitat",
      Rank    = "Rank",
      Taxon   = "Taxon",
      Pct     = "% of total"
    ) %>%
    merge_v(j = c("Site", "Habitat")) %>%
    valign(j = c("Site", "Habitat"), valign = "top", part = "body") %>%
    align(j = c("Rank", "Value", "Pct"), align = "right", part = "all") %>%
    align(j = c("Site", "Habitat", "Taxon"), align = "left", part = "all") %>%
    colformat_double(j = "Value", digits = 2, big.mark = ",") %>%
    colformat_double(j = "Pct",   digits = 1) %>%
    theme_booktabs() %>%
    hline(i = ~ !is.na(Site), border = fp_border(color = "grey80", width = 0.5)) %>%
    bold(part = "header") %>%
    set_caption(caption_txt) %>%
    autofit()

  ft
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
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# Word (best for pasting into a thesis chapter, keeps italics)
save_as_docx(
  `Table 1. Top 3 taxa by biomass` = ft_biomass,
  `Table 2. Top 3 taxa by density` = ft_density,
  path = "outputs/figures/Top3_taxa_tables.docx",
  pr_section = prop_section(
    page_size = page_size(orient = "portrait"),
    type = "continuous"
  )
)

# HTML (quick preview / repository browsing)
save_as_html(
  `Top 3 taxa by biomass` = ft_biomass,
  `Top 3 taxa by density` = ft_density,
  path = "outputs/figures/Top3_taxa_tables.html"
)

# PNG (needs webshot2 + chromote installed)
# save_as_image(ft_biomass, path = "outputs/figures/Table1_top3_biomass.png")
# save_as_image(ft_density, path = "outputs/figures/Table2_top3_density.png")

# PowerPoint, if your figure repo takes editable vector objects
# save_as_pptx(`Top3 biomass` = ft_biomass,
#              `Top3 density` = ft_density,
#              path = "outputs/figures/Top3_taxa_tables.pptx")
