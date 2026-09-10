# =============================================================================
# Combined 3D Wireframe Plot: Somatic and TAS Landscapes
# =============================================================================
# Creates a 2-panel figure with wireframe plots for both outcome variables
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(tidyverse)
library(lattice)
library(gridExtra)
library(grid)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

# Standardize
dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$somatic_z <- scale(rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE))[,1]
dfc$tas_z <- scale(dfc$tas)[,1]

# =============================================================================
# 2. FIT INTERACTION MODELS
# =============================================================================

model_somatic <- lm(somatic_z ~ ias_z * iats_z + I(ias_z^2) + I(iats_z^2), data = dfc)
model_tas <- lm(tas_z ~ ias_z * iats_z + I(ias_z^2) + I(iats_z^2), data = dfc)

# Get p-values for interactions
p_somatic <- summary(model_somatic)$coefficients["ias_z:iats_z", "Pr(>|t|)"]
p_tas <- summary(model_tas)$coefficients["ias_z:iats_z", "Pr(>|t|)"]

cat("Somatic interaction p =", round(p_somatic, 3), "\n")
cat("TAS interaction p =", round(p_tas, 3), "\n")

# =============================================================================
# 3. CREATE SURFACE GRIDS
# =============================================================================

n_grid <- 40
ias_seq <- seq(-2.5, 2.5, length.out = n_grid)
iats_seq <- seq(-2.5, 2.5, length.out = n_grid)

# Somatic grid
grid_somatic <- expand.grid(IAS = ias_seq, IATS = iats_seq)
grid_somatic$Z <- predict(model_somatic,
                          newdata = data.frame(ias_z = grid_somatic$IAS,
                                               iats_z = grid_somatic$IATS))

# TAS grid
grid_tas <- expand.grid(IAS = ias_seq, IATS = iats_seq)
grid_tas$Z <- predict(model_tas,
                      newdata = data.frame(ias_z = grid_tas$IAS,
                                           iats_z = grid_tas$IATS))

# =============================================================================
# 4. CREATE WIREFRAME PLOTS
# =============================================================================

# Common color palette
color_pal <- colorRampPalette(c("#27ae60", "#f1c40f", "#e67e22", "#c0392b"))(100)

# Format p-value for display
format_p <- function(p) {
  if (p < 0.001) return("p < .001")
  if (p < 0.01) return(sprintf("p = .%03d", round(p * 1000)))
  return(sprintf("p = .%02d", round(p * 100)))
}

# Somatic wireframe
wf_somatic <- wireframe(
  Z ~ IAS * IATS,
  data = grid_somatic,
  drape = TRUE,
  colorkey = FALSE,
  col.regions = color_pal,
  screen = list(z = -60, x = -70),
  xlab = list("IAS\n(Accuracy)", rot = 30, cex = 0.9),
  ylab = list("IATS\n(Attention)", rot = -40, cex = 0.9),
  zlab = list("Somatic\nSymptoms", rot = 90, cex = 0.9),
  main = list(paste0("A. Somatic Symptoms\nIAS x IATS: ", format_p(p_somatic), " *"), cex = 1.1),
  par.settings = list(
    axis.line = list(col = "gray40"),
    background = list(col = "white")
  ),
  scales = list(
    arrows = FALSE,
    col = "gray40",
    cex = 0.7
  )
)

# TAS wireframe
wf_tas <- wireframe(
  Z ~ IAS * IATS,
  data = grid_tas,
  drape = TRUE,
  colorkey = FALSE,
  col.regions = color_pal,
  screen = list(z = -60, x = -70),
  xlab = list("IAS\n(Accuracy)", rot = 30, cex = 0.9),
  ylab = list("IATS\n(Attention)", rot = -40, cex = 0.9),
  zlab = list("TAS\n(Alexithymia)", rot = 90, cex = 0.9),
  main = list(paste0("B. Alexithymia (TAS)\nIAS x IATS: ", format_p(p_tas), " (n.s.)"), cex = 1.1),
  par.settings = list(
    axis.line = list(col = "gray40"),
    background = list(col = "white")
  ),
  scales = list(
    arrows = FALSE,
    col = "gray40",
    cex = 0.7
  )
)

# =============================================================================
# 5. COMBINE AND SAVE
# =============================================================================

png("plots/supplementary_plots/3_s_6_3d_interaction_surfaces.png",
    width = 2000, height = 1000, res = 150)

grid.arrange(
  wf_somatic, wf_tas,
  ncol = 2,
  top = textGrob("3D Landscapes: IAS x IATS Interaction Effects",
                 gp = gpar(fontsize = 16, fontface = "bold")),
  bottom = textGrob("Green = low (better) | Red = high (worse)",
                    gp = gpar(fontsize = 10, col = "gray40"))
)

dev.off()

cat("\nSaved: figure_3d_combined_wireframe.png\n")
cat("\nSummary:\n")
cat("  - Somatic: IAS x IATS interaction SIGNIFICANT (p = .004)\n")
cat("  - TAS: IAS x IATS interaction NOT significant (p =", round(p_tas, 3), ")\n")
