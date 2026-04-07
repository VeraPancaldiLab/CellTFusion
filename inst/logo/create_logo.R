#!/usr/bin/env Rscript
# Create CellTFusion hex sticker logo
# Dense cell clusters with nuclei and TF network edges

library(ggplot2)
library(hexSticker)
library(showtext)
library(sysfonts)
library(ggforce)

font_add_google("Montserrat", "montserrat")
font_add("LiberationMono", regular = "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf")
showtext_auto()

set.seed(42)

# --- Dense cell clusters filling the entire hex ---

# All cells placed in y range [0.0, 0.75] to leave top clear for the name

# Cluster 1 (left) - immune cells (blue)
n1 <- 16
c1 <- data.frame(
  x = rnorm(n1, 0.18, 0.10), y = rnorm(n1, 0.55, 0.08),
  r = runif(n1, 0.028, 0.060), type = "cell_immune"
)

# Cluster 2 (right) - tumor cells (red)
n2 <- 15
c2 <- data.frame(
  x = rnorm(n2, 0.82, 0.09), y = rnorm(n2, 0.20, 0.07),
  r = runif(n2, 0.028, 0.060), type = "cell_tumor"
)

# Cluster 3 (centre) - fusion zone (purple)
n3 <- 14
c3 <- data.frame(
  x = rnorm(n3, 0.50, 0.08), y = rnorm(n3, 0.38, 0.07),
  r = runif(n3, 0.025, 0.055), type = "cell_mixed"
)

# Cluster 4 (upper-right) - satellite immune
n4 <- 10
c4 <- data.frame(
  x = rnorm(n4, 0.75, 0.06), y = rnorm(n4, 0.58, 0.05),
  r = runif(n4, 0.022, 0.048), type = "cell_immune"
)

# Cluster 5 (lower-left) - satellite tumor
n5 <- 10
c5 <- data.frame(
  x = rnorm(n5, 0.25, 0.06), y = rnorm(n5, 0.18, 0.05),
  r = runif(n5, 0.022, 0.048), type = "cell_tumor"
)

# Cluster 6 (far left mid) - edge immune
n6 <- 8
c6 <- data.frame(
  x = rnorm(n6, 0.05, 0.04), y = rnorm(n6, 0.38, 0.07),
  r = runif(n6, 0.020, 0.042), type = "cell_immune"
)

# Cluster 7 (far right mid) - edge tumor
n7 <- 8
c7 <- data.frame(
  x = rnorm(n7, 0.95, 0.04), y = rnorm(n7, 0.38, 0.07),
  r = runif(n7, 0.020, 0.042), type = "cell_tumor"
)

# Cluster 8 (upper band) - mixed filling under name
n8a <- 8
c8a <- data.frame(
  x = rnorm(n8a, 0.50, 0.12), y = rnorm(n8a, 0.65, 0.04),
  r = runif(n8a, 0.020, 0.040), type = "cell_mixed"
)

# Cluster 9 (bottom centre) - mixed near bottom
n8b <- 8
c8b <- data.frame(
  x = rnorm(n8b, 0.50, 0.10), y = rnorm(n8b, 0.08, 0.04),
  r = runif(n8b, 0.020, 0.040), type = "cell_mixed"
)

# Cluster 10 (centre-right bridge)
n9 <- 7
c9 <- data.frame(
  x = rnorm(n9, 0.65, 0.05), y = rnorm(n9, 0.48, 0.04),
  r = runif(n9, 0.018, 0.040), type = "cell_immune"
)

# Cluster 11 (centre-left bridge)
n10 <- 7
c10 <- data.frame(
  x = rnorm(n10, 0.35, 0.05), y = rnorm(n10, 0.28, 0.04),
  r = runif(n10, 0.018, 0.040), type = "cell_tumor"
)

# Scattered filler cells (below y=0.70)
n_fill <- 20
c_fill <- data.frame(
  x = runif(n_fill, -0.02, 1.02), y = runif(n_fill, 0.02, 0.68),
  r = runif(n_fill, 0.014, 0.032),
  type = sample(c("cell_immune", "cell_tumor", "cell_mixed"), n_fill, replace = TRUE)
)

cells <- rbind(c1, c2, c3, c4, c5, c6, c7, c8a, c8b, c9, c10, c_fill)

# Nuclei
nuclei <- data.frame(
  x = cells$x + runif(nrow(cells), -0.005, 0.005),
  y = cells$y + runif(nrow(cells), -0.005, 0.005),
  r = cells$r * runif(nrow(cells), 0.32, 0.48),
  type = paste0(cells$type, "_nuc")
)

# TF-network edges
set.seed(77)
n_edges <- 40
idx_from <- sample(seq_len(nrow(cells)), n_edges, replace = TRUE)
idx_to   <- sample(seq_len(nrow(cells)), n_edges, replace = TRUE)
edges <- data.frame(
  x = cells$x[idx_from], y = cells$y[idx_from],
  xend = cells$x[idx_to], yend = cells$y[idx_to]
)
edges <- edges[!(edges$x == edges$xend & edges$y == edges$yend), ]

# --- Build subplot ---
p <- ggplot() +
  geom_segment(
    data = edges, aes(x = x, y = y, xend = xend, yend = yend),
    color = "#A0C4FF", alpha = 0.15, linewidth = 0.18
  ) +
  geom_circle(
    data = cells, aes(x0 = x, y0 = y, r = r, fill = type),
    color = "white", linewidth = 0.2, alpha = 0.75
  ) +
  geom_circle(
    data = nuclei, aes(x0 = x, y0 = y, r = r, fill = type),
    color = NA, alpha = 0.85
  ) +
  scale_fill_manual(values = c(
    "cell_immune"     = "#4DA6FF",
    "cell_tumor"      = "#FF6B6B",
    "cell_mixed"      = "#9B59B6",
    "cell_immune_nuc" = "#1A5FA0",
    "cell_tumor_nuc"  = "#A83232",
    "cell_mixed_nuc"  = "#5B2D80"
  )) +
  coord_fixed(xlim = c(-0.08, 1.08), ylim = c(-0.05, 0.80)) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

# --- Create hexSticker ---
sticker(
  subplot    = p,
  s_x        = 1.0,
  s_y        = 0.85,
  s_width    = 1.75,
  s_height   = 1.25,
  package    = "CellTFusion",
  p_size     = 15,
  p_y        = 1.48,
  p_x        = 1.0,
  p_color    = "#FFFFFF",
  p_family   = "LiberationMono",
  p_fontface = "bold",
  h_fill     = "#1B2A4A",
  h_color    = "#4DA6FF",
  h_size     = 1.8,
  filename   = "man/figures/logo.png",
  dpi        = 300
)

cat("Logo created at man/figures/logo.png\n")
