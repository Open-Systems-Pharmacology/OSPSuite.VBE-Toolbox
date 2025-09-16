library(ggplot2)

plasmaDf <- refAndTestSimulationsInVirtualPopulation
plasmaDf$conc <- toUnit(
  quantityOrDimension = ospDimensions$`Concentration (mass)`,
  values = plasmaDf$conc,
  targetUnit = ospUnits$`Concentration [mass]`$`µg/l`
)

mmdf <- NULL

q05 <- function(x) {
  quantile(x = x, 0.05)
}
q95 <- function(x) {
  quantile(x = x, 0.95)
}
fnVec <- c("median" = median, "min" = q05, "max" = q95)
for (fn in names(fnVec)) {
  df <- aggregate(
    x = plasmaDf$conc,
    by = list(
      plasmaDf$time,
      plasmaDf$drug
    ),
    fnVec[[fn]]
  )
  names(df) <- c("time", "drug", fn)
  mmdf[["time"]] <- df[["time"]]
  mmdf[["drug"]] <- df[["drug"]]
  mmdf[[fn]] <- df[[fn]]
}
mmdf <- as.data.frame(mmdf)
mmdf$drug <- as.character(mmdf$drug)
mmdf$drug[mmdf$drug == "R"] <- "Petrolatum"
mmdf$drug[mmdf$drug == "T1"] <- "Ethylene glycol"

plt <- ggplot() +
  geom_ribbon(
    data = mmdf[mmdf$drug == "Petrolatum", ],
    mapping = aes(
      x = time / 60,
      y = (median) * 1e5,
      ymin = (min) * 1e5,
      ymax = (max) * 1e5,
      fill = drug
    ), 
    alpha = 0.25
  ) +
  theme_minimal()

pkData <- read.csv(pkDataFilePath)
pkData$drug <- "Petrolatum"
pkData$id <- as.factor(pkData$id)
plt <- plt + geom_point(data = pkData, mapping = aes(
  x = time, y = (outputValues) * 1e5, color = id,
  shape = "Individual"
))

plt <- plt + xlab("Time (h)") + ylab(expression(Testosterone ~ plasma ~ concentration %*% 10^5 ~ (µg / l)))
plt <- plt + guides(color = "none")
plt <- plt + theme(
  panel.background = element_rect(fill = "white", colour = "white"),
  plot.background = element_rect(fill = "white", colour = "white"),
  legend.position = "bottom", legend.direction = "horizontal",
  axis.text = element_text(size = 12), axis.title = element_text(size = 14),
  legend.text = element_text(size = 12), legend.title = element_text(size = 14),
  strip.text = element_text(size = 14)
)

plt <- plt + scale_fill_manual(name = "Petrolatum\nDose: 2 µg/cm²", labels = "Simulated", values = "#ff0000")
plt <- plt + scale_shape_manual(name = " \n ", labels = "Observed", values = 16)
plt <- plt + guides(
  fill = guide_legend(order = 1), # Set the order for color legend
  shape = guide_legend(order = 2) # Set the order for linetype legend
)
methods::show(plt)


plasmaDf <- refAndTestSimulationsInVirtualPopulation
plasmaDf$conc <- toUnit(
  quantityOrDimension = ospDimensions$`Concentration (mass)`,
  values = plasmaDf$conc,
  targetUnit = ospUnits$`Concentration [mass]`$`µg/l`
)

vbemmdf <- mmdf
vbemmdf$drug <- factor(vbemmdf$drug, levels = c("Petrolatum", "Ethylene glycol"))
vbeplt <- ggplot() +
  geom_ribbon(
    data = vbemmdf,
    mapping = aes(
      x = time / 60,
      y = (median) * 1e5,
      ymin = (min) * 1e5,
      ymax = (max) * 1e5,
      fill = drug
    ), 
    alpha = 0.25
  )

vbeplt <- vbeplt + geom_line(
  data = vbemmdf,
  mapping = aes(
    x = time / 60,
    y = (median) * 1e5,
    color = drug
  ),
  linewidth = 2
) + theme_minimal()

vbeplt <- vbeplt + theme(
  panel.background = element_rect(fill = "white", colour = "white"),
  plot.background = element_rect(fill = "white", colour = "white"),
  legend.position = "bottom", legend.direction = "horizontal",
  axis.text = element_text(size = 12), axis.title = element_text(size = 14),
  legend.text = element_text(size = 12), legend.title = element_text(size = 14),
  strip.text = element_text(size = 14)
)

vbeplt <- vbeplt + xlab("Time (h)") + ylab(expression(Testosterone ~ plasma ~ concentration %*% 10^5 ~ (µg / l)))

vbeplt <- vbeplt + scale_fill_manual(name = "Formulation", labels = c("Petrolatum\nDose: 2 µg/cm²", "Ethylene glycol\nDose: 2 µg/cm²"), values = c("#ff0000", "#1e90ff"))
vbeplt <- vbeplt + scale_color_manual(name = "Formulation", labels = c("Petrolatum\nDose: 2 µg/cm²", "Ethylene glycol\nDose: 2 µg/cm²"), values = c("#ff0000", "#1e90ff"))
methods::show(vbeplt)
