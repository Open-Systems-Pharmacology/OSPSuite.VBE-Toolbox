library(ggplot2)

mmdf <- NULL
q05 <- function(x) {
  quantile(x = x, 0.05)
}
q95 <- function(x) {
  quantile(x = x, 0.95)
}
fnVec <- c("median" = median, "min" = min, "max" = max)
for (fn in names(fnVec)) {
  df <- aggregate(
    x = refAndTestSimulationsInVirtualPopulation$conc,
    by = list(
      refAndTestSimulationsInVirtualPopulation$time,
      refAndTestSimulationsInVirtualPopulation$drug
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

#################

pkData <- read.csv(pkDataFilePath)
pkData$drug <- "R"
pkData$id <- as.factor(pkData$id)

plt <- ggplot() +
  geom_ribbon(
    data = mmdf[mmdf$drug == "R", ],
    mapping = aes(
      x = time / 60,
      y = (median),
      ymin = (min),
      ymax = (max),
      fill = drug
    ),
    alpha = 0.25
  ) +
  theme_minimal()
plt <- plt + geom_point(data = na.omit(pkData), mapping = aes(
  x = time / 60,
  y = (outputValues),
  color = id,
  shape = "Individual"
))
plt <- plt + xlab("Time (h)") + ylab(expression(Bupropion ~ plasma ~ concentration ~ (µmol / l)))
plt <- plt + guides(color = "none")
plt <- plt + theme(
  legend.position = "bottom", legend.direction = "horizontal",
  panel.background = element_rect(colour = "white", fill = "white"),
  plot.background = element_rect(colour = "white", fill = "white"),
  axis.text = element_text(size = 12), axis.title = element_text(size = 14),
  legend.text = element_text(size = 12), legend.title = element_text(size = 14),
  strip.text = element_text(size = 14)
)

plt <- plt + scale_fill_manual(name = "Sustained release", labels = "Simulated", values = "#ff0000")
plt <- plt + scale_shape_manual(name = "", labels = "Observed", values = 16)
plt <- plt + guides(
  fill = guide_legend(order = 1), # Set the order for color legend
  shape = guide_legend(order = 2) # Set the order for linetype legend
)
methods::show(plt)

#################

TpkData <- read.csv(file.path(subfolder, "pk_data_test.csv"))
TpkData$id <- as.numeric(as.factor(TpkData$id))
TpkData <- TpkData[TpkData$frm == "ER", ]
TpkData$time_h <- TpkData$time_h
TpkData$id <- as.numeric(as.factor(TpkData$id))
TpkData$frm <- "T"
names(TpkData) <- c("time", "id", "outputValues", "drug")
TpkData$id <- as.factor(TpkData$id)

Tplt <- ggplot() +
  geom_ribbon(
    data = mmdf[mmdf$drug == "T1", ],
    mapping = aes(
      x = time / 60,
      y = (median),
      ymin = (min),
      ymax = (max),
      fill = drug
    ),
    alpha = 0.25
  ) +
  theme_minimal()
Tplt <- Tplt + geom_point(data = na.omit(TpkData), mapping = aes(x = time, y = (outputValues), color = id, shape = "Individual"))
Tplt <- Tplt + xlab("Time (h)") + ylab(expression(Bupropion ~ plasma ~ concentration ~ (µmol / l)))
Tplt <- Tplt + guides(color = "none")
Tplt <- Tplt + theme(
  panel.background = element_rect(fill = "white", colour = "white"),
  plot.background = element_rect(fill = "white", colour = "white"),
  legend.position = "bottom", legend.direction = "horizontal",
  axis.text = element_text(size = 12), axis.title = element_text(size = 14),
  legend.text = element_text(size = 12), legend.title = element_text(size = 14),
  strip.text = element_text(size = 14)
)

Tplt <- Tplt + scale_fill_manual(name = "Extended release", labels = "Simulated", values = "#1e90ff")
Tplt <- Tplt + scale_shape_manual(name = "", labels = "Observed", values = 16)
Tplt <- Tplt + guides(
  fill = guide_legend(order = 1), # Set the order for color legend
  shape = guide_legend(order = 2) # Set the order for linetype legend
)
methods::show(Tplt)

###################

vbemmdf <- mmdf
vbemmdf$drug <- factor(vbemmdf$drug, levels = c("R", "T1"))
vbeplt <- ggplot() +
  geom_ribbon(
    data = vbemmdf,
    mapping = aes(
      x = time / 60,
      y = (median),
      ymin = (min),
      ymax = (max),
      fill = drug
    ), alpha = 0.25
  )
vbeplt <- vbeplt + geom_line(
  data = vbemmdf,
  mapping = aes(
    x = time / 60,
    y = (median),
    color = drug
  ),
  linewidth = 2
) + theme_minimal()
vbeplt <- vbeplt + theme(
  legend.position = "bottom", legend.direction = "horizontal",
  panel.background = element_rect(color = "white", fill = "white"),
  plot.background = element_rect(color = "white", fill = "white"),
  axis.text = element_text(size = 12), axis.title = element_text(size = 14),
  legend.text = element_text(size = 12), legend.title = element_text(size = 14),
  strip.text = element_text(size = 14)
)
vbeplt <- vbeplt + xlab("Time (h)") + ylab(expression(Bupropion ~ plasma ~ concentration ~ (µmol / l)))
vbeplt <- vbeplt + scale_fill_manual(name = "Formulation", labels = c("R" = "Sustained release", "T1" = "Extended release"), values = c("#ff0000", "#1e90ff"))
vbeplt <- vbeplt + scale_color_manual(name = "Formulation", labels = c("R" = "Sustained release", "T1" = "Extended release"), values = c("#ff0000", "#1e90ff"))
methods::show(vbeplt)
