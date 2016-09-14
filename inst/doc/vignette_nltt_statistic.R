## ------------------------------------------------------------------------
library(nLTT) # nolint

## ------------------------------------------------------------------------
options(warn = 2)

## ------------------------------------------------------------------------
set.seed(42)

## ------------------------------------------------------------------------
species_tree <- ape::rcoal(n = 10)

## ----fig.width = 7, fig.height = 4---------------------------------------
par(mfrow = c(1, 2))
ape::plot.phylo(species_tree)
nLTT::nltt_plot(species_tree)
par(mfrow = c(1, 1))

## ------------------------------------------------------------------------
posteriors <- c(
  ape::rcoal(n = 10), 
  ape::rcoal(n = 10), 
  species_tree, 
  ape::rcoal(n = 10)
)

## ----fig.width = 7, fig.height = 10--------------------------------------
par(mfrow = c(4, 2))
for (p in posteriors) {
  ape::plot.phylo(p, cex = 1)
  nLTT::nltt_plot(p)
}
par(mfrow = c(1, 1))

## ----fig.width = 7, fig.height = 7---------------------------------------
# Set the plot area
nLTT::nltt_plot(species_tree)
# Set the background to light grey
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = colors()[240])
# Draw the species tree
nLTT::nltt_lines(species_tree, lwd = 5)
# Draw the posteriors
col_index <- 1
for (p in posteriors) {
  nLTT::nltt_lines(p, lwd = 3, lty = 2, col = grDevices::heat.colors(n = 4)[col_index])
  col_index <- col_index + 1
}

## ------------------------------------------------------------------------
nltt_stats_exact  <- rep(x = 0, times = length(posteriors))
nltt_stats_approx <- rep(x = 0, times = length(posteriors))
i <- 1
for (p in posteriors) {
  nltt_stats_exact[i]  <- nLTTstat_exact(species_tree, p)
  nltt_stats_approx[i] <- nLTTstat(species_tree, p)
  i <- i + 1
}

## ------------------------------------------------------------------------
nltt_stats <- data.frame(
  id = seq(1, length(nltt_stats_exact)),
  nltt_stat_exact = nltt_stats_exact, 
  nltt_stat_approx = nltt_stats_approx
)
knitr::kable(nltt_stats)

## ------------------------------------------------------------------------
df <- reshape2::melt(
  nltt_stats, 
  id.vars = c("id"), 
  measure.vars = c( "nltt_stat_exact", "nltt_stat_approx")
)
names(df) <- c("id", "method", "nltt_stat")
df$id <- as.factor(df$id)
df$method <- plyr::revalue(
  df$method, 
  c("nltt_stat_exact" = "exact", "nltt_stat_approx" = "approx")
)

## ------------------------------------------------------------------------
ggplot2::ggplot(
  data = df, 
  ggplot2::aes(x = df$method, y = df$nltt_stat)
) + ggplot2::geom_boxplot(
) + ggplot2::geom_point(color = df$id
) + ggplot2::scale_y_continuous(name = "nLTT statistic"
) + ggplot2::scale_x_discrete(name = "Method"
) + ggplot2::ggtitle("Posterior nLTT statistic distribution")

