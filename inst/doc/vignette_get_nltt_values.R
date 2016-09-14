# ------------------------------------------------------------------------
library(ape)
library(Hmisc)
library(nLTT) # nolint
library(ggplot2)
library(testit)

# ------------------------------------------------------------------------
newick1 <- "((A:1,B:1):2,C:3);"
newick2 <- "((A:2,B:2):1,C:3);"
phylogeny1 <- ape::read.tree(text = newick1)
phylogeny2 <- ape::read.tree(text = newick2)
phylogenies <- c(phylogeny1, phylogeny2)

# ------------------------------------------------------------------------
plot(phylogeny1)
add.scale.bar() #nolint

# ------------------------------------------------------------------------
nltt_plot(phylogeny1, ylim = c(0, 1))

# ------------------------------------------------------------------------
t <- nLTT::get_phylogeny_nltt_matrix(phylogeny1)
knitr::kable(t)

# ------------------------------------------------------------------------
df <- as.data.frame(nLTT::get_phylogeny_nltt_matrix(phylogeny1))
qplot(
  time, N, data = df, geom = "step", ylim = c(0, 1), direction = "vh",
  main = "NLTT plot of phylogeny 1"
)

# ------------------------------------------------------------------------
plot(phylogeny2)
add.scale.bar() #nolint

# ------------------------------------------------------------------------
nltt_plot(phylogeny2, ylim = c(0, 1))

# ------------------------------------------------------------------------
t <- nLTT::get_phylogeny_nltt_matrix(phylogeny2)
knitr::kable(t)

# ------------------------------------------------------------------------
df <- as.data.frame(nLTT::get_phylogeny_nltt_matrix(phylogeny2))
qplot(
  time, N, data = df, geom = "step", ylim = c(0, 1), direction = "vh",
  main = "NLTT plot of phylogeny 2"
)

# ------------------------------------------------------------------------
t <- nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny1), dt = 0.20, step_type = "upper"
)
knitr::kable(t)

# ------------------------------------------------------------------------
t <- nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny2), dt = 0.20, step_type = "upper"
)
knitr::kable(t)

# ------------------------------------------------------------------------
t <- nLTT::get_average_nltt_matrix(phylogenies, dt = 0.20)
knitr::kable(t)

# ------------------------------------------------------------------------
nLTT::get_average_nltt(phylogenies, dt = 0.20, plot_nltts = TRUE)

# ------------------------------------------------------------------------
t <- nLTT::get_nltt_values(c(phylogeny1, phylogeny2), dt = 0.2)
knitr::kable(t)

# ------------------------------------------------------------------------
df <- nLTT::get_nltt_values(c(phylogeny1, phylogeny2), dt = 0.01)

# ------------------------------------------------------------------------
qplot(
  t, nltt, data = df, geom = "point", ylim = c(0, 1),
  main = "Average nLTT plot of phylogenies", color = id, size = I(0.1)
) +
  stat_summary(fun.data = "mean_cl_boot", color = "red", geom = "smooth")

# ------------------------------------------------------------------------
qplot(t, nltt, data = df, geom = "blank", ylim = c(0, 1),
  main = "Average nLTT plot of phylogenies"
) +
  stat_summary(fun.data = "mean_cl_boot", color = "red", geom = "smooth")

# ------------------------------------------------------------------------
newick1 <- "((A:1,B:1):1,(C:1,D:1):1);"
newick2 <- paste0("((((XD:1,ZD:1):1,CE:2):1,(FE:2,EE:2):1):4,((AE:1,BE:1):1,",
  "(WD:1,YD:1):1):5);"
)
phylogeny1 <- ape::read.tree(text = newick1)
phylogeny2 <- ape::read.tree(text = newick2)
phylogenies <- c(phylogeny1, phylogeny2)

# ------------------------------------------------------------------------
plot(phylogeny1)
add.scale.bar() #nolint

# ------------------------------------------------------------------------
nltt_plot(phylogeny1, ylim = c(0, 1))

# ------------------------------------------------------------------------
t <- nLTT::get_phylogeny_nltt_matrix(phylogeny2)
knitr::kable(t)

# ------------------------------------------------------------------------
plot(phylogeny2)
add.scale.bar() #nolint

# ------------------------------------------------------------------------
nltt_plot(phylogeny2, ylim = c(0, 1))

# ------------------------------------------------------------------------
t <- nLTT::get_phylogeny_nltt_matrix(phylogeny2)
knitr::kable(t)

# ------------------------------------------------------------------------
t <- nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny1), dt = 0.20, step_type = "upper"
)
knitr::kable(t)

# ------------------------------------------------------------------------
t <- nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny2), dt = 0.20, step_type = "upper"
)
knitr::kable(t)

# ------------------------------------------------------------------------
t <- nLTT::get_average_nltt_matrix(phylogenies, dt = 0.20)
knitr::kable(t)

# ------------------------------------------------------------------------
t <- nLTT::get_nltt_values(c(phylogeny1, phylogeny2), dt = 0.2)
knitr::kable(t)

# ------------------------------------------------------------------------
df <- nLTT::get_nltt_values(c(phylogeny1, phylogeny2), dt = 0.01)

# ------------------------------------------------------------------------
qplot(
  t, nltt, data = df, geom = "point", ylim = c(0, 1),
  main = "Average nLTT plot of phylogenies", color = id, size = I(0.1)
) +
  stat_summary(fun.data = "mean_cl_boot", color = "red", geom = "smooth")

# ------------------------------------------------------------------------
qplot(t, nltt, data = df, geom = "blank", ylim = c(0, 1), main = "Average nLTT plot of phylogenies") + 
  stat_summary(fun.data = "mean_cl_boot", color = "red", geom = "smooth")

# ------------------------------------------------------------------------
set.seed(42)
phylogeny1 <- rcoal(10)
phylogeny2 <- rcoal(20)
phylogeny3 <- rcoal(30)
phylogeny4 <- rcoal(40)
phylogeny5 <- rcoal(50)
phylogeny6 <- rcoal(60)
phylogeny7 <- rcoal(70)
phylogenies <- c(phylogeny1, phylogeny2, phylogeny3, 
  phylogeny4, phylogeny5, phylogeny6, phylogeny7)

# ------------------------------------------------------------------------
t <- nLTT::get_nltt_values(phylogenies, dt = 0.2)
knitr::kable(t)

# ------------------------------------------------------------------------
qplot(t, nltt, data = df, geom = "point", ylim = c(0, 1),
  main = "Average nLTT plot of phylogenies", color = id, size = I(0.1)
) +
  stat_summary(fun.data = "mean_cl_boot", color = "red", geom = "smooth")

# ------------------------------------------------------------------------
qplot(t, nltt, data = df, geom = "blank", ylim = c(0, 1), main = "Average nLTT plot of phylogenies") + 
  stat_summary(fun.data = "mean_cl_boot", color = "red", geom = "smooth")

