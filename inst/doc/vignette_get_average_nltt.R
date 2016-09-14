# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
newick1 <- "((A:1,B:1):2,C:3);"
newick2 <- "((A:2,B:2):1,C:3);"
phylogeny1 <- ape::read.tree(text = newick1)
phylogeny2 <- ape::read.tree(text = newick2)

# ------------------------------------------------------------------------
ape::plot.phylo(phylogeny1)
ape::add.scale.bar() #nolint

# ------------------------------------------------------------------------
nLTT::nltt_plot(phylogeny1, ylim = c(0, 1))

# ------------------------------------------------------------------------
nLTT::get_phylogeny_nltt_matrix(phylogeny1)

# ------------------------------------------------------------------------
ape::plot.phylo(phylogeny2)
ape::add.scale.bar() #nolint

# ------------------------------------------------------------------------
nLTT::nltt_plot(phylogeny2, ylim = c(0, 1))

# ------------------------------------------------------------------------
nLTT::get_phylogeny_nltt_matrix(phylogeny2)

# ------------------------------------------------------------------------
nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny1),
  dt = 0.20,
  step_type = "upper"
)

# ------------------------------------------------------------------------
nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny2), 
  dt = 0.20,
  step_type = "upper"
)

# ------------------------------------------------------------------------
nLTT::get_average_nltt_matrix(c(phylogeny1, phylogeny2), dt = 0.20)

# ------------------------------------------------------------------------
nLTT::get_average_nltt(c(phylogeny1, phylogeny2), dt = 0.20, plot_nltts = TRUE)

# ------------------------------------------------------------------------
newick1 <- "((A:1,B:1):1,(C:1,D:1):1);"
newick2 <- paste0("((((XD:1,ZD:1):1,CE:2):1,(FE:2,EE:2):1):4,((AE:1,BE:1):1,",
  "(WD:1,YD:1):1):5);"
)
phylogeny1 <- ape::read.tree(text = newick1)
phylogeny2 <- ape::read.tree(text = newick2)

# ------------------------------------------------------------------------
ape::plot.phylo(phylogeny1)
ape::add.scale.bar() #nolint

# ------------------------------------------------------------------------
nLTT::nltt_plot(phylogeny1, ylim = c(0, 1))

# ------------------------------------------------------------------------
nLTT::get_phylogeny_nltt_matrix(phylogeny2)

# ------------------------------------------------------------------------
ape::plot.phylo(phylogeny2)
ape::add.scale.bar() #nolint

# ------------------------------------------------------------------------
nLTT::nltt_plot(phylogeny2, ylim = c(0, 1))

# ------------------------------------------------------------------------
nLTT::get_phylogeny_nltt_matrix(phylogeny2)

# ------------------------------------------------------------------------
nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny1),
  dt = 0.20,
  step_type = "upper"
)

# ------------------------------------------------------------------------
nLTT::stretch_nltt_matrix(
  nLTT::get_phylogeny_nltt_matrix(phylogeny2), 
  dt = 0.20,
  step_type = "upper"
)

# ------------------------------------------------------------------------
nLTT::get_average_nltt_matrix(c(phylogeny1, phylogeny2), dt = 0.20)

# ------------------------------------------------------------------------
nLTT::get_average_nltt(
  c(phylogeny1, phylogeny2), 
  dt = 0.20,
  plot_nltts = TRUE
)

