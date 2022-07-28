###########################################################

test_that(desc = "testing mspace, general behavior for landmark data", code = {
  data("tails")
  shapes <- tails$shapes

  msp <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- class(msp) == "mspace"
  result2 <- msp$datype == "landm"
  result3 <- all(names(msp) %in% c("x", "rotation","center","datype","ordtype","plotinfo"))

  msp2 <- mspace(geomorph::two.d.array(shapes), p = nrow(shapes), k = ncol(shapes),
                 mag = 0.7, axes = c(1,2), plot = FALSE)

  result4 <- nrow(msp$x) == dim(shapes)[3]
  result5 <- nrow(msp$x) == nrow(msp2$x)

  mspace(shapes, mag = 0.7, axes = c(1,2))
  result6 <- .Device != "null device"

  expect_true(all(result1,result2,result3,result4,result5,result6))
  dev.off()
})


test_that(desc = "testing mspace, general behavior for Fourier data", code = {
  data("shells")
  shapes <- shells$shapes

  msp <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- class(msp) == "mspace"
  result2 <- msp$datype == "fcoef"
  result3 <- all(names(msp) %in% c("x", "rotation","center","datype","ordtype","plotinfo"))

  msp2 <- mspace(shapes$coe, mag = 0.7, axes = c(1,2), plot = FALSE)

  result4 <- nrow(msp$x) == dim(shapes)[1]
  result5 <- nrow(msp$x) == nrow(msp2$x)

  mspace(shapes, mag = 0.7, axes = c(1,2))
  result6 <- .Device != "null device"

  expect_true(all(result1,result2,result3,result4,result5,result6))
  dev.off()
})


test_that(desc = "testing mspace, prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  links <- tails$links
  tree <- tails$tree
  template <- tails$template

  pca <- prcomp(geomorph::two.d.array(shapes))
  bgpca <- bg_prcomp(geomorph::two.d.array(shapes), species)
  phypca <- phy_prcomp(geomorph::two.d.array(expected_shapes(shapes, species)), tree)
  pls <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = model.matrix(~species)[,-1])

  msp1 <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordtype == "prcomp"
  result2 <- all(msp1$x == pca$x)
  result3 <- all(msp1$rotation == pca$rotation)
  result4 <- all(msp1$center == pca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  pca <- prcomp(geomorph::two.d.array(shapes))

  msp1 <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordtype == "prcomp"
  result2 <- all(msp1$x == pca$x)
  result3 <- all(msp1$rotation == pca$rotation)
  result4 <- all(msp1$center == pca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, bg_prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  bgpca <- bg_prcomp(geomorph::two.d.array(shapes), species)

  msp1 <- mspace(shapes, FUN = bg_prcomp, groups = species, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordtype == "bg_prcomp"
  result2 <- all(msp1$x == bgpca$x)
  result3 <- all(msp1$rotation == bgpca$rotation)
  result4 <- all(msp1$center == bgpca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, phy_prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  tree <- tails$tree
  phypca <- phy_prcomp(geomorph::two.d.array(expected_shapes(shapes, species)), tree)

  msp1 <- mspace(expected_shapes(shapes, species), FUN = phy_prcomp, tree = tree, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordtype == "phy_prcomp"
  result2 <- all(msp1$x == phypca$x)
  result3 <- all(msp1$rotation == phypca$rotation)
  result4 <- all(msp1$center == phypca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, pls_shapes", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  pls <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = model.matrix(~species)[,-1])

  msp1 <- mspace(shapes, FUN = pls_shapes, X = model.matrix(~species)[,-1],
                 mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordtype == "pls_shapes"
  result2 <- all(msp1$x == pls$x)
  result3 <- all(msp1$rotation == pls$rotation)
  result4 <- all(msp1$center == pls$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, phylogenetic pls_shapes", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  tree <- tails$tree
  sp_sizes <- cbind(tapply(tails$sizes,species,mean))
  phypls <- pls_shapes(shapes = geomorph::two.d.array(expected_shapes(shapes, species)),
                       X = sp_sizes, tree = tree)

  msp1 <- mspace(expected_shapes(shapes,species), FUN = pls_shapes, X = sp_sizes, tree = tree,
                 mag = 0.7, axes = 1, plot = FALSE)
  result1 <- msp1$ordtype == "phy_pls_shapes"
  result2 <- all(msp1$x == phypls$x)
  result3 <- all(msp1$rotation == phypls$rotation)
  result4 <- all(msp1$center == phypls$center)

  expect_true(all(result1,result2,result3,result4))

})


###########################################################

test_that(desc = "testing proj_shapes", code = {
  data("tails")
  shapes <- tails$shapes

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_shapes(shapes = shapes)

  result1 <- all(names(msp1) %in% c("x", "rotation","center","datype","ordtype","plotinfo"))

  n <- 10
  msp2 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_shapes(shapes = shapes[,,1:n])

  result2 <- nrow(msp2$x) == n

  expect_true(all(result1,result2))

})


###########################################################

test_that(desc = "testing proj_consensus", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  sp_shapes <- expected_shapes(shapes, species)

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_consensus(shapes = sp_shapes)

  result1 <- all(names(msp1) %in% c("x", "rotation","center","datype","ordtype","plotinfo","gr_centroids"))

  sc1 <- msp1$gr_centroids
  sc2 <- apply(msp1$x,2,tapply,species,mean)
  result2 <- all(round(sc1,10) == round(sc2,10))

  sc3 <- proj_consensus(msp1, shapes = sp_shapes, pipe = FALSE)
  result3 <- all(round(sc3,10) == round(sc2,10))

  expect_true(all(result1,result2,result3))
})


###########################################################

test_that(desc = "testing proj_groups, hulls", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(groups = species)

  result1 <- all(names(msp1) %in% c("x", "rotation","center","datype","ordtype","plotinfo","gr_class"))
  result2 <- all(msp1$gr_class == species)

  expect_true(all(result1,result2))
})


test_that(desc = "testing proj_groups", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(groups = species, ellipse = TRUE)

  result1 <- all(names(msp1) %in% c("x", "rotation","center","datype","ordtype","plotinfo","gr_class"))
  result2 <- all(msp1$gr_class == species)

  expect_true(all(result1,result2))
})


###########################################################

test_that(desc = "testing proj_phylogeny", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  sp_shapes <- expected_shapes(shapes, species)
  tree <- tails$tree

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_consensus(shapes = sp_shapes) %>%
    proj_phylogeny(tree = tree)

  result1 <- all(names(msp1) %in% c("x", "rotation","center","datype","ordtype","plotinfo","gr_centroids","phylo_scores","phylo"))

  sc1 <- round(msp1$phylo_scores[rownames(msp1$phylo_scores),],10)
  sc2 <- round(rbind(msp1$gr_centroids, apply(msp1$gr_centroids,2,phytools::fastAnc, tree=tree)),10)[rownames(msp1$phylo_scores),]
  result2 <- all(round(sc1,10) == round(sc2,10))

  sc3 <- proj_phylogeny(mspace = msp1, tree = tree, pipe = FALSE)[rownames(msp1$phylo_scores),]
  result3 <- all(round(sc3,10) == round(sc2,10))

  result4 <- all(msp1$phylo$tip.labels == tree$tip.labels)

  expect_true(all(result1,result2,result3,result4))
})


###########################################################

test_that(desc = "testing proj_axis", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  sp_shapes <- expected_shapes(shapes, species)
  pca <- prcomp(geomorph::two.d.array(shapes))

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_axis(obj = pca)

  result1 <- all(names(msp1) %in% c("x", "rotation","center","datype","ordtype","plotinfo","shape_axis"))

  sc1 <- proj_axis(obj = pca, axis = 1, mspace = msp1, pipe = FALSE)
  sc2 <- range(msp1$x[,1])
  result2 <- all(round(sc1[,1],10) == round(sc2,10))

  expect_true(all(result1,result2))

})

