###########################################################

test_that(desc = "testing mspace, general behavior for landmark data", code = {
  data("tails")
  shapes <- tails$shapes

  msp <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- class(msp) == "mspace"
  result2 <- msp$ordination$datype == "landm"
  result3 <- all(names(msp) %in% c("ordination", "projected", "plotinfo"))
  result4 <- all(c("x", "rotation","center","datype","ordtype") %in% names(msp$ordination))
  result5 <- all(names(msp$projected) %in% c("shapemodels"))

  msp2 <- mspace(geomorph::two.d.array(shapes), p = nrow(shapes), k = ncol(shapes),
                 mag = 0.7, axes = c(1,2), plot = FALSE)

  result6 <- nrow(msp$ordination$x) == dim(shapes)[3]
  result7 <- nrow(msp$ordination$x) == nrow(msp2$ordination$x)

  mspace(shapes, mag = 0.7, axes = c(1,2))
  result8 <- .Device != "null device"

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,result8))
  dev.off()
})


test_that(desc = "testing mspace, general behavior for Fourier data", code = {
  data("shells")
  shapes <- shells$shapes

  msp <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- class(msp) == "mspace"
  result2 <- msp$datype == "fcoef"
  result3 <- all(names(msp) %in% c("ordination", "projected", "plotinfo"))
  result4 <- all(c("x", "rotation","center","datype","ordtype") %in% names(msp$ordination))
  result5 <- all(names(msp$projected) %in% c("shapemodels"))

  msp2 <- mspace(shapes$coe, mag = 0.7, axes = c(1,2), plot = FALSE)

  result6 <- nrow(msp$ordination$x) == dim(shapes)[1]
  result7 <- nrow(msp$ordination$x) == nrow(msp2$ordination$x)

  mspace(shapes, mag = 0.7, axes = c(1,2))
  result8 <- .Device != "null device"

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,result8))
  dev.off()
})


test_that(desc = "testing mspace, prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  pca <- prcomp(geomorph::two.d.array(shapes))

  msp1 <- mspace(shapes, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordination$ordtype == "prcomp"
  result2 <- all(msp1$ordination$x == pca$x)
  result3 <- all(msp1$ordination$rotation == pca$rotation)
  result4 <- all(msp1$ordination$center == pca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, bg_prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  bgpca <- bg_prcomp(geomorph::two.d.array(shapes), species)

  msp1 <- mspace(shapes, FUN = bg_prcomp, groups = species, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordination$ordtype == "bg_prcomp"
  result2 <- all(msp1$ordination$x == bgpca$x)
  result3 <- all(msp1$ordination$rotation == bgpca$rotation)
  result4 <- all(msp1$ordination$center == bgpca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, phy_prcomp", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  tree <- tails$tree
  phypca <- phy_prcomp(geomorph::two.d.array(expected_shapes(shapes, species)), tree)

  msp1 <- mspace(expected_shapes(shapes, species), FUN = phy_prcomp, tree = tree, mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordination$ordtype == "phy_prcomp"
  result2 <- all(msp1$ordination$x == phypca$x)
  result3 <- all(msp1$ordination$rotation == phypca$rotation)
  result4 <- all(msp1$ordination$center == phypca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, phyalign_comp", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  tree <- tails$tree
  paca <- phyalign_comp(geomorph::two.d.array(expected_shapes(shapes, species)), tree)

  msp1 <- mspace(expected_shapes(shapes, species), FUN = phyalign_comp, tree = tree,
                 mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordination$ordtype == "phyalign_comp"
  result2 <- all(msp1$ordination$x == paca$x)
  result3 <- all(msp1$ordination$rotation == paca$rotation)
  result4 <- all(msp1$ordination$center == paca$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, pls_shapes", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  pls <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = model.matrix(~species)[,-1])

  msp1 <- mspace(shapes, FUN = pls_shapes, X = model.matrix(~species)[,-1],
                 mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordination$ordtype == "pls_shapes"
  result2 <- all(msp1$ordination$x == pls$x)
  result3 <- all(msp1$ordination$rotation == pls$rotation)
  result4 <- all(msp1$ordination$center == pls$center)

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
  result1 <- msp1$ordination$ordtype == "phy_pls_shapes"
  result2 <- all(msp1$ordination$x == phypls$x)
  result3 <- all(msp1$ordination$rotation == phypls$rotation)
  result4 <- all(msp1$ordination$center == phypls$center)

  expect_true(all(result1,result2,result3,result4))

})


test_that(desc = "testing mspace, burnaby", code = {
  data("tails")
  shapes <- tails$shapes
  sizes <- log(tails$sizes)
  burn <- burnaby(x = geomorph::two.d.array(shapes), vars = sizes)

  msp1 <- mspace(shapes, FUN = burnaby, vars = sizes,
                 mag = 0.7, axes = c(1,2), plot = FALSE)
  result1 <- msp1$ordination$ordtype == "burnaby"
  result2 <- all(msp1$ordination$x == burn$x)
  result3 <- all(msp1$ordination$rotation == burn$rotation)
  result4 <- all(msp1$ordination$center == burn$center)


  sizeax <- lm(geomorph::two.d.array(shapes) ~ sizes)$coef[2,]
  cent <- colMeans(lm(geomorph::two.d.array(shapes) ~ sizes)$fitted)
  burn2 <- burnaby(x = geomorph::two.d.array(shapes), axmat = sizeax, center = cent)

  msp2 <- mspace(shapes, FUN = burnaby, axmat = sizeax, center = cent,
                 mag = 0.7, axes = c(1,2), plot = FALSE)
  result5 <- msp2$ordination$ordtype == "burnaby"
  result6 <- all(msp2$ordination$x == burn2$x)
  result7 <- all(msp2$ordination$rotation == burn2$rotation)
  result8 <- all(msp2$ordination$center == burn2$center)

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,
                  result8))
})


test_that(desc = "testing mspace, phylogenetic burnaby", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  tree <- tails$tree
  sp_sizes <- cbind(tapply(tails$sizes,species,mean))
  phyburn <- burnaby(x = geomorph::two.d.array(expected_shapes(shapes, species)),
                       vars = sp_sizes, tree = tree)

  msp1 <- mspace(expected_shapes(shapes,species), FUN = burnaby, vars = sp_sizes, tree = tree,
                 mag = 0.7, axes = 1, plot = FALSE)
  result1 <- msp1$ordination$ordtype == "phy_burnaby"
  result2 <- all(msp1$ordination$x == phyburn$x)
  result3 <- all(msp1$ordination$rotation == phyburn$rotation)
  result4 <- all(msp1$ordination$center == phyburn$center)

  sp_sizeax <- lm(geomorph::two.d.array(expected_shapes(shapes, species)) ~ sp_sizes)$coef[2,]
  sp_cent <- colMeans(lm(geomorph::two.d.array(expected_shapes(shapes, species)) ~ sp_sizes)$fitted)
  burn2 <- suppressWarnings(burnaby(x = geomorph::two.d.array(expected_shapes(shapes, species)),
                   axmat = sp_sizeax, center = cent, tree = tree))

  msp2 <- suppressWarnings(mspace(expected_shapes(shapes, species), FUN = burnaby, axmat = sp_sizeax,
                                  center = cent, tree = tree, mag = 0.7, axes = c(1,2), plot = FALSE))
  result5 <- msp2$ordination$ordtype == "phy_burnaby"
  result6 <- all(msp2$ordination$x == burn2$x)
  result7 <- all(msp2$ordination$rotation == burn2$rotation)
  result8 <- all(msp2$ordination$center == burn2$center)

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,
                  result8))
})

###########################################################

test_that(desc = "testing proj_shapes, general behavior", code = {
  data("tails")
  shapes <- tails$shapes

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_shapes(shapes = shapes)
  result1 <- c("scores") %in% names(msp1$projected)
  result2 <- nrow(msp1$projected$scores) == dim(shapes)[3]

  plotinfo1 <- msp1$plotinfo
  result3 <- all(length(plotinfo1$col.points) == dim(shapes)[3],
                 length(plotinfo1$pch.points) == dim(shapes)[3],
                 length(plotinfo1$bg.points) == dim(shapes)[3],
                 length(plotinfo1$cex.points) == dim(shapes)[3])
  result4 <- all(inherits(plotinfo1$col.points, "character"),
                 inherits(plotinfo1$pch.points, "numeric"),
                 inherits(plotinfo1$cex.points, "numeric"),
                 inherits(plotinfo1$bg.points, "character"))

  msp2 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_shapes(shapes = shapes, col = "red", col = 2, bg = "transparent", pch = 16)
  plotinfo2 <- msp2$plotinfo
  result5 <- all(length(plotinfo2$col.points) == dim(shapes)[3],
                 length(plotinfo2$pch.points) == dim(shapes)[3],
                 length(plotinfo2$bg.points) == dim(shapes)[3],
                 length(plotinfo2$cex.points) == dim(shapes)[3])
  result6 <- all(inherits(plotinfo2$col.points, "character"),
                 inherits(plotinfo2$pch.points, "numeric"),
                 inherits(plotinfo2$cex.points, "numeric"),
                 inherits(plotinfo2$bg.points, "character"))

  expect_true(all(result1,result2,result3,result4,result5))

})


test_that(desc = "testing proj_shapes, stacking behavior", code = {
  data("tails")

  library(geomorph)

  shapes <- tails$shapes
  index <- tails$data$type == "DF"

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_shapes(shapes = shapes[,,1])
  result1 <- nrow(msp1$projected$scores) == 1

  msp2 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_shapes(shapes = shapes[,,index], pch = 1) %>%
    proj_shapes(shapes = shapes[,,!index], pch = 2, cex = 2)
  result2 <- nrow(msp2$projected$scores) == dim(shapes)[3]
  result3 <- length(msp2$plotinfo$pch.points) == dim(shapes)[3]
  result4 <- all(msp2$plotinfo$pch.points == c(rep(1, sum(index)), rep(2, sum(!index))))
  result5 <- all(msp2$plotinfo$cex.points == c(rep(1, sum(index)), rep(2, sum(!index))))

  scores <- msp2$projected$scores
  data2d <- geomorph::two.d.array(shapes)
  dec <- 5
  index_x_in_sc <- as.numeric(unlist
                              (apply(round(scores,dec), 1, \(x, y) {
                                which(apply(y, 1, \(z, x){
                                  all(z == x)}, x))},
                                round(stats::prcomp(data2d)$x,dec))))
  result6 <- all(index_x_in_sc == c(which(index), which(!index)))


  expect_true(all(result1,result2,result3,result4,result5,result6))
})

###########################################################

test_that(desc = "testing proj_groups, general behavior", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species

  msp0 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups()
  result1 <- all(c("gr_scores","gr_class") %in% names(msp0$projected))
  result2 <- all(msp0$projected$gr_class == rep(1, dim(shapes)[3]))
  result3 <- all(msp0$projected$gr_scores == msp0$ordination$x)

  plotinfo0 <- msp0$plotinfo
  result4 <- is.null(plotinfo0$lwd.groups)
  result5 <- all(length(plotinfo0$lty.groups) == 1,
                 length(plotinfo0$bg.groups) == 1,
                 length(plotinfo0$col.groups) == 1)

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(groups = species, lwd = 2)

  result6 <- all(c("gr_scores","gr_class") %in% names(msp1$projected))
  result7 <- all(msp1$gr_class == species)
  result8 <- all(msp1$gr_scores == prcomp(geomorph::two.d.array(shapes))$x)

  plotinfo1 <- msp1$plotinfo
  result9 <- all(length(plotinfo1$lwd.groups) == 1,
                  length(plotinfo1$lty.groups) == nlevels(species),
                  length(plotinfo1$bg.groups) == nlevels(species),
                  length(plotinfo1$col.groups) == nlevels(species))


  msp2 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(groups = species, ellipse = TRUE)

  result10 <- all(c("gr_scores","gr_class") %in% names(msp2$projected))
  result11 <- all(msp2$projected$gr_class == species)
  result12 <- all(msp2$gr_scores == prcomp(geomorph::two.d.array(shapes))$x)

  plotinfo2 <- msp2$plotinfo
  result13 <- is.null(plotinfo2$lwd.groups)
  result14 <- all(length(plotinfo2$lty.groups) == nlevels(species),
                  length(plotinfo2$bg.groups) == nlevels(species),
                  length(plotinfo2$col.groups) == nlevels(species))

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,
                  result8,result9,result10,result12,result13,result14))
})


test_that(desc = "testing proj_groups, stacking behavior", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  index <- tails$data$type == "DF"

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(shapes = shapes[,,!index], groups = factor(species[!index]),
                col = "red", lty = 1) %>%
    proj_groups(shapes = shapes[,,index], groups = factor(species[index]),
                col = "blue", lty = 2)

  dec <- 5
  index_x_in_sc <- as.numeric(unlist
                              (apply(round(msp1$projected$gr_scores,dec), 1, \(x, y) {
                                which(apply(y, 1, \(z, x){
                                  all(z == x)}, x))},
                                round(prcomp(geomorph::two.d.array(shapes))$x,dec))))

  result1 <- all(index_x_in_sc == c(which(!index), which(index)))
  result2 <- all(as.character(msp1$projected$gr_class) == as.character(c(species[!index], species[index])))

  plotinfo1 <- msp1$plotinfo
  result3 <- all(plotinfo1$col.groups == c(rep(col2hex("red"), nlevels(factor(species[!index]))),
                                           rep(col2hex("blue"), nlevels(factor(species[index])))))
  result4 <- all(plotinfo1$lty.groups == c(rep(1, nlevels(factor(species[!index]))),
                                           rep(2, nlevels(factor(species[index])))))

  msp2 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(shapes = shapes[,,!index], groups = factor(species[!index]),
                col = 1:11, lty = rep(1, 11)) %>%
    proj_groups(shapes = shapes[,,index], groups = factor(species[index]),
                col = 12:13, lty = rep(2, 2))

  plotinfo2 <- msp2$plotinfo
  result5 <- all(plotinfo2$col.groups == col2hex(1:(nlevels(factor(species[!index])) +
                                                      nlevels(factor(species[index])))))
  result6 <- all(plotinfo2$lty.groups == c(rep(1, nlevels(factor(species[!index]))),
                                           rep(2, nlevels(factor(species[index])))))


  msp3 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_groups(shapes = shapes[,,!index], groups = species[!index], lty = 1) %>%
    proj_groups(shapes = shapes[,,index], groups = species[index], lty = 2)

  index_x_in_sc <- as.numeric(unlist
                              (apply(round(msp3$projected$gr_scores,dec), 1, \(x, y) {
                                which(apply(y, 1, \(z, x){
                                  all(z == x)}, x))},
                                round(prcomp(geomorph::two.d.array(shapes))$x,dec))))

  result7 <- all(index_x_in_sc == c(which(!index), which(index)))
  result8 <- all(as.character(msp3$projected$gr_class) == c(as.character(paste0(species[!index], "_bis.")),
                                                            as.character(species[index])))

  plotinfo3 <- msp3$plotinfo
  cols1 <- cols2 <- col2hex(1:13)
  cols1[!levels(species) %in% as.character(unique(species[!index]))] <-
    cols2[!levels(species) %in% as.character(unique(species[index]))] <- "#FFFFFF"
  result9 <- all(plotinfo3$col.groups == c(cols1,cols2))
  result10 <- all(plotinfo3$lty.groups == c(rep(1, nlevels(unique(species[!index]))),
                                            rep(2, nlevels(unique(species[index])))))

  expect_true(all(result1,result2,result3,result4,result5,result6,
                  result7,result8,result9,result10))
})

###########################################################

test_that(desc = "testing proj_phylogeny, general behavior", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  sp_shapes <- expected_shapes(shapes, species)
  tree <- tails$tree

  msp1 <- mspace(sp_shapes, axes = c(1,2), plot = FALSE) %>%
    proj_phylogeny(sp_shapes, tree = tree, col.tips = seq_len(nlevels(species)),
                   col.nodes = c("red", rep(1, tree$Nnode)))

  result1 <- all(c("phylo_scores","phylo") %in% names(msp1$projected))

  tipscores <- prcomp(geomorph::two.d.array(sp_shapes))$x[tree$tip.label,]
  nodescores <- apply(tipscores,2,phytools::fastAnc, tree=tree)
  sc1 <- round(msp1$projected$phylo_scores,10)
  sc2 <- round(rbind(tipscores,nodescores),10)
  result2 <- all(round(sc1,10) == round(sc2,10))

  sc3 <- proj_phylogeny(mspace = msp1, shapes = sp_shapes, tree = tree, pipe = FALSE)[rownames(msp1$projected$phylo_scores),]
  result3 <- all(round(sc3,10) == round(sc2,10))

  result4 <- all(msp1$phylo$tip.labels == tree$tip.labels)

  result5 <- all(length(msp1$plotinfo$col.tips) == nlevels(species),
                 length(msp1$plotinfo$col.nodes) == length(c("red", rep(1, tree$Nnode))))

  msp2 <- mspace(sp_shapes, axes = c(1,2), plot = FALSE) %>%
    proj_phylogeny(sp_shapes, tree = tree)

  result6 <- all(length(msp2$plotinfo$col.tips) == 1,
                 length(msp2$plotinfo$col.nodes) == 1)

  expect_true(all(result1,result2,result3,result4,result5,result6))
})


###########################################################

test_that(desc = "testing proj_axis, general behavior", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  sp_shapes <- expected_shapes(shapes, species)
  pca <- prcomp(geomorph::two.d.array(shapes))

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_axis(obj = pca)

  result1 <- all("shape_axis" %in% names(msp1$projected))

  sc1 <- proj_axis(obj = pca, axis = 1, mspace = msp1, pipe = FALSE)
  sc2 <- range(msp1$ordination$x[,1])
  result2 <- all(round(sc1[,1],10) == round(sc2,10))

  expect_true(all(result1,result2))
})


test_that(desc = "testing proj_axis, stacking behavior", code = {
  data("tails")
  shapes <- tails$shapes
  species <- tails$data$species
  sp_shapes <- expected_shapes(shapes, species)
  pca <- prcomp(geomorph::two.d.array(shapes))

  msp1 <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
    proj_axis(obj = pca, axis = 1, col = "red") %>%
    proj_axis(obj = pca, axis = 2, col = "blue")

  result1 <- length(msp1$projected$shape_axis) == 2
  result2 <- length(msp1$plotinfo$col.axis) == length(msp1$projected$shape_axis)

  expect_true(all(result1,result2))
})

###########################################################

test_that(desc = "testing proj_landscape, general behavior", code = {
  data("tails")
  shapes <- tails$shapes
  type <- tails$data$type

  msp1 <- mspace(shapes, axes = c(1,2)) %>%
    proj_landscape(FUN = morphospace:::computeLD)

  result1 <- "landsc" %in% names(msp1$projected)
  result2 <- all(c("x","z") %in% names(msp1$projected$landsc))
  result3 <- msp1$projected$landsc$type == "theoretical"
  result4 <- all(length(msp1$plotinfo$lty.landsc) == 1,
                 length(msp1$plotinfo$lwd.landsc) == 1,
                 length(msp1$plotinfo$alpha.landsc) == 1,
                 length(msp1$plotinfo$display.landsc) == 1,
                 length(msp1$plotinfo$resolution.landsc) == 1,
                 length(msp1$plotinfo$nlevels.landsc) == 1,
                 length(msp1$plotinfo$nlevels.landsc) == 1)


  msp2 <- mspace(shapes, axes = 1) %>%
    proj_landscape(FUN = morphospace:::computeLD)

  result5 <- "landsc" %in% names(msp2$projected)
  result6 <- all(c("x","z") %in% names(msp2$projected$landsc))
  result7 <- msp2$projected$landsc$type == "theoretical"
  result8 <- all(length(msp2$plotinfo$lty.landsc) == 1,
                 length(msp2$plotinfo$lwd.landsc) == 1,
                 length(msp2$plotinfo$alpha.landsc) == 1,
                 length(msp2$plotinfo$display.landsc) == 1,
                 length(msp2$plotinfo$resolution.landsc) == 1,
                 length(msp2$plotinfo$nlevels.landsc) == 1,
                 length(msp2$plotinfo$nlevels.landsc) == 1)


  msp3 <- mspace(shapes, axes = c(1,2)) %>%
    proj_landscape(shapes[,,type == "NDF"], FUN = morphospace:::computeLD)

  result9 <- "landsc" %in% names(msp3$projected)
  result10 <- all(c("x","z") %in% names(msp3$projected$landsc))
  result11 <- msp3$projected$landsc$type == "empirical"
  result12 <- all(length(msp3$plotinfo$lty.landsc) == 1,
                  length(msp3$plotinfo$lwd.landsc) == 1,
                  length(msp3$plotinfo$alpha.landsc) == 1,
                  length(msp3$plotinfo$display.landsc) == 1,
                  length(msp3$plotinfo$resolution.landsc) == 1,
                  length(msp3$plotinfo$nlevels.landsc) == 1,
                  length(msp3$plotinfo$nlevels.landsc) == 1)


  msp4 <- mspace(shapes, axes = 1) %>%
    proj_landscape(shapes[,,type == "NDF"], FUN = morphospace:::computeLD)

  result13 <- "landsc" %in% names(msp4$projected)
  result14 <- all(c("x","z") %in% names(msp4$projected$landsc))
  result15 <- msp4$projected$landsc$type == "empirical"
  result16 <- all(length(msp4$plotinfo$lty.landsc) == 1,
                  length(msp4$plotinfo$lwd.landsc) == 1,
                  length(msp4$plotinfo$alpha.landsc) == 1,
                  length(msp4$plotinfo$display.landsc) == 1,
                  length(msp4$plotinfo$resolution.landsc) == 1,
                  length(msp4$plotinfo$nlevels.landsc) == 1,
                  length(msp4$plotinfo$nlevels.landsc) == 1)

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,
                  result8,result9,result10,result12,result13,result14,
                  result15,result16))

  dev.off()
})

