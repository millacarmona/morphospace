test_that(desc = "testing expected_shapes, factors, ndims & accuracy, for landmark data", code = {
  data(tails)

  mshape_all <- expected_shapes(tails$shapes)
  ndims <- dim(mshape_all)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == nrow(tails$shapes)
  result3 <- ndims[2] == ncol(tails$shapes)

  testshape1 <- matrix(colMeans(geomorph::two.d.array(tails$shapes)),ncol=2,byrow=TRUE)
  result4 <- all(round(mshape_all,10) == round(testshape1,10))


  mshape_spp <- expected_shapes(tails$shape, tails$data$species)
  ndims <- dim(mshape_spp)

  result5 <- length(ndims) == 3
  result6 <- ndims[1] == nrow(tails$shapes)
  result7 <- ndims[2] == ncol(tails$shapes)
  result8 <- ndims[3] == nlevels(tails$data$species)

  testshapes2 <- geomorph::arrayspecs(
    apply(geomorph::two.d.array(tails$shapes), 2, tapply, tails$data$species, mean),
    k = 2, p = 9)
  result9 <- all(round(mshape_spp,10) == round(testshapes2,10))


  mshape_ndf <- expected_shapes(tails$shapes, tails$data$type, "NDF")
  ndims <- dim(mshape_ndf)

  result10 <- length(ndims) == 2
  result11 <- ndims[1] == nrow(tails$shapes)
  result12 <- ndims[2] == ncol(tails$shapes)

  testshape3 <- matrix(colMeans(
    geomorph::two.d.array(tails$shapes[,,tails$data$type == "NDF"])),ncol=2,byrow=TRUE)
  result13 <- all(round(mshape_ndf,10) == round(testshape3,10))


  expect_true(all(result1, result2, result3, result4, result5,result6,
                  result7, result8, result9, result10, result11, result12, result13))

})

test_that(desc = "testing expected_shapes, numerics, ndims & accuracy, for landmark data", code = {
  data(tails)

  shapes <- tails$shapes
  sizes <- tails$sizes
  predshapes <- expected_shapes(shapes, x = sizes)
  ndims <- dim(predshapes)

  result1 <- length(ndims) == 3
  result2 <- ndims[1] == nrow(tails$shapes)
  result3 <- ndims[2] == ncol(tails$shapes)
  result4 <- ndims[3] == length(sizes)

  testshapes1 <- geomorph::arrayspecs(
    lm(geomorph::two.d.array(shapes) ~ sizes)$fitted,
    k = 2, p = 9)
  result5 <- all(round(predshapes,10) == round(testshapes1,10))


  predshape_large <- expected_shapes(shapes, x = sizes, xvalue = max(sizes))
  ndims <- dim(predshape_large)

  result6 <- length(ndims) == 2
  result7 <- ndims[1] == nrow(tails$shapes)
  result8 <- ndims[2] == ncol(tails$shapes)

  testshape2 <- geomorph::arrayspecs(
    lm(geomorph::two.d.array(shapes) ~ sizes)$fitted,
    k = 2, p = 9)[,,which.max(sizes)]
  result9 <- all(round(predshape_large,10) == round(testshape2,10))

  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8, result9))


})


test_that(desc = "testing expected_shapes, factors, ndims & accuracy, for landmark data, with phylogeny", code = {
  data(tails)
  spp_shapes <- expected_shapes(tails$shape, tails$data$species)
  spp_type <- factor(c(tapply(as.character(tails$data$type), tails$data$species, unique)))
  spp_sizes <- c(tapply(tails$sizes, tails$data$species, mean))

  mshape_all <- expected_shapes(shapes = spp_shapes, tree = tails$tree)
  ndims <- dim(mshape_all)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == nrow(tails$shapes)
  result3 <- ndims[2] == ncol(tails$shapes)

  testshape1 <- matrix(
    apply(geomorph::two.d.array(spp_shapes),2,phytools::fastAnc,tree=tails$tree)[1,],
    ncol=2,byrow=TRUE)
  result4 <- all(round(mshape_all,10) == round(testshape1,10))


  mshape_type <- expected_shapes(spp_shapes, spp_type, tree = tails$tree)
  ndims <- dim(mshape_type)

  result5 <- length(ndims) == 3
  result6 <- ndims[1] == nrow(tails$shapes)
  result7 <- ndims[2] == ncol(tails$shapes)
  result8 <- ndims[3] == nlevels(tails$data$type)


  testshapes2 <- geomorph::arrayspecs(
    unique(
      do.call("cbind",
              lapply(seq_len(ncol(geomorph::two.d.array(spp_shapes))), function(i) {
                cd=caper::comparative.data(
                  data = data.frame(shapes = geomorph::two.d.array(spp_shapes)[tails$tree$tip.label,i],
                                    type = spp_type[tails$tree$tip.label],
                                    n = tails$tree$tip.label),
                  names.col = "n", phy = tails$tree)
                pgls <- caper::pgls(shapes ~ type, data = cd)
                pgls$fitted
              }))), k = 2, p = 9)
  result9 <- all(round(testshapes2,10) == round(mshape_type[,,2:1],10))


  mshape_ndf <- expected_shapes(spp_shapes, spp_type, "NDF", tree = tails$tree)
  ndims <- dim(mshape_ndf)

  result10 <- length(ndims) == 2
  result11 <- ndims[1] == nrow(tails$shapes)
  result12 <- ndims[2] == ncol(tails$shapes)

  testshape3 <- matrix(unlist(lapply(seq_len(ncol(geomorph::two.d.array(spp_shapes))), function(i) {
    cd=caper::comparative.data(
      data = data.frame(shapes = geomorph::two.d.array(spp_shapes)[tails$tree$tip.label,i],
                        type = spp_type[tails$tree$tip.label],
                        n = tails$tree$tip.label),
      names.col = "n", phy = tails$tree)
    pgls <- caper::pgls(shapes ~ type, data = cd)
    sum(pgls[[1]]$coef)
  })), ncol = 2, byrow = TRUE)
  result13 <- all(round(testshape3,10) == round(mshape_ndf,10))

  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8, result9, result10,
                  result11, result12, result13))

})


test_that(desc = "testing expected_shapes, numerics, ndims & accuracy, for landmark data, with phylogeny", code = {
  data(tails)

  spp_shapes <- expected_shapes(tails$shapes, tails$data$species)
  spp_sizes <- c(tapply(tails$sizes, tails$data$species, mean))
  predshapes <- expected_shapes(spp_shapes, x = spp_sizes, tree = tails$tree)
  ndims <- dim(predshapes)

  result1 <- length(ndims) == 3
  result2 <- ndims[1] == nrow(tails$shapes)
  result3 <- ndims[2] == ncol(tails$shapes)
  result4 <- ndims[3] == length(spp_sizes)

  testshapes1 <- geomorph::arrayspecs(
    unique(
      do.call("cbind",
              lapply(seq_len(ncol(geomorph::two.d.array(spp_shapes))), function(i) {
                cd=caper::comparative.data(
                  data = data.frame(shapes = geomorph::two.d.array(spp_shapes)[tails$tree$tip.label,i],
                                    sizes = spp_sizes[tails$tree$tip.label],
                                    n = tails$tree$tip.label),
                  names.col = "n", phy = tails$tree)
                pgls <- caper::pgls(shapes ~ sizes, data = cd)
                pgls$fitted
              }))), k = 2, p = 9)
  result5 <- all(round(testshapes1,10) == round(predshapes[,,tails$tree$tip.label],10))


  predshape_large <- expected_shapes(spp_shapes, x = spp_sizes,
                                     xvalue = max(spp_sizes), tree = tails$tree)
  ndims <- dim(predshape_large)

  result6 <- length(ndims) == 2
  result7 <- ndims[1] == nrow(tails$shapes)
  result8 <- ndims[2] == ncol(tails$shapes)

  testshapes2 <- geomorph::arrayspecs(
    unique(
      do.call("cbind",
              lapply(seq_len(ncol(geomorph::two.d.array(spp_shapes))), function(i) {
                cd=caper::comparative.data(
                  data = data.frame(shapes = geomorph::two.d.array(spp_shapes)[tails$tree$tip.label,i],
                                    sizes = spp_sizes[tails$tree$tip.label],
                                    n = tails$tree$tip.label),
                  names.col = "n", phy = tails$tree)
                pgls <- caper::pgls(shapes ~ sizes, data = cd)
                pgls$fitted
              }))), k = 2, p = 9)[,,which.max(spp_sizes[tails$tree$tip.label])]

  result9 <- all(round(testshapes2,10) == round(predshape_large,10))

  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8, result9))


})



test_that(desc = "testing expected_shapes, factors, ndims & accuracy, for Fourier data", code = {
  data(shells)

  mshape_all <- expected_shapes(shells$shapes)
  ndims <- dim(mshape_all)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == 1
  result3 <- ndims[2] == ncol(shells$shapes)

  testshape1 <- colMeans(shells$shapes$coe)
  result4 <- all(round(mshape_all,10) == round(testshape1,10))


  mshape_spp <- expected_shapes(shells$shape, shells$data$species)
  ndims <- dim(mshape_spp)

  result5 <- length(ndims) == 2
  result6 <- ndims[1] == nlevels(shells$data$species)
  result7 <- ndims[2] == ncol(shells$shapes)

  testshapes2 <- apply(shells$shapes$coe, 2, tapply, shells$data$species, mean)
  result8 <- all(round(mshape_spp,10) == round(testshapes2,10))


  mshape_mula <- expected_shapes(shells$shapes, x = shells$data$locality, xvalue = "Agua de la Mula")
  ndims <- dim(mshape_mula)

  result9 <- length(ndims) == 2
  result10 <- ndims[1] == 1
  result11 <- ndims[2] == ncol(shells$shapes)

  testshape3 <- colMeans(shells$shapes$coe[shells$data$locality == "Agua de la Mula",])
  result12 <- all(round(mshape_mula,10) == round(testshape3,10))


  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8, result9, result10,
                  result11, result12))


})


test_that(desc = "testing expected_shapes, numerics, ndims & accuracy, for Fourier data", code = {
  data(shells)

  shapes <- shells$shapes$coe
  sizes <- shells$sizes
  predshapes <- expected_shapes(shapes, x = sizes)
  ndims <- dim(predshapes)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == nrow(shells$shapes)
  result3 <- ndims[2] == ncol(shells$shapes)

  testshapes1 <- lm(shapes ~ sizes)$fitted
  result4 <- all(round(predshapes,10) == round(testshapes1,10))


  predshape_large <- expected_shapes(shapes, x = sizes, xvalue = max(sizes))
  ndims <- dim(predshape_large)

  result5 <- length(ndims) == 2
  result6 <- ndims[1] == 1
  result7 <- ndims[2] == ncol(shells$shapes)

  testshape2 <- lm(shapes ~ sizes)$fitted[which.max(sizes),]
  result8 <- all(round(predshape_large,10) == round(testshape2,10))

  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8))


})

#further tests covering accuracy of the actual shape(s) are implemented in tests_internal_builders


###########################################

test_that(desc = "testing detrend_shapes, method = orthogonal + default settings,
          for factors and numerics", code = {
  data(shells)

  shapes <- shells$shapes$coe
  species <- shells$data$species
  logsizes <- log(shells$sizes)

  model1 <- lm(shapes ~ logsizes)
  detr_shapes1 <- detrend_shapes(model1, method = "orthogonal")
  result1 <- all(dim(shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

  pca <- prcomp(detr_shapes1)
  testshape1 <- rev_eigen(colMeans(pca$x), pca$rotation, pca$center)
  detr_mshape1 <- expected_shapes(detr_shapes1)
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, method = "orthogonal")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- rev_eigen(colMeans(pca$x), pca$rotation, pca$center)
  detr_mshape2 <- expected_shapes(detr_shapes2)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})

test_that(desc = "testing detrend_shapes, method = residuals + default settings,
          for factors and numerics", code = {
  data(shells)

  shapes <- shells$shapes$coe
  species <- shells$data$species
  logsizes <- log(shells$sizes)

  model1 <- lm(shapes ~ logsizes)
  detr_shapes1 <- detrend_shapes(model1, method = "residuals")
  result1 <- all(dim(shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

  pca <- prcomp(detr_shapes1)
  testshape1 <- rev_eigen(colMeans(pca$x), pca$rotation, pca$center)
  detr_mshape1 <- expected_shapes(detr_shapes1)
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, method = "residuals")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- rev_eigen(colMeans(pca$x), pca$rotation, pca$center)
  detr_mshape2 <- expected_shapes(detr_shapes2)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})


test_that(desc = "testing detrend_shapes, method = orthogonal, xvalue, for factors and numerics", code = {
  data(shells)

  shapes <- shells$shapes$coe
  species <- shells$data$species
  logsizes <- log(shells$sizes)

  model1 <- lm(shapes ~ logsizes)
  detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "orthogonal")
  result1 <- all(dim(shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

  pca <- prcomp(detr_shapes1)
  testshape1 <- rev_eigen(lm(pca$x ~ logsizes)$fitted[which.max(logsizes),], pca$rotation, pca$center)
  detr_mshape1 <- expected_shapes(detr_shapes1)
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "orthogonal")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- rev_eigen(colMeans(pca$x[species == "koeneni",]), pca$rotation, pca$center)
  detr_mshape2 <- expected_shapes(detr_shapes2)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})

test_that(desc = "testing detrend_shapes, method = residuals, xvalue, for factors and numerics", code = {
  data(shells)

  shapes <- shells$shapes$coe
  species <- shells$data$species
  logsizes <- log(shells$sizes)

  model1 <- lm(shapes ~ logsizes)
  detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "residuals")
  result1 <- all(dim(shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

  pca <- prcomp(detr_shapes1)
  testshape1 <- rev_eigen(lm(pca$x ~ logsizes)$fitted[which.max(logsizes),], pca$rotation, pca$center)
  detr_mshape1 <- expected_shapes(detr_shapes1)
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "residuals")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- rev_eigen(colMeans(pca$x[species == "koeneni",]), pca$rotation, pca$center)
  detr_mshape2 <- expected_shapes(detr_shapes2)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})


test_that(desc = "testing detrend_shapes, method = orthogonal, newdata, for numerics", code = {

  data('shells3D')
  shapes <- shells3D$shapes
  species <- shells3D$data$species
  sizes <- log(shells3D$sizes)

  index1 <- species == levels(species)[7]
  shapes1 <- shapes[,,index1]
  logsizes1 <- sizes[index1]

  index2 <- species == levels(species)[3]
  shapes2 <- shapes[,,index2]
  logsizes2 <- sizes[index2]

  mod1 <- lm(geomorph::two.d.array(shapes1) ~ logsizes1)
  mod2 <- lm(geomorph::two.d.array(shapes2) ~ logsizes2)


  general_space <- prcomp(geomorph::two.d.array(shapes))

  burn1 <- burnaby(x = geomorph::two.d.array(shapes1), vars = logsizes1)
  burn2 <- burnaby(x = geomorph::two.d.array(shapes2), vars = logsizes2)

  detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(sizes[index2]),
                                     newdata = mod2)

  ax <- c(1:30)
  na_shapes2minus1 <- rev_eigen(scores = burn2$x[,ax],
                                vectors = burn2$rotation[,ax] - burn1$rotation[,ax],
                                center = colMeans(detshapes2using1))

  scores2using1 <- proj_eigen(detshapes2using1,
                              vectors = general_space$rotation[,1:2],
                              center = general_space$center)
  slope2using1 <- lm(scores2using1[,2] ~ scores2using1[,1])$coef[2]

  scores2minus1 <- proj_eigen(na_shapes2minus1,
                              vectors = general_space$rotation[,1:2],
                              center = general_space$center)
  slope2minus1 <- lm(scores2minus1[,2] ~ scores2minus1[,1])$coef[2]


  result1 <- round(slope2using1, 3) == round(slope2minus1, 3)
  expect_true(all(result1))

  # refmesh <- shells3D$mesh_meanspec
  # template <- Morpho::tps3d(x = refmesh, refmat = shapes[,,geomorph::findMeanSpec(shapes)], tarmat = expected_shapes(shapes))
  # msp <- mspace(shapes, template = template, bg.models = "gray",
  #               cex.ldm = 0, alpha.models = 0.7, adj_frame = c(0.93,0.93)) %>%
  #   proj_shapes(detshapes2using1, pch = 1, col = 4) %>%
  #   proj_shapes(na_shapes2minus1, pch = 21, bg = 4, cex= 0.5)
})

test_that(desc = "testing detrend_shapes, method = residuals, newdata, for numerics", code = {
  data(shells)

  index1 <- shells$data$species == "koeneni"
  shapes1 <- shells$shapes$coe[index1,]
  logsizes1 <- log(shells$sizes[index1])

  index2 <- shells$data$species == "windhauseni"
  shapes2 <- shells$shapes$coe[index2,]
  logsizes2 <- log(shells$sizes[index2])

  model_for2 <- lm(shapes2 ~ logsizes2)
  model_for1 <- lm(shapes1 ~ logsizes1)
  detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

  slope1<-lm(shapes1~logsizes1)$coefficients[2,]
  slope2<-lm(shapes2~logsizes2)$coefficients[2,]

  slope1minus2 <- lm(shapes1~logsizes1)$coefficients[2,] - lm(shapes2~logsizes2)$coefficients[2,]
  slopedetr1 <- lm(detr_shapes1~logsizes1)$coefficients[2,]

  results <- all(round(slope1minus2,10) == round(slopedetr1,10))
  expect_true(results)

})


test_that(desc = "testing detrend_shapes, method = residuals,
          with phylogeny, for factors and numerics", code = {
  data(tails)

  shapes <- tails$shapes
  species <- tails$data$species
  logsizes <- log(tails$sizes)

  tree <- tails$tree
  sp_shapes <- geomorph::two.d.array(expected_shapes(shapes, species))
  sp_logsizes <- c(tapply(logsizes,species,mean))

  model1 <- lm(sp_shapes ~ sp_logsizes)
  detr_shapes1 <- detrend_shapes(model1, tree = tree, method = "residuals")
  result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  extres <- phytools::phyl.resid(tree = tree, x = sp_logsizes,
                                 Y = sp_shapes)$resid[rownames(sp_shapes),]

  extmean <- apply(sp_shapes, 2, phytools::fastAnc, tree= tree)[1,]
  intres <- t(t(detr_shapes1) - extmean)
  result3 <- all(round(extres,10) == round(intres,10))

  pca <- prcomp(detr_shapes1)
  testshape1 <- matrix(rev_eigen(colMeans(pca$x), pca$rotation, pca$center), ncol = 2, byrow = TRUE)
  detr_mshape1 <- expected_shapes(geomorph::arrayspecs(detr_shapes1, k = 2, p = 9))
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))

  sp_type <- rep("NDF", 13) ; sp_type[c(7,10)] <- "DF"
  sp_type <- factor(sp_type) ; sp_type <- setNames(sp_type, levels(species))
  model2 <- lm(sp_shapes ~ sp_type)
  detr_shapes2 <- detrend_shapes(model2, tree = tree, method = "residuals")
  result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- matrix(rev_eigen(colMeans(pca$x), pca$rotation, pca$center), ncol = 2, byrow = TRUE)
  detr_mshape2 <- expected_shapes(geomorph::arrayspecs(detr_shapes2, k = 2, p = 9))
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  extres <- phytools::phyl.resid(tree = tree, x = sp_type,
                                 Y = sp_shapes)$resid[rownames(sp_shapes),]
  extmean <- apply(sp_shapes, 2, phytools::fastAnc, tree= tree)[1,]
  intres <- t(t(detr_shapes2) - extmean)
  result8 <- all(round(extres,10) == round(intres,10))

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})

test_that(desc = "testing detrend_shapes, method = orthogonal,
          with phylogeny, for factors and numerics", code = {

  data(tails)

  shapes <- tails$shapes
  species <- tails$data$species
  logsizes <- log(tails$sizes)

  tree <- tails$tree
  sp_shapes <- geomorph::two.d.array(expected_shapes(shapes, species))
  sp_logsizes <- c(tapply(logsizes,species,mean))

  model1 <- lm(sp_shapes ~ sp_logsizes)
  detr_shapes1 <- detrend_shapes(model1, tree = tree, method = "orthogonal")
  result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  exp_shapes1 <- expected_shapes(shapes = sp_shapes, x = sp_logsizes,
                                 tree = tree, returnarray = FALSE)
  angle <- Morpho::angle.calc(prcomp(detr_shapes1)$rotation[,1],
                              prcomp(exp_shapes1)$rotation[,1]) * 57.2958
  result3 <- all(round(angle,3) == 90)

  pca <- prcomp(detr_shapes1)
  testshape1 <- matrix(rev_eigen(colMeans(pca$x), pca$rotation, pca$center), ncol = 2, byrow = TRUE)
  detr_mshape1 <- expected_shapes(geomorph::arrayspecs(detr_shapes1, k = 2, p = 9))
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))

  sp_type <- rep("NDF", 13) ; sp_type[c(7,10)] <- "DF"
  sp_type <- factor(sp_type) ; sp_type <- setNames(sp_type, levels(species))
  model2 <- lm(sp_shapes ~ sp_type)
  detr_shapes2 <- detrend_shapes(model2, tree = tree, method = "orthogonal")
  result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- matrix(rev_eigen(colMeans(pca$x), pca$rotation, pca$center), ncol = 2, byrow = TRUE)
  detr_mshape2 <- expected_shapes(geomorph::arrayspecs(detr_shapes2, k = 2, p = 9))
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  exp_shapes2 <- expected_shapes(shapes = sp_shapes, x = sp_type,
                                 tree = tree, returnarray = FALSE)
  angle <- Morpho::angle.calc(prcomp(detr_shapes2)$rotation[,1],
                              prcomp(exp_shapes2)$rotation[,1]) * 57.2958
  result8 <- all(round(angle,3) == 90)


  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})

############################

test_that(desc = "testing rotate_fcoef", code = {
  data(shells)

  shapes_coe <- rbind(shells$shapes[1:10,])
  shapes_rot_coe <- rotate_fcoef(shapes_coe)

  result1 <- all(dim(shapes_coe) == dim(shapes_rot_coe))

  shape1_coe <- rbind(shells$shapes[1,])
  shape1_rot_coe <- rotate_fcoef(shape1_coe)

  result2 <- all(dim(shape1_coe) == dim(shape1_rot_coe))

  shape1_coo <- inv_efourier(shape1_coe, 300)
  shape1_rot_coo <- inv_efourier(shape1_rot_coe, 300)

  testshape <- geomorph::rotate.coords(geomorph::rotate.coords(shape1_coo, type = "rotateCC"), type = "rotateCC")
  result3 <- all(round(testshape[order(testshape[,1])],5) == round(shape1_rot_coo[order(shape1_rot_coo[,1])],5))

  expect_true(result1, result2, result3)
})

#############################


test_that(desc = "testing correct_efourier, automatic behavior", code = {
  data(shells)

  shapes_coe <- shells$shapes
  shapes_rot_coe <- correct_efourier(shapes_coe, index = 1:10)

  result1 <- all(dim(shapes_coe) == dim(shapes_rot_coe))

  shape1_coo <- inv_efourier(shapes_coe$coe[1,], 300)
  shape1_rot_coo <- inv_efourier(shapes_rot_coe$coe[1,], 300)

  testshape <- geomorph::rotate.coords(geomorph::rotate.coords(shape1_coo, type = "rotateCC"), type = "rotateCC")
  result2 <- all(round(testshape[order(testshape[,1])],5) == round(shape1_rot_coo[order(shape1_rot_coo[,1])],5))

  expect_true(result1, result2)
  dev.off()
})


###############################


test_that(desc = "testing ax_transformations, dimensions", code = {
  data(tails)

  shapes <- tails$shapes
  sizes <- tails$sizes
  species <- tails$data$species
  tree <- tails$tree

  model <- lm(geomorph::two.d.array(shapes) ~ sizes)
  pca <- prcomp(geomorph::two.d.array(shapes))
  bgpca <- bg_prcomp(geomorph::two.d.array(shapes), groups = species)
  phypca <- phy_prcomp(geomorph::two.d.array(expected_shapes(shapes, species)), tree = tree)
  pls <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = sizes)

  ext1 <- ax_transformation(obj = model)
  ext2 <- ax_transformation(obj = pca)
  ext3 <- ax_transformation(obj = bgpca)
  ext4 <- ax_transformation(obj = phypca)
  ext5 <- ax_transformation(obj = pls)

  results1 <- all(c(dim(ext1)[1], dim(ext2)[1], dim(ext3)[1], dim(ext4)[1], dim(ext4)[1]) == 2)
  results2 <- all(c(dim(ext1)[2], dim(ext2)[2], dim(ext3)[2], dim(ext4)[2], dim(ext4)[2]) == nrow(shapes) * ncol(shapes))

  expect_true(all(results1, results2))
})



test_that(desc = "testing ax_transformations, accuracy of shapes", code = {
  data("tails")

  shapes <- tails$shapes
  sizes <- tails$sizes
  sex <- tails$data$sex

  model1 <- lm(geomorph::two.d.array(shapes) ~ sizes)
  ext_model <- geomorph::arrayspecs(ax_transformation(obj = model1, mag = 1), k = ncol(shapes), p = nrow(shapes))

  pls <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = sizes)
  ext_pls <- geomorph::arrayspecs(ax_transformation(obj = pls, mag = 1), k = ncol(shapes), p = nrow(shapes))
  ext_pls2 <- geomorph::arrayspecs(rev_eigen(range(lm(pls$x ~ sizes)$fitted), pls$rotation[,1], pls$center),
                                   k = ncol(shapes), p = nrow(shapes))

  result1 <- all(round(ext_model[,,order(ext_model[1,1,])],10) == round(ext_pls2[,,order(ext_pls2[1,1,])],10))

  model2 <- lm(geomorph::two.d.array(shapes) ~ sex)
  ext_model <- geomorph::arrayspecs(ax_transformation(obj = model2, mag = 1), k = ncol(shapes), p = nrow(shapes))

  bgpca <- bg_prcomp(geomorph::two.d.array(shapes), groups = sex)
  ext_bgpca <- geomorph::arrayspecs(ax_transformation(obj = bgpca, mag = 1), k = ncol(shapes), p = nrow(shapes))
  ext_bgpca2 <- geomorph::arrayspecs(rev_eigen(range(lm(bgpca$x[,1] ~ sex)$fitted), bgpca$rotation[,1], bgpca$center),
                                     k = ncol(shapes), p = nrow(shapes))

  result2 <- all(round(ext_model[,,order(ext_model[1,1,])],10) == round(ext_bgpca2[,,order(ext_bgpca2[1,1,])],10))

  expect_true(all(result1, result2))
})

######################################################################################

test_that(desc = "testing extract_shapes, landmark data, non-interactive behavior and accuracy", code = {
  data(wings)

  shapes <- wings$shapes
  temp <- wings$template
  sizes <- wings$sizes
  sex <- wings$data$sex

  msp1 <- mspace(shapes, nh = 8, nv = 6, plot = FALSE)
  es1 <- suppressWarnings(extract_shapes(msp1))

  result1 <- names(es1) == "shapes"
  result2 <- all(dim(es1$shapes) == dim(msp1$projected$shapemodels))
  result3 <- all(round(es1$shapes, 10) == round(msp1$projected$shapemodels, 10))


  nsh <- 5
  ax <- 2
  es2 <- suppressWarnings(extract_shapes(msp1, axis = ax, nshapes = nsh))
  vals <- seq(from = min(msp1$ordination$x[,ax]), to = max(msp1$ordination$x[,ax]),
              length.out = nsh)
  newshps1 <- geomorph::arrayspecs(rev_eigen(scores = vals,
                                             vectors = msp1$ordination$rotation[,ax],
                                             center = msp1$ordination$center), k =2, p = 9)
  result4 <- dim(es2$shapes)[3] == nsh
  result5 <- all(newshps1 == es2$shapes)


  es3 <- suppressWarnings(extract_shapes(msp1, nshapes = 3,
                                         scores = cbind(0,vals)))
  result6 <- dim(es3$shapes)[3] == length(vals)
  result7 <- all(newshps1 == es3$shapes)


  msp2 <- mspace(shapes, template = temp, nh = 8, nv = 6, plot = FALSE)
  es4 <- extract_shapes(msp2)
  result8 <- all(c("templates", "shapes") %in% names(es4))
  result9 <- all(dim(es4$templates)[3] == dim(msp1$projected$shapemodels)[3])
  result10 <- all(dim(es4$templates)[2] + dim(msp1$projected$shapemodels)[2])
  result11 <- all(dim(es4$templates)[1] + dim(msp1$projected$shapemodels)[1] == dim(temp)[1])

  es5 <- extract_shapes(msp2, scores = c(0,0))
  mtemp <- Morpho::tps2d(x = temp[-c(1:9),], refmat = temp[1:9,],
                            tarmat = expected_shapes(shapes))
  result12 <- all(na.omit(round(mtemp,10) == round(es5$templates[,,1],10)))


  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8,
                  result9, result10, result11, result12))
})


test_that(desc = "testing extract_shapes, Fourier data, non-interactive behavior and accuracy", code = {
  data(shells)

  shapes <- shells$shapes$coe
  temp <- shells$template
  sizes <- shells$sizes
  sex <- shells$data$sex

  msp1 <- mspace(shapes, nh = 8, nv = 6, plot = FALSE)
  es1 <- suppressWarnings(extract_shapes(msp1))

  result1 <- names(es1) == "shapes"
  result2 <- all(dim(es1$shapes) == dim(msp1$projected$shapemodels))
  result3 <- all(round(es1$shapes, 10) == round(msp1$projected$shapemodels, 10))


  nsh <- 5
  ax <- 2
  es2 <- suppressWarnings(extract_shapes(msp1, axis = ax, nshapes = nsh))
  vals <- seq(from = min(msp1$ordination$x[,ax]), to = max(msp1$ordination$x[,ax]),
              length.out = nsh)
  newshps1 <- rev_eigen(scores = vals,
                        vectors = msp1$ordination$rotation[,ax],
                        center = msp1$ordination$center)
  result4 <- dim(es2$shapes)[1] == nsh
  result5 <- all(newshps1 == es2$shapes)


  es3 <- suppressWarnings(extract_shapes(msp1, nshapes = 3,
                                         scores = cbind(0,vals)))
  result6 <- dim(es3$shapes)[1] == length(vals)
  result7 <- all(newshps1 == es3$shapes)


  expect_true(all(result1, result2, result3, result4, result5, result6, result7))
})




