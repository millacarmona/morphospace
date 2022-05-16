
##################################

test_that(desc = "testing consensus, dimensions, for landmark data", code = {
  data(tails)

  mshape_all <- consensus(tails$shapes)
  ndims <- dim(mshape_all)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == nrow(tails$shapes)
  result3 <- ndims[2] == ncol(tails$shapes)

  mshape_spp <- consensus(tails$shape, tails$data$species)
  ndims <- dim(mshape_spp)

  result4 <- length(ndims) == 3
  result5 <- ndims[1] == nrow(tails$shapes)
  result6 <- ndims[2] == ncol(tails$shapes)
  result7 <- ndims[3] == nlevels(tails$data$species)

  mshape_ndf <- consensus(tails$shapes, which(tails$data$type == "NDF"))
  ndims <- dim(mshape_ndf)

  result8 <- length(ndims) == 2
  result9 <- ndims[1] == nrow(tails$shapes)
  result10 <- ndims[2] == ncol(tails$shapes)

  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8, result9, result10))


})


test_that(desc = "testing consensus, dimensions, for Fourier data", code = {
  data(shells)

  mshape_all <- consensus(shells$shapes)
  ndims <- dim(mshape_all)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == 1
  result3 <- ndims[2] == ncol(shells$shapes)

  mshape_spp <- consensus(shells$shape, shells$data$species)
  ndims <- dim(mshape_spp)

  result4 <- length(ndims) == 2
  result5 <- ndims[1] == nlevels(shells$data$species)
  result6 <- ndims[2] == ncol(shells$shapes)

  mshape_mula <- consensus(shells$shapes, which(shells$data$locality == "Agua de la Mula"))
  ndims <- dim(mshape_mula)

  result7 <- length(ndims) == 2
  result8 <- ndims[1] == 1
  result9 <- ndims[2] == ncol(shells$shapes)

  expect_true(all(result1, result2, result3, result4, result5,
                  result6, result7, result8, result9))


})


#further tests regarding accuracy of the actual shape(s) are implemented in tests_internal_builders


###########################################

test_that(desc = "testing detrend_shapes, default settings, for factors and numerics", code = {
  data(shells)

  shapes <- shells$shapes$coe
  species <- shells$data$species
  logsizes <- log(shells$sizes)

  model1 <- lm(shapes ~ logsizes)
  detr_shapes1 <- detrend_shapes(model1)
  result1 <- all(dim(shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

  pca <- prcomp(detr_shapes1)
  testshape1 <- rev_eigen(colMeans(pca$x), pca$rotation, pca$center)
  detr_mshape1 <- consensus(detr_shapes1)
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2)
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- rev_eigen(colMeans(pca$x), pca$rotation, pca$center)
  detr_mshape2 <- consensus(detr_shapes2)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})


test_that(desc = "testing detrend_shapes, xvalue, for factors and numerics", code = {
  data(shells)

  shapes <- shells$shapes$coe
  species <- shells$data$species
  logsizes <- log(shells$sizes)

  model1 <- lm(shapes ~ logsizes)
  detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes))
  result1 <- all(dim(shapes) == dim(detr_shapes1))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
  result2 <- (totvar_raw > totvar_detr)

  result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

  pca <- prcomp(detr_shapes1)
  testshape1 <- rev_eigen(lm(pca$x ~ logsizes)$fitted[which.max(logsizes),], pca$rotation, pca$center)
  detr_mshape1 <- consensus(detr_shapes1)
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, "koeneni")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  pca <- prcomp(detr_shapes2)
  testshape2 <- rev_eigen(colMeans(pca$x[species == "koeneni",]), pca$rotation, pca$center)
  detr_mshape2 <- consensus(detr_shapes2)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})



test_that(desc = "testing detrend_shapes, newx/newy, for numerics", code = {
  data(shells)

  index1 <- shells$data$species == "koeneni"
  shapes1 <- shells$shapes$coe[index1,]
  logsizes1 <- log(shells$sizes[index1])

  index2 <- shells$data$species == "windhauseni"
  shapes2 <- shells$shapes$coe[index2,]
  logsizes2 <- log(shells$sizes[index2])

  model_for2 <- lm(shapes2 ~ logsizes2)
  detr_shapes1 <- detrend_shapes(model_for2, newx = logsizes1, newy = shapes1)

  slope1<-lm(shapes1~logsizes1)$coefficients[2,]
  slope2<-lm(shapes2~logsizes2)$coefficients[2,]

  slope1minus2 <- lm(shapes1~logsizes1)$coefficients[2,] - lm(shapes2~logsizes2)$coefficients[2,]
  slopedetr1 <- lm(detr_shapes1~logsizes1)$coefficients[2,]

  results <- all(round(slope1minus2,10) == round(slopedetr1,10))
  expect_true(results)

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


##############################


test_that(desc = "testing expected_shapes", code = {
  data(tails)

  shapes <- tails$shapes
  sizes <- tails$sizes
  model <- lm(geomorph::two.d.array(shapes) ~ sizes)
  predshapes <- expected_shapes(model, xvalue = sizes)

  result <- all(round(model$fitted.values,10) == round(predshapes,10))

  expect_true(result)
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
  phypca <- phy_prcomp(geomorph::two.d.array(consensus(shapes, species)), tree = tree)
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



