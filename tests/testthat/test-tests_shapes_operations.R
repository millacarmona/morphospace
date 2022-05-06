
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



data(tails)
shapes_mat <- geomorph::two.d.array(tails$shapes)

model1 <- lm(shapes_mat ~ tails$sizes)
detr_shapes1 <- geomorph::arrayspecs(detrend_shapes(model1), k = 2, p = 9)



