test_that(desc = "testing expected_shapes, factors, ndims & accuracy, for landmark data", code = {
  data(tails)

  mshape_all <- expected_shapes(tails$shapes)
  ndims <- dim(mshape_all)

  result1 <- length(ndims) == 2
  result2 <- ndims[1] == nrow(tails$shapes)
  result3 <- ndims[2] == ncol(tails$shapes)

  testshape1 <- matrix(colMeans(geomorph::two.d.array(tails$shapes)),ncol=2,byrow=TRUE)
  result4 <- all(round(mshape_all,10) == round(testshape1,10))


  mshape_spp <- expected_shapes(tails$shape, x = tails$data$species)
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
  skip_if_not_installed("caper")


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
  skip_if_not_installed("caper")

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


######################################################################################
######################################################################################

test_that(desc = "testing detrend_shapes, method = orthogonal + default settings,
          for factors and numerics, stats::lm", code = {
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


  detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
  testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, method = "orthogonal")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)


  detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
  testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})

test_that(desc = "testing detrend_shapes, method = orthogonal + default settings,
          for factors and numerics, geomorph::procD.lm", code = {
            data(shells)

            library(geomorph)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, logsizes = logsizes)
            model1 <- geomorph::procD.lm(f1 = shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, method = "orthogonal")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- geomorph::procD.lm(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, method = "orthogonal")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes, method = orthogonal + default settings,
          for factors and numerics, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, logsizes = logsizes)
            model1 <- RRPP::lm.rrpp(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, method = "orthogonal")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- RRPP::lm.rrpp(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, method = "orthogonal")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })

test_that(desc = "testing detrend_shapes, method = orthogonal + default settings,
          for factors and numerics, mvMORPH::mvols", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            model1 <- suppressWarnings(mvMORPH::mvols(shapes ~ logsizes))
            detr_shapes1 <- suppressWarnings(detrend_shapes(model1, method = "orthogonal"))
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            model2 <- suppressWarnings(mvMORPH::mvols(shapes ~ species))
            detr_shapes2 <- suppressWarnings(detrend_shapes(model2, method = "orthogonal"))
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes, method = residuals + default settings,
          for factors and numerics, stats::lm", code = {
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

  detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
  testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, method = "residuals")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
  testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})


test_that(desc = "testing detrend_shapes, method = residuals + default settings,
          for factors and numerics, geomorph::procD.lm", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)


            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = logsizes)
            model1 <- geomorph::procD.lm(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, method = "residuals")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- geomorph::procD.lm(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, method = "residuals")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes, method = residuals + default settings,
          for factors and numerics, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = logsizes)
            model1 <- RRPP::lm.rrpp(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, method = "residuals")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- RRPP::lm.rrpp(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, method = "residuals")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes, method = residuals + default settings,
          for factors and numerics, mvMORPH::mvols", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            model1 <- suppressWarnings(mvMORPH::mvols(shapes ~ logsizes))
            detr_shapes1 <- detrend_shapes(model1, method = "residuals")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            model2 <- suppressWarnings(mvMORPH::mvols(shapes ~ species))
            detr_shapes2 <- detrend_shapes(model2, method = "residuals")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })



test_that(desc = "testing detrend_shapes, method = orthogonal, xvalue,
          for factors and numerics, stats::lm", code = {
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

  detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
  testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                         returnarray = FALSE))
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "orthogonal")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
  testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                returnarray = FALSE)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})

test_that(desc = "testing detrend_shapes, method = orthogonal, xvalue,
          for factors and numerics, geomorph::procD.lm", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, logsizes = logsizes)
            model1 <- geomorph::procD.lm(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "orthogonal")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, x = cbind(logsizes), xvalue = max(logsizes),
                                                   returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- geomorph::procD.lm(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "orthogonal")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                          returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })

test_that(desc = "testing detrend_shapes, method = orthogonal, xvalue,
          for factors and numerics, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, logsizes = logsizes)
            model1 <- RRPP::lm.rrpp(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "orthogonal")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                                   returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- RRPP::lm.rrpp(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "orthogonal")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                          returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })

test_that(desc = "testing detrend_shapes, method = orthogonal, xvalue,
          for factors and numerics, mvMORPH::mvols", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            model1 <- suppressWarnings(mvMORPH::mvols(shapes ~ logsizes))
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "orthogonal")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                                   returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            model2 <- suppressWarnings(mvMORPH::mvols(shapes ~ species))
            detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "orthogonal")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                          returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes, method = residuals, xvalue,
          for factors and numerics, stats::lm", code = {
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

  detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
  testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                         returnarray = FALSE))
  result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



  model2 <- lm(shapes ~ species)
  detr_shapes2 <- detrend_shapes(model2, xvalue = "koeneni", method = "residuals")
  result5 <- all(dim(shapes) == dim(detr_shapes2))

  totvar_raw <- sum(apply(shapes, 2, stats::var))
  totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
  result6 <- (totvar_raw > totvar_detr)

  detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
  testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                returnarray = FALSE)
  result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

  expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

})



test_that(desc = "testing detrend_shapes, method = residuals, xvalue,
          for factors and numerics, geomorph::procD.lm", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, logsizes = logsizes)
            model1 <- geomorph::procD.lm(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "residuals")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                                   returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- geomorph::procD.lm(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "residuals")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                          returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })

test_that(desc = "testing detrend_shapes, method = residuals, xvalue,
          for factors and numerics, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            gdf <- geomorph::geomorph.data.frame(shapes = shapes, logsizes = logsizes)
            model1 <- RRPP::lm.rrpp(shapes ~ logsizes, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "residuals")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                                   returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(shapes = shapes, species = species)
            model2 <- RRPP::lm.rrpp(shapes ~ species, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "residuals")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                          returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })

test_that(desc = "testing detrend_shapes, method = residuals, xvalue,
          for factors and numerics, mvMORPH::mvols", code = {
            data(shells)

            shapes <- shells$shapes$coe
            species <- shells$data$species
            logsizes <- log(shells$sizes)

            model1 <- suppressWarnings(mvMORPH::mvols(shapes ~ logsizes))
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(logsizes), method = "residuals")
            result1 <- all(dim(shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(lm(detr_shapes1 ~ logsizes)$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes, x = logsizes, xvalue = max(logsizes),
                                                   returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            model2 <- suppressWarnings(mvMORPH::mvols(shapes ~ species))
            detr_shapes2 <- detrend_shapes(model2, "koeneni", method = "residuals")
            result5 <- all(dim(shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(shapes, x = species, xvalue = "koeneni",
                                          returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(lm(detr_shapes2 ~ species)$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes, method = orthogonal, newdata,
          for numerics and factors, stats::lm", code = {

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            mod1 <- stats::lm(geomorph::two.d.array(shapes1) ~ logsizes1)
            mod2 <- stats::lm(geomorph::two.d.array(shapes2) ~ logsizes2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)

            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,15) == 0)


            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            expect_true(all(result1, result2, result3))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            mod1 <- stats::lm(geomorph::two.d.array(shapes1) ~ sex1)
            mod2 <- stats::lm(geomorph::two.d.array(shapes2) ~ sex2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            expect_true(all(result4, result5, result6))

})


test_that(desc = "testing detrend_shapes, method = orthogonal, newdata,
          for numerics and factors, geomorph::procD.lm", code = {

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            mod1 <- geomorph::procD.lm(geomorph::two.d.array(shapes1) ~ logsizes1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            mod2 <- geomorph::procD.lm(geomorph::two.d.array(shapes2) ~ logsizes2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,15) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))


            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            mod1 <- geomorph::procD.lm(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            mod2 <- geomorph::procD.lm(geomorph::two.d.array(shapes2) ~ sex2, data = gdf1)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,15) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))

          })


test_that(desc = "testing detrend_shapes, method = orthogonal, newdata,
          for numerics and factors, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            mod1 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes1) ~ logsizes1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            mod2 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes2) ~ logsizes2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,15) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            mod1 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            mod2 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,15) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))

          })

test_that(desc = "testing detrend_shapes, method = orthogonal, newdata,
          for numerics and factors, mvMORPH::mvols", code = {

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            mod1 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes1) ~ logsizes1))
            mod2 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes2) ~ logsizes2))

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,14) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            mod1 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes1) ~ sex1))
            mod2 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes2) ~ sex2))

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))


          })


test_that(desc = "testing detrend_shapes, method = orthogonal, newdata, xvalue,
          for numerics and factors, stats::lm", code = {

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            mod1 <- stats::lm(geomorph::two.d.array(shapes1) ~ logsizes1)
            mod2 <- stats::lm(geomorph::two.d.array(shapes2) ~ logsizes2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,15) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            mod1 <- stats::lm(geomorph::two.d.array(shapes1) ~ sex1)
            mod2 <- stats::lm(geomorph::two.d.array(shapes2) ~ sex2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, xvalue = "H",
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))

          })

test_that(desc = "testing detrend_shapes, method = orthogonal, newdata, xvalue,
          for numerics and factors, geomorph::procD.lm", code = {

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            mod1 <- geomorph::procD.lm(geomorph::two.d.array(shapes1) ~ logsizes1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            mod2 <- geomorph::procD.lm(geomorph::two.d.array(shapes2) ~ logsizes2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,15) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))




            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            mod1 <- geomorph::procD.lm(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            mod2 <- geomorph::procD.lm(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, xvalue = "H",
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))

          })


test_that(desc = "testing detrend_shapes, method = orthogonal, newdata, xvalue,
          for numerics and factors, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            mod1 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes1) ~ logsizes1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            mod2 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes2) ~ logsizes2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,15) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,15) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x =logsizes1, xvalue = max(logsizes2),
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            mod1 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            mod2 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = mod2) %>% arrayspecs(k=2,p=9)
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, xvalue = "H",
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))

          })

test_that(desc = "testing detrend_shapes, method = orthogonal, newdata, xvalue,
          for numerics and factors, mvMORPH::mvols", code = {

            data(shells3D)
            shapes <- shells3D$shapes
            species <- shells3D$data$species
            sizes <- shells3D$sizes
            logsizes <- log(sizes)

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            logsizes1 <- sizes[index1]

            index2 <- species == levels(species)[3]
            shapes2 <- shapes[,,index2]
            logsizes2 <- sizes[index2]

            mod1 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes1) ~ logsizes1))
            mod2 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes2) ~ logsizes2))

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = geomorph::two.d.array(shapes2)) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result1 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = mod2) %>% arrayspecs(k=3,p=90)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result2 <- (round(r,14) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            mod1 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes1) ~ sex1))
            mod2 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes2) ~ sex2))

            detshapes2using1 <- suppressWarnings(
              detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                             newdata = geomorph::two.d.array(shapes2)) %>%
                arrayspecs(k=2,p=9)
            )
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- suppressWarnings(
              detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                             newdata = mod2) %>%
                arrayspecs(k=2,p=9)
            )
            expshapes1 <- expected_shapes(shapes1, x = sex1)
            expshapes2 <- expected_shapes(shapes2, x = sex2)
            r <- cor(prcomp(two.d.array(detshapes2using1))$rotation[,1], prcomp(two.d.array(expshapes1))$rotation[,1])
            result5 <- (round(r,14) == 0)


            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, xvalue = "H",
                                                   returnarray = FALSE))
            result6 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result4, result5, result6))

          })


test_that(desc = "testing detrend_shapes, method = residuals, newdata,
          for numerics and factors, stats::lm", code = {
  data(shells)

  index1 <- shells$data$species == "koeneni"
  shapes1 <- shells$shapes$coe[index1,]
  logsizes1 <- log(shells$sizes[index1])

  index2 <- shells$data$species == "windhauseni"
  shapes2 <- shells$shapes$coe[index2,]
  logsizes2 <- log(shells$sizes[index2])

  model_for1 <- stats::lm(shapes1 ~ logsizes1)
  model_for2 <- stats::lm(shapes2 ~ logsizes2)
  detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

  slope1 <- model_for1$coefficients[2,]
  slope2 <- model_for2$coefficients[2,]

  slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
  slopedetr1 <- stats::lm(detr_shapes1 ~ logsizes1)$coefficients[2,]

  result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

  detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
  testshape1 <- colMeans(expected_shapes(shapes2,
                                         returnarray = FALSE))
  result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
  expect_true(all(result1, result2))



  data(tails)
  shapes <- tails$shapes
  species <- tails$data$species
  sex <- tails$data$sex

  index1 <- species == levels(species)[7]
  shapes1 <- shapes[,,index1]
  sex1 <- sex[index1]

  index2 <- species == levels(species)[10]
  shapes2 <- shapes[,,index2]
  sex2 <- sex[index2]

  model_for1 <- stats::lm(geomorph::two.d.array(shapes1) ~ sex1)
  model_for2 <- stats::lm(geomorph::two.d.array(shapes2) ~ sex2)
  detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

  slope1 <- model_for1$coefficients[2,]
  slope2 <- model_for2$coefficients[2,]

  slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
  slopedetr1 <- stats::lm(detr_shapes1 ~ sex1)$coefficients[2,]

  result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

  detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
  testshape2 <- colMeans(expected_shapes(shapes2,
                                         returnarray = FALSE))


  result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

  expect_true(all(result3, result4))

})

test_that(desc = "testing detrend_shapes, method = residuals, newdata,
          for numerics and factors, geomorph::procD.lm", code = {
            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            model_for1 <- geomorph::procD.lm(shapes1 ~ logsizes1, data = gdf1)
            model_for2 <- geomorph::procD.lm(shapes2 ~ logsizes2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, logsizes1 = logsizes1)
            slopedetr1 <- geomorph::procD.lm(detr_shapes1 ~ logsizes1, data = gdf)$coefficients[2,]

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            model_for1 <- geomorph::procD.lm(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            model_for2 <- geomorph::procD.lm(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sex1 = sex1)
            slopedetr1 <- geomorph::procD.lm(detr_shapes1 ~ sex1, data = gdf)$coefficients[2,]

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2,
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })


test_that(desc = "testing detrend_shapes, method = residuals, newdata,
          for numerics and factors, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            model_for1 <- RRPP::lm.rrpp(shapes1 ~ logsizes1, data = gdf1)
            model_for2 <- RRPP::lm.rrpp(shapes2 ~ logsizes2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

            slope1 <- model_for1$LM$coefficients[2,]
            slope2 <- model_for2$LM$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, logsizes1 = logsizes1)
            slopedetr1 <- RRPP::lm.rrpp(detr_shapes1 ~ logsizes1, data = gdf)$LM$coefficients[2,]

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            model_for1 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            model_for2 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

            slope1 <- model_for1$LM$coefficients[2,]
            slope2 <- model_for2$LM$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sex1 = sex1)
            slopedetr1 <- RRPP::lm.rrpp(detr_shapes1 ~ sex1, data = gdf)$LM$coefficients

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2,
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })


test_that(desc = "testing detrend_shapes, method = residuals, newdata,
          for numerics and factors, mvMORPH::mvols", code = {
            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            model_for1 <- suppressWarnings(mvMORPH::mvols(shapes1 ~ logsizes1))
            model_for2 <- suppressWarnings(mvMORPH::mvols(shapes2 ~ logsizes2))
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            slopedetr1 <- suppressWarnings(mvMORPH::mvols(detr_shapes1 ~ logsizes1)$coefficients[2,])

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            model_for1 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes1) ~ sex1))
            model_for2 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes2) ~ sex2))
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            slopedetr1 <- suppressWarnings(mvMORPH::mvols(detr_shapes1 ~ sex1)$coefficients[2,])

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2,
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })


test_that(desc = "testing detrend_shapes, method = residuals, newdata, xvalue,
          for numerics and factors, stats::lm", code = {
            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            model_for1 <- stats::lm(shapes1 ~ logsizes1)
            model_for2 <- stats::lm(shapes2 ~ logsizes2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, xvalue = max(logsizes1),
                                           method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            slopedetr1 <- stats::lm(detr_shapes1 ~ logsizes1)$coefficients[2,]

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2, x = logsizes2, xvalue = max(logsizes1),
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))


            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            model_for1 <- stats::lm(geomorph::two.d.array(shapes1) ~ sex1)
            model_for2 <- stats::lm(geomorph::two.d.array(shapes2) ~ sex2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, xvalue = "H",
                                           method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            slopedetr1 <- stats::lm(detr_shapes1 ~ sex1)$coefficients[2,]

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2, x = sex2, xvalue = "H",
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })

test_that(desc = "testing detrend_shapes, method = residuals, newdata, xvalue,
          for numerics and factors, geomorph::procD.lm", code = {
            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            model_for1 <- geomorph::procD.lm(shapes1 ~ logsizes1, data = gdf1)
            model_for2 <- geomorph::procD.lm(shapes2 ~ logsizes2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, xvalue = max(logsizes1),
                                           method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, logsizes1 = logsizes1)
            slopedetr1 <- geomorph::procD.lm(detr_shapes1 ~ logsizes1, data = gdf)$coefficients[2,]

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2, x = logsizes2, xvalue = max(logsizes1),
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            model_for1 <- geomorph::procD.lm(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            model_for2 <- geomorph::procD.lm(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, xvalue = "H",
                                           method = "residuals")

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sex1 = sex1)
            slopedetr1 <- geomorph::procD.lm(detr_shapes1 ~ sex1, data = gdf)$coefficients[2,]

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2, x = sex2, xvalue = "H",
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })



test_that(desc = "testing detrend_shapes, method = residuals, newdata, xvalue,
          for numerics and factors, RRPP::lm.rrpp", code = {
            skip_if_not_installed("RRPP")

            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = logsizes1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = logsizes2)
            model_for1 <- RRPP::lm.rrpp(shapes1 ~ logsizes1, data = gdf1)
            model_for2 <- RRPP::lm.rrpp(shapes2 ~ logsizes2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, xvalue = max(logsizes1),
                                           method = "residuals")

            slope1 <- model_for1$LM$coefficients[2,]
            slope2 <- model_for2$LM$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, logsizes1 = logsizes1)
            slopedetr1 <- RRPP::lm.rrpp(detr_shapes1 ~ logsizes1, data = gdf)$LM$coefficients[2,]

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2, x = logsizes2, xvalue = max(logsizes1),
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2)
            model_for1 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes1) ~ sex1, data = gdf1)
            model_for2 <- RRPP::lm.rrpp(geomorph::two.d.array(shapes2) ~ sex2, data = gdf2)
            detr_shapes1 <- detrend_shapes(model_for2, newdata = model_for1, xvalue = "H",
                                           method = "residuals")

            slope1 <- model_for1$LM$coefficients[2,]
            slope2 <- model_for2$LM$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sex1 = sex1)
            slopedetr1 <- RRPP::lm.rrpp(detr_shapes1 ~ sex1, data = gdf)$LM$coefficients

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2, x = sex2, xvalue = "H",
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })


test_that(desc = "testing detrend_shapes, method = residuals, newdata, xvalue,
          for numerics and factors, mvMORPH::mvols", code = {
            data(shells)

            index1 <- shells$data$species == "koeneni"
            shapes1 <- shells$shapes$coe[index1,]
            logsizes1 <- log(shells$sizes[index1])

            index2 <- shells$data$species == "windhauseni"
            shapes2 <- shells$shapes$coe[index2,]
            logsizes2 <- log(shells$sizes[index2])

            model_for1 <- suppressWarnings(mvMORPH::mvols(shapes1 ~ logsizes1))
            model_for2 <- suppressWarnings(mvMORPH::mvols(shapes2 ~ logsizes2))
            detr_shapes1 <- suppressWarnings(
              detrend_shapes(model_for2, newdata = model_for1, xvalue = max(logsizes1),
                                           method = "residuals")
              )

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            slopedetr1 <- suppressWarnings(mvMORPH::mvols(detr_shapes1 ~ logsizes1)$coefficients[2,])

            result1 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes2, x = logsizes2, xvalue = max(logsizes1),
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)
            shapes <- tails$shapes
            species <- tails$data$species
            sex <- tails$data$sex

            index1 <- species == levels(species)[7]
            shapes1 <- shapes[,,index1]
            sex1 <- sex[index1]

            index2 <- species == levels(species)[10]
            shapes2 <- shapes[,,index2]
            sex2 <- sex[index2]

            model_for1 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes1) ~ sex1))
            model_for2 <- suppressWarnings(mvMORPH::mvols(geomorph::two.d.array(shapes2) ~ sex2))
            detr_shapes1 <- suppressWarnings(
              detrend_shapes(model_for2, newdata = model_for1, xvalue = "H",
                                           method = "residuals")
            )

            slope1 <- model_for1$coefficients[2,]
            slope2 <- model_for2$coefficients[2,]

            slope1minus2 <- model_for1$coefficients[2,] - model_for2$coefficients[2,]
            slopedetr1 <- suppressWarnings(mvMORPH::mvols(detr_shapes1 ~ sex1)$coefficients[2,])

            result3 <- all(round(slope1minus2,10) == round(slopedetr1,10))

            detr_mshape2 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes2, x = sex2, xvalue = "H",
                                                   returnarray = FALSE))
            result4 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            expect_true(all(result3, result4))

          })

#####################################################################################


test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal + default settings,
          for factors and numerics, geomorph::procD.pgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- tapply(type, species, unique)
            sp_type[sp_type == 1] <- "DF" ; sp_type[sp_type == 2] <- "NDF"
            sp_type <- factor(sp_type)
            tree <- tails$tree

            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_logsizes = cbind(sp_logsizes),
                                                 phy = tree)
            model1 <- geomorph::procD.pgls(sp_shapes ~ sp_logsizes, phy = phy, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, method = "orthogonal")
            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1,
                                                 sp_logsizes = cbind(sp_logsizes),
                                                 phy = tree)
            result3 <- all(round(geomorph::procD.pgls(detr_shapes1 ~ sp_logsizes, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes,
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_type = sp_type, phy = tree)
            model2 <- geomorph::procD.pgls(sp_shapes ~ sp_type, phy = phy, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, method = "orthogonal")
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)


            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(sp_shapes, x = cbind(sp_type),
                                                   tree = tree, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            gdf <- geomorph::geomorph.data.frame(detr_shapes2 = detr_shapes2,
                                                 sp_type = cbind(sp_type),
                                                 phy = tree)
            result8 <- all(round(geomorph::procD.pgls(detr_shapes2 ~ sp_type, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal + default settings,
          for factors and numerics, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- tapply(type, species, unique)
            sp_type[sp_type == 1] <- "DF" ; sp_type[sp_type == 2] <- "NDF"
            sp_type <- factor(sp_type)
            tree <- tails$tree


            model1 <- mvMORPH::mvgls(sp_shapes ~ sp_logsizes, tree = tree, model = "BM")
            detr_shapes1 <- detrend_shapes(model1, method = "orthogonal")
            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(mvMORPH::mvgls(detr_shapes1 ~ sp_logsizes, tree = tree, model = "BM")$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes,
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            model2 <- mvMORPH::mvgls(sp_shapes ~ sp_type, tree = tree, model = "BM")
            detr_shapes2 <- detrend_shapes(model2, method = "orthogonal")
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(sp_shapes, x = cbind(sp_type),
                                                   tree = tree, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))

            result8 <- all(round(mvMORPH::mvgls(detr_shapes2 ~ sp_type, tree = tree, model = "BM")$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals + default settings,
          for factors and numerics, geomorph::procD.pgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_logsizes = cbind(sp_logsizes), phy = tree)
            model1 <- geomorph::procD.pgls(sp_shapes ~ sp_logsizes, phy = phy, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, method = "residuals")
            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sp_logsizes = cbind(sp_logsizes), phy = tree)
            result3 <- all(round(geomorph::procD.pgls(detr_shapes1 ~ sp_logsizes, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)


            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes,
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))


            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_type = sp_type, phy = tree)
            model2 <- geomorph::procD.pgls(sp_shapes ~ sp_type, phy = phy, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, method = "residuals")
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(sp_shapes, x = cbind(sp_type),
                                                   tree = tree, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            gdf <- geomorph::geomorph.data.frame(detr_shapes2 = detr_shapes2, sp_type = sp_type, phy = tree)
            result8 <- all(round(geomorph::procD.pgls(detr_shapes2 ~ sp_type, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals + default settings,
          for factors and numerics, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            model1 <- mvMORPH::mvgls(sp_shapes ~ sp_logsizes, tree = tree, model = "BM")
            detr_shapes1 <- detrend_shapes(model1, method = "residuals")
            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(mvMORPH::mvgls(detr_shapes1 ~ sp_logsizes, tree = tree, model = "BM")$coefficients[2,], 10) == 0)


            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes,
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            model2 <- mvMORPH::mvgls(sp_shapes ~ sp_type, tree = tree, model = "BM")
            detr_shapes2 <- detrend_shapes(model2, method = "residuals")
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(sp_shapes, x = cbind(sp_type),
                                                   tree = tree, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            result8 <- all(round(mvMORPH::mvgls(detr_shapes2 ~ sp_type, tree = tree, model = "BM")$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })

test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal, xvalue,
          for factors and numerics, geomorph::procD.pgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_logsizes = cbind(sp_logsizes), phy = tree)
            model1 <- geomorph::procD.pgls(sp_shapes ~ sp_logsizes, phy = phy, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(sp_logsizes),
                                           evmodel = "BM", tree = tree, method = "orthogonal")
            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sp_logsizes = cbind(sp_logsizes), phy = tree)
            result3 <- all(round(geomorph::procD.pgls(detr_shapes1 ~ sp_logsizes, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes, xvalue = max(sp_logsizes),
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_type = sp_type, phy = tree)
            model2 <- geomorph::procD.pgls(sp_shapes ~ sp_type, phy = phy, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, xvalue = "1", method = "orthogonal",
                                           evmodel = "BM", tree = tree)
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)


            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(sp_shapes, x = sp_type, xvalue = "1",
                                          tree = tree, returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            gdf <- geomorph::geomorph.data.frame(detr_shapes2 = detr_shapes2, sp_type = sp_type, phy = tree)
            result8 <- all(round(geomorph::procD.pgls(detr_shapes2 ~ sp_type, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })



test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal, xvalue,
          for factors and numerics, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            model1 <- mvMORPH::mvgls(sp_shapes ~ sp_logsizes, tree = tree, model = "BM")
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(sp_logsizes),
                                           evmodel = "BM", tree = tree, method = "orthogonal")
            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(mvMORPH::mvgls(detr_shapes1 ~ sp_logsizes, tree = tree, model = "BM")$coefficients[2,], 10) == 0)

            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes, xvalue = max(sp_logsizes),
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            model2 <- mvMORPH::mvgls(sp_shapes ~ sp_type, tree = tree, model = "BM")
            detr_shapes2 <- suppressWarnings(
              detrend_shapes(model2, xvalue = "1", method = "orthogonal",
                             evmodel = "BM", tree = tree)
            )
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)


            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- expected_shapes(sp_shapes, x = sp_type, xvalue = "1",
                                          tree = tree, returnarray = FALSE)
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            result8 <- all(round(mvMORPH::mvgls(detr_shapes2 ~ sp_type, tree = tree, model = "BM")$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))
          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals, xvalue,
          for factors and numerics, geomorph::procD.pgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_logsizes = cbind(sp_logsizes), phy = tree)
            model1 <- geomorph::procD.pgls(sp_shapes ~ sp_logsizes, phy = phy, data = gdf)
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(sp_logsizes), method = "residuals",
                                           tree = tree, evmodel = "BM")

            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            gdf <- geomorph::geomorph.data.frame(detr_shapes1 = detr_shapes1, sp_logsizes = cbind(sp_logsizes), phy = tree)
            result3 <- all(round(geomorph::procD.pgls(detr_shapes1 ~ sp_logsizes, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)


            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes, xvalue = max(sp_logsizes),
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            gdf <- geomorph::geomorph.data.frame(sp_shapes = sp_shapes, sp_type = sp_type, phy = tree)
            model2 <- geomorph::procD.pgls(sp_shapes ~ sp_type, phy = phy, data = gdf)
            detr_shapes2 <- detrend_shapes(model2, xvalue = "1", method = "residuals", tree = tree, evmodel = "BM")
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(sp_shapes, x = sp_type, xvalue = "1",
                                                   tree = tree, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            gdf <- geomorph::geomorph.data.frame(detr_shapes2 = detr_shapes2, sp_type = sp_type, phy = tree)
            result8 <- all(round(geomorph::procD.pgls(detr_shapes2 ~ sp_type, phy = phy, data = gdf)$pgls.coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals, xvalue,
          for factors and numerics, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            logsizes <- log(tails$sizes)
            type <- tails$data$type

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            model1 <- mvMORPH::mvgls(sp_shapes ~ sp_logsizes, tree = tree, model = "BM")
            detr_shapes1 <- detrend_shapes(model1, xvalue = max(sp_logsizes), method = "residuals",
                                           tree = tree, evmodel = "BM")

            result1 <- all(dim(sp_shapes) == dim(detr_shapes1))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes1, 2, stats::var))
            result2 <- (totvar_raw > totvar_detr)

            result3 <- all(round(mvMORPH::mvgls(detr_shapes1 ~ sp_logsizes, tree = tree, model = "BM")$coefficients[2,], 10) == 0)


            detr_mshape1 <- expected_shapes(detr_shapes1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(sp_shapes, x = sp_logsizes, xvalue = max(sp_logsizes),
                                                   tree = tree, returnarray = FALSE))
            result4 <- all(round(testshape1, 10) == round(detr_mshape1,10))



            model2 <- mvMORPH::mvgls(sp_shapes ~ sp_type, tree = tree, model = "BM")
            suppressWarnings(
              detr_shapes2 <- detrend_shapes(model2, xvalue = "1", method = "residuals", tree = tree, evmodel = "BM")
            )
            result5 <- all(dim(sp_shapes) == dim(detr_shapes2))

            totvar_raw <- sum(apply(sp_shapes, 2, stats::var))
            totvar_detr <- sum(apply(detr_shapes2, 2, stats::var))
            result6 <- (totvar_raw > totvar_detr)

            detr_mshape2 <- expected_shapes(detr_shapes2, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(sp_shapes, x = sp_type, xvalue = "1",
                                                   tree = tree, returnarray = FALSE))
            result7 <- all(round(testshape2, 10) == round(detr_mshape2,10))


            result8 <- all(round(mvMORPH::mvgls(detr_shapes2 ~ sp_type, tree = tree, model = "BM")$coefficients[2,], 10) == 0)

            expect_true(all(result1, result2, result3, result4, result5, result6, result7, result8))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal, newdata,
          for numerics, geomorph::procD.pgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = cbind(logsizes1), phy = tree1)
            mod1 <- geomorph::procD.pgls(shapes1 ~ logsizes1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = cbind(logsizes2), phy = tree2)
            mod2 <- geomorph::procD.pgls(shapes2 ~ logsizes2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = shapes2)

            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result1 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result2 <- (round(r,14) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, tree = tree1,
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1, phy = tree1)
            mod1 <- geomorph::procD.pgls(shapes1 ~ sex1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2, phy = tree2)
            mod2 <- geomorph::procD.pgls(shapes2 ~ sex2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = shapes2)

            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result4 <- (round(r,14) == 0)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2)
            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, tree = tree1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result4, result5, result6))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal, newdata,
          for numerics, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            mod1 <- mvMORPH::mvgls(shapes1 ~ logsizes1, tree = tree1, model = "BM")
            mod2 <- mvMORPH::mvgls(shapes2 ~ logsizes2, tree = tree2, model = "BM")

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = shapes2)

            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result1 <- (round(r,14) == 0)


            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                               newdata = mod2)
            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result2 <- (round(r,14) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, tree = tree1,
                                                   returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))




            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            mod1 <- mvMORPH::mvgls(shapes1 ~ sex1, tree = tree1, model = "BM")
            mod2 <- mvMORPH::mvgls(shapes2 ~ sex2, tree = tree2, model = "BM")

            detshapes2using1 <- suppressWarnings(
              detrend_shapes(mod1, method = "orthogonal",
                             newdata = shapes2)
            )

            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result4 <- (round(r,14) == 0)

            suppressWarnings(
              detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal",
                                                 newdata = mod2)
            )

            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, tree = tree1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result4, result5, result6))
          })

test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal, newdata, xvalue,
          for numerics and factors, geomorph::procD.pgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = cbind(logsizes1), phy = tree1)
            mod1 <- geomorph::procD.pgls(shapes1 ~ logsizes1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = cbind(logsizes2), phy = tree2)
            mod2 <- geomorph::procD.pgls(shapes2 ~ logsizes2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = shapes2, tree = tree1, evmodel = "BM")

            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result1 <- (round(r,14) == 0)



            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = mod2, tree = tree1, evmodel = "BM")
            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result2 <- (round(r,14) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   tree = tree1, returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1, phy = tree1)
            mod1 <- geomorph::procD.pgls(shapes1 ~ sex1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2, phy = tree2)
            mod2 <- geomorph::procD.pgls(shapes2 ~ sex2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = shapes2, tree = tree1, evmodel = "BM")

            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result4 <- (round(r,14) == 0)


            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                               newdata = mod2, tree = tree1, evmodel = "BM")
            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, tree = tree1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result4, result5, result6))
          })


test_that(desc = "testing detrend_shapes with phylogeny, method = orthogonal, newdata, xvalue,
          for numerics and factors, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            mod1 <- mvMORPH::mvgls(shapes1 ~ logsizes1, tree = tree1, model = "BM")
            mod2 <- mvMORPH::mvgls(shapes2 ~ logsizes2, tree = tree2, model = "BM")

            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = shapes2, tree = tree1, evmodel = "BM")

            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result1 <- (round(r,14) == 0)


            detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = max(logsizes2),
                                               newdata = mod2, tree = tree1, evmodel = "BM")
            expshapes1 <- expected_shapes(shapes1, x = logsizes1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = logsizes2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result2 <- (round(r,14) == 0)

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   tree = tree1, returnarray = FALSE))
            result3 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result1, result2, result3))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            mod1 <- mvMORPH::mvgls(shapes1 ~ sex1, tree = tree1, model = "BM")
            mod2 <- mvMORPH::mvgls(shapes2 ~ sex2, tree = tree2, model = "BM")

            detshapes2using1 <- suppressWarnings(
              detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                             newdata = shapes2, tree = tree1, evmodel = "BM")
            )

            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result4 <- (round(r,14) == 0)

            suppressWarnings(
              detshapes2using1 <- detrend_shapes(mod1, method = "orthogonal", xvalue = "H",
                                                 newdata = mod2, tree = tree1, evmodel = "BM")
            )

            expshapes1 <- expected_shapes(shapes1, x = sex1, tree = tree1, returnarray = FALSE)
            expshapes2 <- expected_shapes(shapes2, x = sex2, tree = tree2, returnarray = FALSE)
            r <- cor(prcomp(detshapes2using1)$rotation[,1], prcomp(expshapes1)$rotation[,1])
            result5 <- (round(r,14) == 0)

            detr_mshape2 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape2 <- colMeans(expected_shapes(shapes1, x = sex1, tree = tree1,
                                                   returnarray = FALSE))
            result6 <- all(round(testshape1, 10) == round(detr_mshape1,10))

            expect_true(all(result4, result5, result6))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals, newdata,
          for numerics and factors, geomorph::procD.pgls", code = {

            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = cbind(logsizes1), phy = tree1)
            model_for1 <- geomorph::procD.pgls(shapes1 ~ logsizes1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = cbind(logsizes2), phy = tree2)
            model_for2 <- geomorph::procD.pgls(shapes2 ~ logsizes2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals",
                                               newdata = model_for2)


            slope2minus1 <- model_for2$pgls.coefficients[2,] - model_for1$pgls.coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detshapes2using1 = detshapes2using1,
                                                 logsizes2 = cbind(logsizes2), phy = tree2)
            slopedetr2 <- geomorph::procD.pgls(detshapes2using1 ~ logsizes2, phy = phy, data = gdf)$pgls.coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, tree = tree1,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = cbind(sex1), phy = tree1)
            model_for1 <- geomorph::procD.pgls(shapes1 ~ sex1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = cbind(sex2), phy = tree2)
            model_for2 <- geomorph::procD.pgls(shapes2 ~ sex2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals",
                                               newdata = model_for2)

            slope2minus1 <- model_for2$pgls.coefficients[2,] - model_for1$pgls.coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detshapes2using1 = detshapes2using1,
                                                 logsizes2 = cbind(sex2), phy = tree2)
            slopedetr2 <- geomorph::procD.pgls(detshapes2using1 ~ sex2, phy = phy, data = gdf)$pgls.coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = sex1, tree = tree1,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals, newdata,
          for numerics and factors, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            model_for1 <- mvMORPH::mvgls(shapes1 ~ logsizes1, tree = tree1, model = "BM")
            model_for2 <- mvMORPH::mvgls(shapes2 ~ logsizes2, tree = tree2, model = "BM")

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals",
                                               newdata = model_for2)

            slope2minus1 <- model_for2$coefficients[2,] - model_for1$coefficients[2,]
            slopedetr2 <- mvMORPH::mvgls(detshapes2using1 ~ logsizes2, tree = tree2, model = "BM")$coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, tree = tree1,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            model_for1 <- mvMORPH::mvgls(shapes1 ~ sex1, tree = tree1, model = "BM")
            model_for2 <- mvMORPH::mvgls(shapes2 ~ sex2, tree = tree2, model = "BM")

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals",
                                               newdata = model_for2)


            slope2minus1 <- model_for2$pgls.coefficients[2,] - model_for1$pgls.coefficients[2,]
            slopedetr2 <- mvMORPH::mvgls(detshapes2using1 ~ sex2, tree = tree2, model = "BM")$coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = sex1, tree = tree1,
                                                   returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals, newdata, xvalue,
          for numerics and factors, geomorph::procD.pgls", code = {

            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, logsizes1 = cbind(logsizes1), phy = tree1)
            model_for1 <- geomorph::procD.pgls(shapes1 ~ logsizes1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, logsizes2 = cbind(logsizes2), phy = tree2)
            model_for2 <- geomorph::procD.pgls(shapes2 ~ logsizes2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals", xvalue = max(logsizes2),
                                               newdata = model_for2, tree = tree1, evmodel = "BM")


            slope2minus1 <- model_for2$pgls.coefficients[2,] - model_for1$pgls.coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detshapes2using1 = detshapes2using1,
                                                 logsizes2 = cbind(logsizes2), phy = tree2)
            slopedetr2 <- geomorph::procD.pgls(detshapes2using1 ~ logsizes2, phy = phy, data = gdf)$pgls.coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   tree = tree1, returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            gdf1 <- geomorph::geomorph.data.frame(shapes1 = shapes1, sex1 = sex1, phy = tree1)
            model_for1 <- geomorph::procD.pgls(shapes1 ~ sex1, phy = phy, data = gdf1)
            gdf2 <- geomorph::geomorph.data.frame(shapes2 = shapes2, sex2 = sex2, phy = tree2)
            model_for2 <- geomorph::procD.pgls(shapes2 ~ sex2, phy = phy, data = gdf2)

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals", xvalue = "H",
                                               newdata = model_for2, tree = tree1, evmodel = "BM")

            slope2minus1 <- model_for2$pgls.coefficients[2,] - model_for1$pgls.coefficients[2,]
            gdf <- geomorph::geomorph.data.frame(detshapes2using1 = detshapes2using1,
                                                 logsizes2 = cbind(sex2), phy = tree2)
            slopedetr2 <- geomorph::procD.pgls(detshapes2using1 ~ sex2, phy = phy, data = gdf)$pgls.coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = sex1, xvalue = "H",
                                                   tree = tree1, returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))

          })


test_that(desc = "testing detrend_shapes with phylogeny, method = residuals, newdata, xvalue,
          for numerics and factors, mvMORPH::mvgls", code = {
            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type
            logsizes <- log(tails$sizes)

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type <- factor(type)

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            sp_logsizes <- tapply(logsizes, species, mean)
            sp_type <- factor(tapply(type, species, unique))
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            logsizes1 <- sp_logsizes[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            logsizes2 <- sp_logsizes[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            model_for1 <- mvMORPH::mvgls(shapes1 ~ logsizes1, tree = tree1, model = "BM")
            model_for2 <- mvMORPH::mvgls(shapes2 ~ logsizes2, tree = tree2, model = "BM")

            detshapes2using1 <- detrend_shapes(model_for1, method = "residuals", xvalue = max(logsizes2),
                                               newdata = model_for2, tree = tree1, evmodel = "BM")

            slope2minus1 <- model_for2$coefficients[2,] - model_for1$coefficients[2,]
            slopedetr2 <- mvMORPH::mvgls(detshapes2using1 ~ logsizes2, tree = tree2, model = "BM")$coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = logsizes1, xvalue = max(logsizes2),
                                                   tree = tree1, returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))



            data(tails)

            shapes <- geomorph::two.d.array(tails$shapes)
            species <- tails$data$species
            type <- tails$data$type

            type <- as.character(type)
            type[which(species == "T. vociferans")] <- "DF"
            type[which(species == "T. crassirostris")] <- "DF"
            type <- factor(type)

            sp_type <- factor(tapply(type, species, unique))
            sp_sex <- factor(rep(c("H","M"), length.out = 13))

            sp_shapes <- expected_shapes(shapes, species, returnarray = FALSE)
            tree <- tails$tree


            index1 <- sp_type == levels(sp_type)[1]
            shapes1 <- sp_shapes[index1,]
            sex1 <- sp_sex[index1]

            index2 <- sp_type == levels(sp_type)[2]
            shapes2 <- sp_shapes[index2,]
            sex2 <- sp_sex[index2]

            tree1 <- ape::drop.tip(phy = tree, tip = rownames(shapes2))
            tree2 <- ape::drop.tip(phy = tree, tip = rownames(shapes1))


            model_for1 <- mvMORPH::mvgls(shapes1 ~ sex1, tree = tree1, model = "BM")
            model_for2 <- mvMORPH::mvgls(shapes2 ~ sex2, tree = tree2, model = "BM")

            suppressWarnings(
              detshapes2using1 <- detrend_shapes(model_for1, method = "residuals", xvalue = "H",
                                                 newdata = model_for2, tree = tree1, evmodel = "BM")
            )

            slope2minus1 <- model_for2$pgls.coefficients[2,] - model_for1$pgls.coefficients[2,]
            slopedetr2 <- mvMORPH::mvgls(detshapes2using1 ~ sex2, tree = tree2, model = "BM")$coefficients[2,]

            result1 <- all(round(slope2minus1,10) == round(slopedetr2,10))

            detr_mshape1 <- expected_shapes(detshapes2using1, returnarray = FALSE)
            testshape1 <- colMeans(expected_shapes(shapes1, x = sex1, xvalue = "H",
                                                   tree = tree1, returnarray = FALSE))
            result2 <- all(round(testshape1, 10) == round(detr_mshape1,10))
            expect_true(all(result1, result2))

          })


######################################################################################
######################################################################################

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

######################################################################################
######################################################################################


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


######################################################################################
######################################################################################

test_that(desc = "testing ax_transformations, dimensions", code = {
  data(tails)

  shapes <- tails$shapes
  sizes <- tails$sizes
  species <- tails$data$species
  tree <- tails$tree

  model1 <- stats::lm(geomorph::two.d.array(shapes) ~ sizes)
  gdf <- geomorph::geomorph.data.frame(shapes = geomorph::two.d.array(shapes), sizes = sizes)
  model2 <- geomorph::procD.lm(shapes ~ sizes, data= gdf)
  model3 <- suppressWarnings(
    mvMORPH::mvols(geomorph::two.d.array(shapes) ~ sizes)
  )
  model4 <- RRPP::lm.rrpp(shapes ~ sizes, data= gdf)
  gdf2 <- geomorph::geomorph.data.frame(shapes = expected_shapes(shapes, species, returnarray = FALSE),
                                        sizes = cbind(tapply(sizes,species,mean)), phy = tree)
  model5 <- geomorph::procD.pgls(shapes ~ sizes, phy = phy, data= gdf2)
  model6 <- mvMORPH::mvgls(expected_shapes(shapes, species, returnarray = FALSE) ~ cbind(tapply(sizes,species,mean)),
                           model = "BM", tree = tree)

  pca1 <- stats::prcomp(geomorph::two.d.array(shapes))
  pca2 <- geomorph::gm.prcomp(geomorph::two.d.array(shapes))
  pca3 <- suppressWarnings(
    mvMORPH::mvgls.pca(mvMORPH::mvols(
      geomorph::two.d.array(shapes) ~ sizes), plot = FALSE)
  )
  pca3$center <- colMeans(geomorph::two.d.array(shapes))
  bgpca1 <- bg_prcomp(geomorph::two.d.array(shapes), groups = species)
  bgpca2 <- Morpho::groupPCA(geomorph::two.d.array(shapes), groups = species)
  phypca1 <- phytools::phyl.pca(geomorph::two.d.array(expected_shapes(shapes, species)), tree = tree)
  phypca2 <- geomorph::gm.prcomp(geomorph::two.d.array(expected_shapes(shapes, species)), phy = tree)
  phypca3 <- mvMORPH::mvgls.pca(mvMORPH::mvgls(
    expected_shapes(shapes, species, returnarray = FALSE) ~ cbind(tapply(sizes,species,mean)),
    model = "BM", tree = tree), plot = FALSE)
  phypca3$center <- colMeans(mvMORPH::mvgls(
    expected_shapes(shapes, species, returnarray = FALSE) ~ cbind(tapply(sizes,species,mean)),
    model = "BM", tree = tree)$fitted)
  pls1 <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = sizes)
  pls2 <- Morpho::pls2B(y = geomorph::two.d.array(shapes), x = sizes)
  pls3 <- geomorph::two.b.pls(A2 = geomorph::two.d.array(shapes), A1 = sizes)
  paca <- geomorph::gm.prcomp(geomorph::two.d.array(expected_shapes(shapes, species)),
                              phy = tree, align.to.phy = TRUE)


  ext1 <- ax_transformation(obj = model1)
  ext2 <- ax_transformation(obj = model2)
  ext3 <- ax_transformation(obj = model3)
  ext4 <- ax_transformation(obj = model4)
  ext5 <- ax_transformation(obj = model5)
  ext6 <- ax_transformation(obj = model5)

  ext7 <- ax_transformation(obj = pca1)
  ext8 <- ax_transformation(obj = pca2)
  ext9 <- ax_transformation(obj = pca3)
  ext10 <- ax_transformation(obj = bgpca1)
  ext11 <- ax_transformation(obj = bgpca2)
  ext12 <- ax_transformation(obj = phypca1)
  ext13 <- ax_transformation(obj = phypca2)
  ext14 <- ax_transformation(obj = phypca3)
  ext15 <- ax_transformation(obj = pls1)
  ext16 <- ax_transformation(obj = pls2)
  ext17 <- ax_transformation(obj = pls3)
  ext18 <- ax_transformation(obj = paca)

  results1 <- all(c(dim(ext1)[1], dim(ext2)[1], dim(ext3)[1], dim(ext4)[1], dim(ext5)[1],
                    dim(ext6)[1], dim(ext7)[1], dim(ext8)[1], dim(ext9)[1], dim(ext10)[1], dim(ext11)[1],
                    dim(ext12)[1], dim(ext13)[1], dim(ext14)[1], dim(ext15)[1], dim(ext16)[1], dim(ext17)[1],
                    dim(ext18)[1]) == 2)
  results2 <- all(c(dim(ext1)[2], dim(ext2)[2], dim(ext3)[2], dim(ext4)[2], dim(ext5)[2],
                    dim(ext6)[2], dim(ext7)[2], dim(ext8)[2], dim(ext9)[2], dim(ext10)[2], dim(ext11)[2],
                    dim(ext12)[2], dim(ext13)[2], dim(ext14)[2], dim(ext15)[2], dim(ext16)[2], dim(ext17)[2],
                    dim(ext18)[2]) == nrow(shapes) * ncol(shapes))
  expect_true(all(results1, results2))
})


test_that(desc = "testing ax_transformations, accuracy of shapes (basic)", code = {
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


test_that(desc = "testing ax_transformations, accuracy of shapes (extended)", code = {
  data(tails)

  shapes <- tails$shapes
  sizes <- tails$sizes
  species <- tails$data$species
  tree <- tails$tree

  model1 <- stats::lm(geomorph::two.d.array(shapes) ~ sizes)
  gdf <- geomorph::geomorph.data.frame(shapes = geomorph::two.d.array(shapes), sizes = sizes)
  model2 <- geomorph::procD.lm(shapes ~ sizes, data= gdf)
  model3 <- suppressWarnings(
    mvMORPH::mvols(geomorph::two.d.array(shapes) ~ sizes)
  )
  model4 <- RRPP::lm.rrpp(shapes ~ sizes, data= gdf)
  gdf2 <- geomorph::geomorph.data.frame(shapes = expected_shapes(shapes, species, returnarray = FALSE),
                                        sizes = cbind(tapply(sizes,species,mean)), phy = tree)
  model5 <- geomorph::procD.pgls(shapes ~ sizes, phy = phy, data= gdf2)
  model6 <- mvMORPH::mvgls(expected_shapes(shapes, species, returnarray = FALSE) ~ cbind(tapply(sizes,species,mean)),
                           model = "BM", tree = tree)


  pca1 <- stats::prcomp(geomorph::two.d.array(shapes))
  pca2 <- geomorph::gm.prcomp(geomorph::two.d.array(shapes))
  pca3 <- suppressWarnings(
    mvMORPH::mvgls.pca(mvMORPH::mvols(
      geomorph::two.d.array(shapes) ~ 1), plot = FALSE)
  )
  pca3$center <- colMeans(geomorph::two.d.array(shapes))

  bgpca1 <- bg_prcomp(x = geomorph::two.d.array(shapes), groups = species)
  bgpca1$x <- bgpca1$x * -1 ; bgpca1$rotation <- bgpca1$rotation * -1
  bgpca2 <- Morpho::groupPCA(dataarray = geomorph::two.d.array(shapes), groups = species)

  phypca1 <- phytools::phyl.pca(geomorph::two.d.array(expected_shapes(shapes, species)), tree = tree)
  phypca2 <- geomorph::gm.prcomp(geomorph::two.d.array(expected_shapes(shapes, species)), phy = tree, GLS = TRUE)
  phypca3 <- mvMORPH::mvgls.pca(mvMORPH::mvgls(
    expected_shapes(shapes, species, returnarray = FALSE) ~ 1,
    model = "BM", tree = tree), plot = FALSE)
  phypca3$center <- colMeans(mvMORPH::mvgls(
    expected_shapes(shapes, species, returnarray = FALSE) ~ 1,
    model = "BM", tree = tree)$fitted)

  pls1 <- pls_shapes(shapes = geomorph::two.d.array(shapes), X = sizes)
  pls2 <- Morpho::pls2B(y = geomorph::two.d.array(shapes), x = sizes)
  pls3 <- geomorph::two.b.pls(A2 = geomorph::two.d.array(shapes), A1 = sizes)


  ext_lm <- ax_transformation(obj = model1)
  ext_procD.lm <- ax_transformation(obj = model2)
  ext_mvols <- ax_transformation(obj = model3)
  ext_lm.rrpp <- ax_transformation(obj = model4)

  ext_procD.pgls <- ax_transformation(obj = model5)
  ext_mvgls <- ax_transformation(obj = model6)

  ext_prcomp <- ax_transformation(obj = pca1)
  ext_gm.prcomp1 <- ax_transformation(obj = pca2)
  ext_mvgls.pca1 <- ax_transformation(obj = pca3)

  ext_bg_prcomp <- ax_transformation(obj = bgpca1)
  ext_groupPCA <- ax_transformation(obj = bgpca2)

  ext_phyl.pca <- ax_transformation(obj = phypca1)
  ext_gm.prcomp2 <- ax_transformation(obj = phypca2)
  ext_mvgls.pca2 <- ax_transformation(obj = phypca3)

  ext_pls2b <- ax_transformation(obj = pls1)
  ext_pls2B <- ax_transformation(obj = pls2)
  ext_two.b.pls <- ax_transformation(obj = pls3)


  result1 <- all(round(arrayspecs(ext_lm, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_procD.lm, k = 2, p = 9)[,,1],5))
  result2 <- all(round(arrayspecs(ext_lm, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_mvols, k = 2, p = 9)[,,1],5))
  result3 <- all(round(arrayspecs(ext_lm, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_lm.rrpp, k = 2, p = 9)[,,1],5))

  result4 <- all(round(arrayspecs(ext_procD.pgls, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_mvgls, k = 2, p = 9)[,,1],5))

  result5 <- all(round(arrayspecs(ext_prcomp, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_gm.prcomp1, k = 2, p = 9)[,,1],5))
  result6 <- all(round(arrayspecs(ext_prcomp, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_mvgls.pca1, k = 2, p = 9)[,,1],5))

  result7 <- all(round(arrayspecs(ext_bg_prcomp, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_groupPCA, k = 2, p = 9)[,,1],5))

  result8 <- all(round(arrayspecs(ext_phyl.pca, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_gm.prcomp2, k = 2, p = 9)[,,1],5))
  result9 <- all(round(arrayspecs(ext_phyl.pca, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_mvgls.pca2, k = 2, p = 9)[,,1],5))

  result10 <- all(round(arrayspecs(ext_pls2b, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_pls2b, k = 2, p = 9)[,,1],5))
  result11 <- all(round(arrayspecs(ext_pls2b, k = 2, p = 9)[,,1],5) == round(arrayspecs(ext_two.b.pls, k = 2, p = 9)[,,1],5))

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,result8,result9,result10,result11))



  result1 <- all(round(arrayspecs(ext_lm, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_procD.lm, k = 2, p = 9)[,,2],5))
  result2 <- all(round(arrayspecs(ext_lm, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_mvols, k = 2, p = 9)[,,2],5))
  result3 <- all(round(arrayspecs(ext_lm, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_lm.rrpp, k = 2, p = 9)[,,2],5))

  result4 <- all(round(arrayspecs(ext_procD.pgls, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_mvgls, k = 2, p = 9)[,,2],5))

  result5 <- all(round(arrayspecs(ext_prcomp, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_gm.prcomp1, k = 2, p = 9)[,,2],5))
  result6 <- all(round(arrayspecs(ext_prcomp, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_mvgls.pca1, k = 2, p = 9)[,,2],5))

  result7 <- all(round(arrayspecs(ext_bg_prcomp, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_groupPCA, k = 2, p = 9)[,,2],5))

  result8 <- all(round(arrayspecs(ext_phyl.pca, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_gm.prcomp2, k = 2, p = 9)[,,2],5))
  result9 <- all(round(arrayspecs(ext_phyl.pca, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_mvgls.pca2, k = 2, p = 9)[,,2],5))

  result10 <- all(round(arrayspecs(ext_pls2b, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_pls2b, k = 2, p = 9)[,,2],5))
  result11 <- all(round(arrayspecs(ext_pls2b, k = 2, p = 9)[,,2],5) == round(arrayspecs(ext_two.b.pls, k = 2, p = 9)[,,2],5))

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,result8,result9,result10,result11))

})


test_that(desc = "addendum: testing ax_transformations for Momocs::PCA,
          dimensions and accuracyof shapes", code = {
  data(shells)
  shapes <- shells$shapes

  pca0 <- prcomp(shapes$coe)
  pca1 <- Momocs::PCA(shapes)

  ext0 <- ax_transformation(pca0)
  ext1 <- ax_transformation(pca1)

  result1 <- all(c(dim(ext0)[1], dim(ext1)[1]) == 2)
  result2 <- all(c(dim(ext0)[2], dim(ext1)[2]) == ncol(shapes$coe))

  result3 <- all(ext0 == ext1)

  expect_true(all(result1,result2,result3))

})


######################################################################################
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




