#################

test_that(desc = "testing rev_eigen, entire empiric data set", code = {
  data(tails)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))

  backshapes_mat <- rev_eigen(scores = pca$x,
                              vectors = pca$rotation, center = pca$center)
  backshapes_arr <- geomorph::arrayspecs(backshapes_mat, k = 2, p = 9)

  result <- all(round(tails$shapes, 10) == round(backshapes_arr, 10))
  expect_true(result)
})


test_that(desc = "testing rev_eigen, single empiric case", code = {
  data(tails)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))

  backshapes_mat <- rev_eigen(scores = pca$x[1,],
                              vectors = pca$rotation, center = pca$center)
  backshapes_arr <- matrix(backshapes_mat, ncol = 2, byrow = TRUE)

  result <- all(round(tails$shapes[,,1] ,10) == round(backshapes_arr ,10))
  expect_true(result)
})


test_that(desc = "testing rev_eigen (and consensus!), single theoretic case", code = {
  data(tails)
  meanshape <- consensus(tails$shapes)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))

  backshapes_mat <- rev_eigen(scores = rep(0, length.out = ncol(pca$x)),
                              vectors = pca$rotation, center = pca$center)
  backshapes_arr <- matrix(backshapes_mat, ncol = 2, byrow = TRUE)

  result <- all(round(meanshape, 10) == round(backshapes_arr, 10))
  expect_true(result)
})


test_that(desc = "testing rev_eigen (and consensus!), single vector", code = {
  data(tails)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))
  mshape <- consensus(tails$shapes)

  newshape_mat <- rev_eigen(scores = 0, vectors = pca$rotation[,1], center = pca$center)
  newshape_arr <- matrix (newshape_mat, ncol = 2, byrow = TRUE)

  result <- all(round(newshape_arr, 10) == round(mshape, 10))
  expect_true(result)
})



#################

test_that(desc = "testing proj_eigen, entire empiric data set", code = {
  data(tails)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))
  newscores <- proj_eigen(x = geomorph::two.d.array(tails$shapes),
                          vectors = pca$rotation, center = pca$center)

  result <- all(round(pca$x, 10) == round(newscores, 10))
  expect_true(result)
})


test_that(desc = "testing proj_eigen, single empiric case", code = {
  data(tails)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))
  newscore <- proj_eigen(x = geomorph::two.d.array(tails$shapes)[1,],
                         vectors = pca$rotation, center = pca$center)

  result <- all(round(pca$x[1,], 10) == round(newscore, 10))
  expect_true(result)
})


test_that(desc = "testing proj_eigen (and consensus!), single theoretic case", code = {
  data(tails)
  meanshape <- consensus(tails$shapes)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))

  centroid <- proj_eigen(x = matrix(t(meanshape), nrow = 1),
                         vectors = pca$rotation, center = pca$center)

  result <- all(round(centroid, 10) == round(colMeans(pca$x), 10))
  expect_true(result)
})


test_that(desc = "testing proj_eigen (and consensus!), single vector", code = {
  data(tails)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))
  mshape <- consensus(tails$shapes)

  newscore <- proj_eigen(x = matrix(t(mshape), nrow = 1, byrow = TRUE), vectors = pca$rotation[,1], center = pca$center)

  result <- round(newscore, 10) == 0
  expect_true(result)
})


#######################

test_that(desc = "testing inv_efourier, dimensions", code = {
  data(shells)
  shape_coe <- shells$shapes$coe[1,]

  n <- 50
  shape_xy <- inv_efourier(shape_coe, nb.pts = n)

  result1 <- nrow(shape_xy) == n
  result2 <- ncol(shape_xy) == 2

  expect_true(all(result1, result2))
})


######################

test_that(desc = "testing shapes_mat, fourier data", code = {
  data("shells")
  type <- shapes_mat(shells$shapes)$datype
  data2d <- shapes_mat(shells$shapes)$data2d

  result1 <- all(dim(data2d) == dim(shells$shapes))
  result2 <- type == "fcoef"


  data("shells")
  type <- shapes_mat(shells$shapes$coe)$datype
  data2d <- shapes_mat(shells$shapes$coe)$data2d

  result3 <- all(dim(data2d) == dim(shells$shapes))
  result4 <- type == "fcoef"

  expect_true(all(result1, result2, result3, result4))

})


test_that(desc = "testing shapes_mat, landmark data", code = {
  data("tails")
  type <- shapes_mat(tails$shapes)$datype
  data2d <- shapes_mat(tails$shapes)$data2d

  result1 <- all(dim(data2d) == c(dim(tails$shapes)[3],
                                  dim(tails$shapes)[1]*dim(tails$shapes)[2]))
  result2 <- type == "landm"

  tails_mat <- geomorph::two.d.array(tails$shapes)
  type <- shapes_mat(tails_mat)$datype
  data2d <- shapes_mat(tails_mat)$data2d

  result3 <- all(dim(data2d) == dim(tails_mat))
  result4 <- type == "landm"

  expect_true(all(result1, result2, result3, result4))
})


#################

test_that(desc = "testing morphogrid, dimensions", code = {
  data("tails")
  nv = 10
  nh = 20
  p = nrow(tails$shapes)
  k = ncol(tails$shapes)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))
  shapes_grid <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
                            k = k, p = p, nh = nh, nv = nv)

  result1 <- dim(shapes_grid$models_mat)[1] == nv * nh * p
  result2 <- all(dim(shapes_grid$models_arr)[1] == p,
                 dim(shapes_grid$models_arr)[2] == k,
                 dim(shapes_grid$models_arr)[3] == nh * nv)

  expect_true(all(result1, result2))
})


test_that(desc = "testing morphogrid, frame of the plot", code = {
  data("tails")
  nv = 2
  nh = 2
  p = nrow(tails$shapes)
  k = ncol(tails$shapes)
  pca <- prcomp(geomorph::two.d.array(tails$shapes))
  shapes_grid <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
                            k = k, p = p, nh = nh, nv = nv)

  result1 <- abs(diff(range(shapes_grid$models_mat[,1]))) >= abs(diff(range(pca$x[,1])))
  result2 <- abs(diff(range(shapes_grid$models_mat[,2]))) >= abs(diff(range(pca$x[,2])))
  expect_true(all(result1, result2))
})




