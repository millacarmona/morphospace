
<!-- README.md is generated from README.Rmd. Please edit that file -->

# morphospace <img src="man/figures/morphosp_hex.png" align="right" width="200"/>

The goal of `morphospace` is to enhance representation and heuristic
exploration of multivariate ordinations of shape data. This package can
handle the most common types of shape data working in integration with
other widely used R packages such as `Morpho` (Schlager 2017),
`geomorph` (Adams et al. 2021), `shapes` (Dryden 2019), `Momocs`
(Bonhomme et al. 2014) and `mvMORPH` (Clavel et al.2015), which cover
other more essential steps in the geometric morphometrics pipeline
(e.g. importation, normalization, statistical analysis, phylogenetic
modeling).

Below there is broad-strokes account of the `morphospace` capacities;
for more specific guidance, refer to [General
usage](https://millacarmona.github.io/morphospace/articles/General-usage.html)
and [Worked
examples](https://millacarmona.github.io/morphospace/articles/Worked-examples.html).

## Installation

You can install the development and CRAN versions of `morphospace` from
[GitHub](https://github.com/) with:

``` r
# Development version
# install.packages("devtools")
devtools::install_github("millacarmona/morphospace")

# CRAN version
# install.packages("morphospace") #not yet ready
```

## Concept

The basic idea behind `morphospace` is to build empirical morphospaces
using multivariate ordination methods, then use the resulting ordination
as a reference frame in which elements representing different aspects of
morphometric variation are projected. These elements are added to both
graphic representations and objects as consecutive ‘layers’ and list
slots, respectively, using the `%>%` pipe operator from `magrittr`
(Bache & Wickham 2022).

The starting point of the `morphospace` workflow is a set of shapes
(i.e. morphometric data that is already free of variation due to
differences in orientation, position and scale). These are fed to the
`mspace` function, which generates a morphospace using a variety of
multivariate methods related to Principal Component Analysis. This
general workflow is outlined below using the `tails` data set from
Fasanelli et al. (2022), which contains tail shapes from 281 specimens
belonging to 13 species of the genus *Tyrannus*.

``` r
library(morphospace)
library(geomorph)
library(Morpho)
library(Momocs)
library(magrittr)
library(rgl)
```

``` r
# Load tail data
data("tails")

shapes <- tails$shapes
spp <- tails$data$species
wf <- tails$links
phy <- tails$tree

# Generate morphospace
mspace(shapes, links = wf, cex.ldm = 5)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

The ordination produced by `mspace` is used as a reference frame in
which scatter points, convex hulls / confidence ellipses, a phylogeny, a
set of morphometric axes or a landscape surface can be projected using
the `proj_*` functions:

``` r
# Get mean shapes of each species
spp_shapes <- expected_shapes(shapes = tails$shapes, x = tails$data$species)

# Generate morphospace and project:
msp <- mspace(shapes = shapes, links = wf, cex.ldm = 5) %>% 
  # scatter points
  proj_shapes(shapes = shapes, col = spp) %>% 
  # convex hulls enclosing groups
  proj_groups(shapes = shapes, groups = spp, alpha = 0.5) %>%
  # phylogenetic relationships
    proj_phylogeny(shapes = spp_shapes, tree = phy, lwd = 1.5, 
                 col.tips = match(phy$tip.label, levels(spp)))
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Once the `"mspace"` object has been created, the `plot_mspace` function
can be used to either regenerate/modify the plot, add a legend, or to
combine morphometric axes with other non-shape variables to produce
‘hybrid’ morphospaces. For example, PC1 can be plotted against size to
explore allometric patterns.

``` r
# Plot PC1 against log-size, add legend
plot_mspace(msp, x = tails$sizes, axes = 1, nh = 6, nv = 6, cex.ldm = 4, 
            alpha.groups = 0.5, col.points = spp, col.groups = 1:nlevels(spp), 
            phylo = FALSE, xlab = "Log-size", legend = TRUE)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Or against taxonomic classification to assess patterns of intra- and/or
interspecific variation.

``` r
# Plot PC1 against species classification
plot_mspace(msp, x = spp, axes = 1, nh = 6, nv = 6, cex.ldm = 4, boxplot.groups = TRUE,
            alpha.groups = 0.5, col.points = spp, col.groups = 1:nlevels(spp), 
            phylo = FALSE, xlab = "Log-size", legend = TRUE)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Alternatively, ordination axes can be combined with a phylogenetic tree
to create a phenogram:

``` r
# Plot vertical phenogram using PC1, add a legend
plot_mspace(msp, y = phy, axes = 1, nh = 6, nv = 6, cex.ldm = 4, 
            col.groups = 1:nlevels(spp), ylab = "Time", legend = TRUE)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

`morphospace` can also handle closed outlines (in the form of elliptic
Fourier coefficients) and 3D landmark data, as shown below briefly using
the `shells` and `shells3D` data sets:

``` r
# Load data
data("shells")

shapes <- shells$shapes
spp <- shells$data$species

# Generate morphospace
mspace(shapes, mag = 1, nh = 5, nv = 4, bg.model = "light gray") %>%
  proj_shapes(shapes = shapes, col = spp) %>%
  proj_groups(shapes = shapes, groups = spp, alpha = 0.5, ellipse = TRUE)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
# Load data
data("shells3D")

shapes <- shells3D$shapes
spp <- shells3D$data$species
mesh_meanspec <- shells3D$mesh_meanspec

# Generate surface mesh template
meanspec_shape <- shapes[,,findMeanSpec(shapes)]
meanmesh <- tps3d(x = mesh_meanspec, 
                  refmat = meanspec_shape, 
                  tarmat = expected_shapes(shapes))

# Generate morphospace
mspace(shapes, mag = 1, bg.model = "gray", cex.ldm = 0, template = meanmesh, 
       adj_frame = c(0.9, 0.85)) %>%
  proj_shapes(shapes = shapes, col = spp, pch = 16) %>%
  proj_groups(shapes = shapes, groups = spp, alpha = 0.3)
#> Preparing for snapshot: rotate mean shape to the desired orientation
#>  (don't close or minimize the rgl device).Press <Enter> in the console to continue:
#> 
#> This can take a few seconds...
#> DONE.
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Aside from working with these types of morphometric data, `morphospace`
provides functions to perform some useful shape operations, use TPS
interpolation of curves/meshes to improve visualizations, and supports a
variety of multivariate methods (bgPCA, phylogenetic PCA, PLS,
phylogenetic PLS) to produce ordinations. Integration with linear models
and ordinations generated with other packages (`geomorph`, `Morpho`,
`Momocs`, `mvMORPH`, `RRPP`, `phytools`) is also possible. For these and
other options and details, go to [General
usage](https://millacarmona.github.io/morphospace/articles/General-usage.html)
and [Worked
examples](https://millacarmona.github.io/morphospace/articles/Worked-examples.html).

## Update 1 (August 2022)

- Different behavior for `proj_shapes` (now replaces `mspace$x` with the
  actual scores being projected) and `proj_axis` (now adds one or more
  axes into an `mspace$shape_axis`).

- New `ellipses_by_groups_2D` (uses `car::ellipse`) function as an
  option for `proj_groups` and `plot_mspace`.

- Morphospaces without background shape models are now an option (for
  both `mspace` and `plot_mspace`).

- `plot_mspace` now regenerates the original mspace plot by default
  (`proj_*` functions were modified such that all the relevant graphical
  parameters are inherited downstream to `plot_mspace`), has further
  flexibility regarding hybrid morphospaces (`plot_phenogram` has been
  updated) and allows adding a legend (and some various bugs were fixed
  as well).

- Univariate morphospaces and associated density distributions are now
  an option (all the `mspace` workflow functions have been modified
  accordingly, especially `proj_shapes` and `proj_groups`).

- `consensus` and `expected_shapes` have been merged in a single
  function (the name `expected_shapes` was retained as the former was
  clashing with `ape::consensus`), which can handle both factors and
  numerics.

- Both `detrend_shapes` and `expected_shapes` can now calculate
  phylogenetically-corrected coefficients for interspecific data sets
  (Revell 2009).

## Update 2 (August 2023)

- The structure of `"mspace"` objects has been reorganized and now
  contain 3 main slots: `$ordination` (multivariate ordination details),
  `$projected` (elements added using `proj_*` functions) and `$plotinfo`
  (used for regeneration using `plot_mspace`). This has been
  complemented with a `print` method for the `"mspace"` class.

- New `proj_landscape` function has been added to represent adaptive
  surfaces interpolated from functional or performance indices (although
  can be used for any numerical variable).

- `proj_consensus` has been removed.

- New `extract_shapes` function for extracting synthetic shapes from
  `"mspace"` objects (background shape models, shapes along ordination
  axes, or specific coordinates selected interactively).

- New `burnaby` function, implementing Burnaby’s approach for
  standardization of morphometric data by computing a shape subspace
  orthogonal to an arbitrary vector or variable

- New `phyalign_comp` function, implementing Phylogenetically aligned
  component analysis, which finds the linear combination of variables
  maximizing covariation between trait variation and phylogenetic
  structure (Collyer & Adams 2021). Still a work in progress.

- Several internal adjustments have been introduced to the `mspace`,
  `proj_*` and `plot_mspace` functions in order to improve visualization
  and make the workflow more flexible.

- Legends created using `plot_mspace` have been improved, and scale bars
  for interpreting landscapes have also been made available.

## Update 3 (February 2024)

Significant changes aimed at enhancing integration with other GM/MV R
packages and improving procedures involving phylogenetic data, as well
as a couple of new features:

- The internal behavior of the `mspace` workflow has been modified so
  objects containing multivariate ordinations produced by `geomorph`,
  `Morpho`, `mvMORPH`, `Momocs` and `phytools` can be now used as input.

- `mspace` can now be fed with objects containing multivariate
  ordinations directly. This is implemented through an alternative
  combination of arguments for the `mspace` function (`ord` + `datype`
  as an alternative to `shapes` + `FUN` + `...`).

- `ax_transformation` (and by extension `proj_axis`) and
  `detrend_shapes` can now be fed with objects containing a linear model
  fitted using functions from `geomorph`, `RRPP` and `mvMORPH`. Internal
  phylogenetic correction of linear coefficients in `detrend_shapes` has
  been abandoned, relying now on the results of the (phylogenetic)
  linear model provided by the user.

- `phy_prcomp` and the experimental `phyalign_comp` have been removed
  (users interested in these methods should refer to
  `phytools::phyl.pca`, `geomorph::gm.prcomp` and/or
  `mvMORPH::mvgls.pca`).

- Estimation of ancestral shapes (performed internally by
  `expected_shapes`, `detrend_shapes`, `proj_phylogeny`, `pls2b`,
  `pls_shapes`, `burnaby` and `plot_mspace`) now relies on
  `mvMORPH::mvgls`, which has the advantage of allowing estimation under
  evolutionary models other than Brownian motion. In addition, ancestral
  character estimation of discrete non-shape variables (attempted
  internally by the phylogenetic version of `pls2b` (and by extension
  `pls_shapes`) and `plot_mspace` in certain situations) is performed
  under a simple Equal rates model via `ape::ace`.

- Introduction of boxplots and violin plots for combining shape
  ordination axes with categorical variables via `plot_mspace`.

- Tip and node labels can now be included in phylomorphospaces,
  phenograms and hybrid phylomorphospaces.

## Update 4 (February 2025)

A few relevant changes:

- `morphospace` can now generate Pareto rank ratio surfaces (Deakin et
  al. 2022), thanks to the code provided by Will Deakin.
  `proj_landscape` has been modified accordingly and can now take either
  two or more functions or two or more functional metrics (arguments
  `FUN` and `X`, respectively; see `?proj_landscape`).

- new `proj_pfront` function for projecting Pareto fronts (i.e., the
  subset of shapes optimizing two or more antagonistic measures of
  performance).

- `proj_landscape` can now accept surfaces created with `Morphoscape`
  (Dickson et al. 2023) through the argument `obj`.

- options for controlling the aspect ratio of morphospaces have been
  included (argument `asp` in `mspace`).

If you find any bugs please send me an email at `pablomillac@gmail.com`.
Thanks!!

## References

Adams D.C., Collyer M.L., Kaliontzopoulou A., & Baken E.K. (2021).
*geomorph: Software for geometric morphometric analyses*. R package
version 4.0.2. <https://cran.r-project.org/package=geomorph>.

Bache S.F., & Wickham H. (2022). *magrittr: A Forward-Pipe Operator for
R*. R package version 2.0.3.
<https://CRAN.R-project.org/package=magrittr>.

Bonhomme V., Picq S., Gaucherel C., & Claude J. (2014). *Momocs: Outline
Analysis Using R*. Journal of Statistical Software, 56(13), 1-24.
<https://www.jstatsoft.org/v56/i13/>.

Clavel, J., Escarguel, G., & Merceron, G. (2015). *mvMORPH: an R package
for fitting multivariate evolutionary models to morphometric data*.
Methods in Ecology and Evolution, 6(11), 1311-1319.
<https://doi.org/10.1111/2041-210X.12420>.

Collyer, M. L., & Adams, D. (2021). *Phylogenetically aligned component
analysis*. Methods in Ecology and Evolution, 12(2), 359-372.
<https://doi.org/10.1111/2041-210X.13515>.

Deakin, W. J., Anderson, P. S., den Boer, W., Smith, T. J., Hill, J. J.,
Rücklin, M., Donoghue, P. C. J. & Rayfield, E. J. (2022). *Increasing
morphological disparity and decreasing optimality for jaw speed and
strength during the radiation of jawed vertebrates*. Science Advances,
8(11), eabl3644. <https://doi.org/10.1126/sciadv.abl3644>.

Dickson, B. V., Pierce, S., & Greifer, N. 2023. *Morphoscape:
computation and visualization of adaptive landscapes*. R package version
1.0.2. <https://CRAN.R-project.org/package=Morphoscape>

Dryden, I.L. (2019). *shapes: statistical shape analysis*. R package
version 1.2.5. <https://CRAN.R-project.org/package=shapes>.

Fasanelli M.N., Milla Carmona P.S., Soto I.M., & Tuero, D.T. (2022).
*Allometry, sexual selection and evolutionary lines of least resistance
shaped the evolution of exaggerated sexual traits within the genus*
Tyrannus. Journal of Evolutionary Biology, in press.
<https://doi.org/10.1111/jeb.14000>.

Revell, L.J. (2009). *Size-correction and principal components for
interspecific comparative studies*. Evolution, 63, 3258-3268
<https://doi.org/10.1111/j.1558-5646.2009.00804.x>.

Schlager S. (2017). *Morpho and Rvcg - Shape Analysis in R*. In Zheng
G., Li S., Szekely G. (eds.), *Statistical Shape and Deformation
Analysis*, 217-256. Academic Press.
<https://doi.org/10.1016/B978-0-12-810493-4.00011-0>.
