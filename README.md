
<!-- README.md is generated from README.Rmd. Please edit that file -->

# morphospace

<!-- badges: start -->
<!-- badges: end -->

A package for people who can’t visualize morphospaces good and wanna
learn to do other geometric morphometrics stuff good too.

## Installation

You can install the development version of morphospace from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("millacarmona/morphospace")
```

This package is intended to aid exploration and depiction of
multivariate ordinations of shape data obtained through geometric
morphometrics analyses. The functions from `morphospace` have been
designed to work in an integrated workflow along other widely used
geometric morphometrics R packages such as `Morpho`, `geomorph`, and
`Momocs`, whose functions cover more essential steps of morphometric
analysis (e.g. import, normalization, and statistical analysis).

Below, the general concept and capabilities of `morphospace` are
displayed using two data sets representing different types of geometric
morphometric data. The first data set, taken from from Fasanelli et
al. (2022), contains a sample of tail shapes from the 13 species of the
genus *Tyrannus*, two of which (*T. savana* and *T. forficatus*) display
exaggeratedly tails, as well as a considerable allometric variation and
sexual dimorphism in tail shape. The `tails` data set contains landmark
data and centroid sizes from the tails of 281 specimens, the
classification fo each specimen to species and sex, and the phylogenetic
relationships between *Tyrannus* species (Harvey et al. 2020). Also
included are the links between landmarks to aid visualization.

``` r
library(morphospace)
library(geomorph)
#> Loading required package: RRPP
#> Loading required package: rgl
#> Loading required package: Matrix

# Load tail data and extract shapes, centroid sizes, classification of sex and species,
# links between landmarks, and phylogenetic tree

data("tails")
shapes <- tails$shapes
species <- tails$data$species
sizes <- tails$sizes
sex <- tails$data$sex
links <- tails$links
tree <- tails$tree
```

Morphometric variation is assumed to be already free of differences of
orientation, position and scale. This standardization can be readily
performed using functions from the aforementioned R packages.

<img src="man/figures/README-pressure-1.png" width="100%" />

## Shapes operations

This package provide some functions that perform basic operations with
shape variables, such as the calculation of mean shapes or the
analytical removal of undesired sources of variation.

``` r
# Remove variation associated to sexual dimorphism and compute the consensus shape
# of each species
cons_shapes <- consensus(shapes = shapes, index = species)
detr_shapes <- arrayspecs(detrend_shapes(model = lm(two.d.array(shapes) ~ sex)),
                          p = nrow(shapes), k = ncol(shapes))
detr_cons_shapes <- consensus(shapes = detr_shapes, index = species)
```

## Workflow

The besic idea behind `morphospace` is to build a reference ordination
by synthesizing morphometric variation with multivariate methods, and
then use the resulting

The basic idea behind the `morphospace` workflow is to build (empiric)
morphospaces using multivariate methods, then use the resulting
synthesis as a reference in which to project different elements. These
elements are added both to the plot and the `"mspace"` object as
succesive ‘layers’ or list slots, respectively, using the `%>%` pipe
operator from `magrittr`.

``` r
# Create and plot morphospace using detrended shapes, and project specimens
morphospace1 <- mspace(detr_shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
  proj_shapes(shapes = detr_shapes)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
names(morphospace1)
#> [1] "x"        "rotation" "center"   "datype"   "ordtype"  "plotinfo"

# Plot morphospace, project specimens and delimit species' range of variation
#using convex hulls
morphospace2 <- mspace(detr_shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
  proj_shapes(shapes = detr_shapes, col = species) %>%
  proj_groups(shapes = detr_shapes, groups = species)
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

``` r
names(morphospace2)
#> [1] "x"        "rotation" "center"   "datype"   "ordtype"  "plotinfo" "gr_class"

# Plot morphospace, project each species' mean shape and range of variation
morphospace3 <- mspace(detr_shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
  proj_consensus(shapes = detr_cons_shapes, pch = 21, bg = 1:13, cex = 1.2) %>%
  proj_groups(shapes = detr_shapes, groups = species)
```

<img src="man/figures/README-unnamed-chunk-4-3.png" width="100%" />

``` r
names(morphospace3)
#> [1] "x"            "rotation"     "center"       "datype"       "ordtype"     
#> [6] "plotinfo"     "gr_centroids" "gr_class"

# Plot morphospace, project species' meanshapes and phylogenetic structure
#(requires a phy object)
morphospace4 <- mspace(detr_shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
  proj_consensus(shapes = detr_cons_shapes, pch = 21, bg = 1:13, cex = 1.2) %>%
  proj_phylogeny(tree = tree, pch = 16)
```

<img src="man/figures/README-unnamed-chunk-4-4.png" width="100%" />

``` r
names(morphospace4)
#> [1] "x"            "rotation"     "center"       "datype"       "ordtype"     
#> [6] "plotinfo"     "gr_centroids" "phylo_scores" "phylo"
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
