
################################################################################

#' \emph{Drosophila} wing shape data set
#'
#' @description Sample of intra- and interspecific wing shapes, quantified
#'   using 9 landmarks, from two species of cactophilic \emph{Drosophila}
#'   (Diptera: Drosophilidae).
#'
#' @format A list containing:
#' \describe{
#'   \item{$shapes:}{ a \code{9 x 2 x 263} array with the superimposed configurations
#'   of landmarks marking veins intersections on the right wing of each specimen.}
#'   \item{$sizes:}{ a vector of length 263 containing the centroid size of
#'   each specimen's wing.}
#'   \item{$data: }{ a 3-column data frame with information about the taxonomic
#'   classification, host cactus, sex, line and replica of each specimen.}
#'   \item{$links: }{ a list indicating how to connect the landmarks to create a
#'   wireframe.}
#'   \item{$template: }{ a 2-column matrix containing the coordinates of the
#'   9 landmarks, followed by the coordinates defining 9 curves curve describing
#'   wing morphology, to be warped using TPS interpolation.}
#'   }
#'
#' @details This data set includes wing shape data from 263 experimentally bred specimens
#'   belonging to two recently diverged cactophilic *Drosophila* species, *D. koepferae*
#'   and *D. buzzattii* from northwestern Argentina. These two species show preference to
#'   different host cacti (columnar cacti from the genus *Trichocereus* and prickly pears
#'   from the genus *Opuntia*, respectively) and are considered morphologically cryptic
#'   (usually, they can be recognized only from genital morphology). Accompanying data
#'   include the centroid size of each wing, as well as information on the host cacti each
#'   specimen was bred, its species, sex, isofemale line, and replica. Also included is a
#'   template containing a series of curves describing wing outline and veins (created using
#'   `build_template2d`), which will be warped using TPS interpolation.
#'
#' @references
#'   Soto, E. M., Goenaga, J., Hurtado, J. P., & Hasson, E. (2012). \emph{Oviposition
#'   and performance in natural hosts in cactophilic} Drosophila. Evolutionary Ecology,
#'   26(4), 975-990.
"wings"


################################################################################

#' \emph{Ptychomya} shell outlines data set
#'
#' @description Sample of intra- and interspecific shell shapes, quantified
#'   as closed outlines using 7 harmonics and elliptic Fourier analysis, from
#'   4 species of the extinct genus \emph{Ptychomya} (Bivalvia: Crassatellidae).
#'
#' @format A list containing:
#' \describe{
#'   \item{$shapes:}{ an \code{"OutCoe"} object containing a \code{137 x 28}
#'   matrix with the normalized Fourier coefficients describing the outline
#'   of the shell of each specimen in left lateral view.}
#'   \item{$sizes:}{ a vector of length 137 containing the centroid size of each
#'   fossil shell.}
#'   \item{$data: }{ a 4-column data frame with information about the taxonomic
#'   classification, ammonoid biozone, estimated geochronologic age (in million
#'   years), and geographic provenance of each fossil shell.}
#'   }
#'
#' @details The \code{shells} data contain data from 137 specimens belonging to
#'   4 species of the extinct bivalve genus \emph{Ptychomya}, tracking their
#'   morphological changes through a 5 million years interval from the Lower
#'   Cretaceous of west-central Argentina (approximately 140 million years ago).
#'   The data set includes the information about the external shell shape (measured
#'   using Fourier coefficients from 7 harmonics), centroid size, age (both relative,
#'   taken from ammonoid biozones, and absolute, estimated using a combination
#'   of the former, absolute dates, and stratigraphic data), geographic provenance,
#'   and taxonomic classification of each fossil specimen.
#'
#' @references
#'   Milla Carmona, P. S., Lazo, D. G., & Soto, I. M. (2018). \emph{Morphological
#'   evolution of the bivalve} Ptychomya \emph{through the Lower Cretaceous of
#'   Argentina}. Paleobiology, 44(1), 101-117.
"shells"


################################################################################

#' \emph{Steinmanella} shell surfaces data set
#'
#' @description Sample of ontogenetic longitudinal shell shapes, quantified using
#'   90 surface semilandmarks, from 7 species of the extinct genus \emph{Steinmanella}
#'   (Bivalvia: Trigoniida).
#'
#' @format A list containing:
#' \describe{
#'   \item{$shapes:}{ a \code{90 x 3 x 278} array with the slid, superimposed
#'   configurations of semilandmarks describing surfaces representing right
#'   valves (3-5 per specimen).}
#'   \item{$sizes:}{ a vector of length 278 containing the centroid size of
#'   each shell surface}
#'   \item{$data: }{ a 4-column data frame with information about the individual
#'   id, taxonomic classification, ammonoid biozone and geographic provenance
#'   of each shell surface.}
#'   \item{$mesh_meanspec: }{ a 3D surface mesh stored as a \code{"mesh3d"} object
#'   (from \code{rgl}) corresponding to the shell shape closest to the sample's
#'   mean (found using [geomorph::findMeanSpec()]), to be warped using TPS
#'   interpolation.}
#'   }
#'
#' @details The \code{shells3D} data contain data from 67 specimens belonging to
#'   7 species of the extinct bivalve genus \emph{Steinmanella}, which inhabited
#'   the Early Cretaceous of west-central Argentina during a ~7 million years interval.
#'   Shapes were measured in each individual using growth lines marking the form of the
#'   animal at different ages, making a total of 278 shell surfaces describing
#'   intra- and interspecific ontogenetic variation. This data includes the
#'   semilandmarks describing these surfaces' shapes (digitally processed to maximize
#'   signal), their centroid sizes, the individual to which they belong, as well as
#'   their relative age (ammonoid biozones), geographic provenance, and taxonomic
#'   classification. In order to improve visualization, the surface mesh of the
#'   specimen closest to the consensus of the sample, to be warped using TPS
#'   interpolation, has also been included.
#'
#' @references
#'   Milla Carmona, P. S., Lazo, D. G., & Soto, I. M. (2021). \emph{Ontogeny in
#'   the steinmanellines (Bivalvia: Trigoniida): an intra- and interspecific
#'   appraisal using the Early Cretaceous faunas from the Neuquen Basin as a
#'   case study}. Paleobiology, in press.
"shells3D"


################################################################################

#' \emph{Tyrannus} tail shape data set
#'
#' @description Sample of intra- and interspecific tail shapes, quantified
#'   using 9 landmarks, from the 13 species of the genus \emph{Tyrannus}
#'   (Aves: Tyrannidae).
#'
#' @format A list containing:
#' \describe{
#'   \item{$shapes:}{ a \code{9 x 2 x 281} array with the symmetric component
#'   of the superimposed configurations of landmarks marking the base of the tail,
#'   as well as the tips of the two outermost rectrices of each specimen.}
#'   \item{$sizes:}{ a vector of length 281 containing the centroid size of
#'   each specimen's tail.}
#'   \item{$data: }{ a 3-column data frame with information about the taxonomic
#'   classification, sex and type (deep forked, DF or non-deep forked, NDF.}
#'   \item{$tree: }{ a \code{"phy"} object containing the phylogenetic tree
#'   of the 13 \emph{Tyrannus} species.}
#'   \item{$links: }{ a list indicating how to connect the landmarks to create a
#'   wireframe.}
#'   }
#'
#' @details This data set contains a sample of tail shapes from the 13 species
#'   of the genus \emph{Tyrannus}, two of which (\emph{T. savana} and
#'   \emph{T. forficatus}) display exaggeratedly elongated tails, as well as a
#'   considerable allometric variation and sexual dimorphism. The \code{tails}
#'   data set  contains landmark data and centroid sizes from the tails of 281
#'   specimens, their classification to species and sex, and the phylogenetic
#'   relationships between \emph{Tyrannus} species (see Fasanelli et al. 2022
#'   and references therein). To further help visualization of shapes, the links
#'   between landmarks have also been included.
#'
#'   Have in mind that since the asymmetric component has been removed from
#'   shape variation, this half the rank of a typical shape data matrix.
#'
#' @references
#'   Fasanelli M. N., Milla Carmona P. S., Soto I. M., & Tuero, D.T . (2022).
#'   \emph{Allometry, sexual selection and evolutionary lines of least
#'   resistance shaped the evolution of exaggerated sexual traits within the
#'   genus} Tyrannus. Journal of Evolutionary Biology, 35(5), 669 - 679.
"tails"
