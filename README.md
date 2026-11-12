# MSstatsTMTAP personal fork

The GitHub repository is `MSstatsTMT-AP`; the installed R package is
`MSstatsTMTAP`. R package names cannot contain hyphens. The AP packages install
alongside the original packages without replacing them. Function names remain
unchanged; use explicit namespace calls when both versions are attached.

Install the renamed forks in dependency order from their feature branches
(after these changes have been committed and pushed):

```r
remotes::install_github("andrew2165/MSstats-AP@bait-normalization-release-3.23")
remotes::install_github("andrew2165/MSstatsTMT-AP@backport/bait-normalization-3.23")
remotes::install_github("andrew2165/MSstatsShiny-AP@feature/label-free-bait-normalization")

library(MSstatsAP)
library(MSstatsTMTAP)
library(MSstatsShinyAP)
MSstatsShinyAP::launch_MSstatsShiny()

# Check the actual installed package locations and GitHub source commits.
for (pkg in c("MSstatsAP", "MSstatsTMTAP", "MSstatsShinyAP")) {
  print(find.package(pkg))
  print(packageDescription(pkg)$RemoteSha)
}
```

Other MSstats family dependencies retain their original package names and may
also load the upstream MSstats namespace. AP analyses use the AP imports and
namespace calls. The upstream overview below is retained for background.
The linked upstream websites describe the original packages.

---

# MSstatsTMT

MSstatsTMT is an R-based package for detecting differentially abundant proteins in shotgun mass spectrometry-based proteomic experiments with tandem mass tag (TMT) labeling. It is applicable to isobaric labeling quantitative proteomics, including iTRAQ and TMT. The official webpage is http://msstats.org/msstatstmt/.
