
# Notes

## General

We use conda for version control for command line programs and python pacakges.

## R

CURRENTLY THIS PIPELINE REQUIRES NO R PACKAGES. THE BELOW CODE SHOULD NOT BE RUN.

Conda, however, does not work well with R packages. One could use conda to install a base version of R. Alternatively, one can pull from from a rocker image (e.g., `rocker/r-ver:3.6.2` instead of `buildpack-deps:curl` or a minimal install like `alpine` ).

For R package management we use `renv`.

```R
# initialize renv
renv::init()

# update renv after intalling any packages
renv::settings$snapshot.type("simple")
renv::snapshot()

# add renv.lock, .Rprofile, and renv/activate.R to git
```

## Key R packages:

```
- r-optparse=1.6.6
- bioconductor-mast=1.14.0
- bioconductor-fgsea=1.14.0
- r-dirichletreg=0.7_0
```

## R package requirements: MAST

The lme4 package is needed to model random effects (de_method = "glmer").

```
- r-lme4=1.1_18_1
```


## R package requirements: RNA-enrich

```R
install.packages("annotate")
install.packages("stats")
install.packages("GO.db")
install.packages("org.Hs.eg.db")
install.packages("KEGG.db")
install.packages("mgcv")
install.packages("rms")

library("annotate")
library("stats")
library("GO.db")
library("org.Hs.eg.db")
library("KEGG.db")
library("mgcv")
library("rms")
```

```
- bioconductor-annotate=1.66.0
- bioconductor-go.db=3.11.4
- bioconductor-org.hs.eg.db=3.11.4
- bioconductor-kegg.db=3.2.4
- r-mgcv=1.8_31
- r-rms=6.0_1
```

Note: missing IDEA
