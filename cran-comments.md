## Release Summary

This release fixes a *name* change in CompQuadForm and re-exports c.seqMeta.

## Test environments
* win-builder (devel and release)
* ubuntu 12.04 (on travis-ci) R 3.3.1
* Windows (on AppVeyor), 3.3.1
* local OSX install, R 3.3.1

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE on Win-Builder:

* checking CRAN incoming feasibility ... NOTE
  Possibly mis-spelled words in DESCRIPTION:
    SKAT (14:43)

SKAT is not mis-spelled it is a commonly used acronym in the genetics field.

## Downstream dependencies
There are no downstream dependencies.