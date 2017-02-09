## Release Summary

This release fixes CRAN Windows old-release build error.
## Test environments
* win-builder (devel and release)
* Debian Linux, R-devel "Unsuffered Consequences" (on RHub)
* Windows old-release (on RHub), 3.2.5
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