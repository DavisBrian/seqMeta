## Release Summary

This release fixes a couple edge cases which were handled incorrectly.  The minimum 'coxme' package version was bumped to 2.2-4 as it changed the definition of variance component estimates in lmekin. Minor changes to user documatation and examples were also made.

## Test environments
* win-builder (devel and release)
* ubuntu 12.04 (on travis-ci) R 3.2.2, 3.2.3, 
* OSX (on travis-ci) R 3.2.2, 3.2.3, 
* Windows (on AppVeyor), R-devel
* local Windows 7 install, R 3.2.2

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

New maintainer:
  Brian Davis <Brian.Davis281@gmail.com>
  
Old maintainer(s):
  Arie Voorman <arie.voorman@gmail.com>
  
On Win-builder along with this NOTE was 


Possibly mis-spelled words in DESCRIPTION:
  SKAT (14:43)
  
SKAT is not mis-spelled it is a comonly used acronym in the genetics field.

## Downstream dependencies
There are no downstream dependencies.