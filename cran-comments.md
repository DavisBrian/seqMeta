## Release Summary

New maintainer for this package.  Arie (the old maintainer) is no longer active in this field.  A copy of our email exchange is at the bottom.  I have asked Arie to email CRAN@R-project.org to verify the change in maintainer. As the new maintainer I have read and agree to the CRAN policies.

This release fixes a couple edge cases which were handled incorrectly.  One new function was added.  The minimum R version was bumped to 3.1.0 so that we could take advantage of the 'anyNA' function.

## Test environments
* local Windows 7 install, R 3.2.1
* local OS X install, R 3.2.1
* ubuntu 12.04 (on travis-ci) R 3.2.1
* win-builder (on AppVeyor), R-devel

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

New maintainer:
  Brian Davis <Brian.Davis281@gmail.com>
  
Old maintainer(s):
  Arie Voorman <arie.voorman@gmail.com>

## Downstream dependencies
There are no downstream dependencies.