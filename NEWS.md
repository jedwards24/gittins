# gittins 0.2.0

Breaking changes:

* `bmab_gi_multiple()` (breaking changes):
  * Now accepts alpha/beta or Sigma/n state inputs.
  * Therefore `bmab_gi_multiple_ab()` is no longer needed and so has been removed.
  * Output is now a data frame with a row for each starting state. The original matrix output is 
  now attached as attributes.
* `nmab_gi_multiple()` give output as a data frame.
* The ordering of some function arguments has been changed.

Other changes:

* Add checks that parameters take appropriate values. 
* The `tol` argument now has a default value (5e-4).
* Code tidying and style changes.
* Documentation edits. 
* Add GPL3 licence.
* Added `@keywords internal` to value functions and `calibrate_arm()` since they will not normally 
be used directly (they are still exported).
* Add package help file with useful links (see `?gittins`).
* Rewrite README.

# gittins 0.1.0

This is the original working version from 2017. It is available on Github under releases/tags.
