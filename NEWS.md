# gittins (development version)

Breaking changes:

* `bmab_gi_multiple()` (breaking changes):
  * Now accepts alpha/beta or Sigma/n state inputs.
  * Therefore `bmab_gi_multiple_ab()` is no longer needed and so has been removed.
  * Output is now a data frame with a row for each starting state. The original matrix output is 
  now attached as attributes.
* The ordering of some function arguments has been changed.

Other changes:

* The `tol` argument now has a default value (5e-4).
* Code tidying and style changes.
* Minor documentation edits. 
* Add GPL3 licence.

# gittins 0.1.0

This is the original working version from 2017. It is available on Github under releases/tags.
