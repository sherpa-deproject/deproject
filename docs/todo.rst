
.. include:: references.rst

To Do
=====

 - Use the Python logging module to produce output and allow for a verbosity
   setting.  [Easy]
 - Move from using obsid to a more-generic label for the multi-data-set
   case (i.e. when there are multiple data sets in a given annulus)
 - Allow the theta value to vary between data set in each annulus (for the
   multi-data-set case)
 - Support elliptical annuli
 - Create and use more generalized ``ModelStack`` and ``DataStack`` classes
   to allow for general mixing models.  [Hard]
 - Add summary methods (similar to Sherpa's show series of commands)
 - Add a specialized version of `plot_source_component` which takes a
   general model expression (i.e. without the `_annulus` suffix) and
   displays all combos per annulus. A bit tricky to get right.
   There's also a `plot_model_component`.
