Healpix
-----------------------

HEALPix is an acronym for Hierarchical Equal Area isoLatitude Pixelization of a sphere. As suggested in the name, this pixelization produces a subdivision of a spherical surface in which each pixel covers the same surface area as every other pixel. The figure below shows the partitioning of a sphere at progressively higher resolutions, from left to right. The green sphere represents the lowest resolution possible with the HEALPix base partitioning of the sphere surface into 12 equal sized pixels. The yellow sphere has a HEALPix grid of 48 pixels, the red sphere has 192 pixels, and the blue sphere has a grid of 768 pixels (~7.3 degree resolution).

![M45](http://healpix.jpl.nasa.gov/images/healpixGridRefinement.jpg)

This Healpix version greatly improves speed of original implementation. We were first to introduce compressed Range Sets into Healpix and reducing speed from N*N to N. This improvement was adopted by Gaia mission to map galaxy and was also ported to C++ code. 

Goal of this project is to be best Healpix implementation from point of code quality, simplicity and testability. We use RangeSets whenever possible (including triangle queries). We also promote more 'functional' programming style, using immutable data structures (sets, vectors..).

Links
-----------------------

This project is based on code from [Nikolay Kuropatkin](http://home.fnal.gov/~kuropat/HEALPIX/PixTools.html), which is Java port of original K.M. Gorski Fortran code.

There is  [official Source Forge site](http://sourceforge.net/projects/healpix/) with implementations in many languages.

There is also abadonded [Scala version](https://github.com/jankotek/asterope-scala-legacy/tree/master/src/org/asterope/healpix)) of this library. We decided it would be more usefull to continue this project in Java.