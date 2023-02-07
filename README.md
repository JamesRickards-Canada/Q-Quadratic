# Q Quadratic v1.0
A PARI/GP package for integral binary quadratic forms and quaternion algebras over Q, with an emphasis on indefinite quadratic forms and indefinite quaternion algebras.

At the moment, this package is compatible with versions 2.15 and 2.16. If you are on Windows, then you must use Windows Subsystem for Linux (WSL).
If your PARI/GP is stored in /usr/local, then you are all set! Call "gp qquadratic" to open gp and begin coding.
If this is not the case, then you must remake the package. Modify the first line of the Makefile to the correct installation location, and call make to remake it.

This package is developed with version 2.16, so the 2.15 library file will not necessarily be up to date (though I don't expect to work on this much anymore anyway). If in doubt, remake the package.

The file "QQuadratic_GP_guide" is a guide to the GP-accessible methods. It may be slightly out of date.
The file "QQuadratic_PARI_guide" is a guide to using this package in library mode. It too may be slightly out of date.

If you want to test the examples of intersection series from my paper "Hecke operators acting on optimal embeddings in indefinite quaternion algebras", use the file "intersectionseries.gp".
