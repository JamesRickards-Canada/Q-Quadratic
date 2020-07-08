# Q Quadratic; current version: 0.2
A PARI/GP package for integral binary quadratic forms (and coming soon: quaternion algebras) over Q, with an emphasis on indefinite quadratic forms (and coming soon: indefinite quaternion algebras).
Currently only the binary quadratic form methods have been implemented; the quaternion algebra methods should be ready by October 2020.

If you would like to use this package in GP, please consult the "Qquadratic_GP_guide", and download the files "libqquadratic.so", and "qquadratic.gp".

If you would like to use this package with PARI in library mode, please consult the "Qquadratic_PARI_guide"

If you plan on modifying these methods, you need the files "c_base.c", "c_bqf.c", and "qquadraticdecl.h"

If you plan on using and not moddifying these methods, you instead need "libqquadratic.so", "qquadraticdecl.h", and "qquadratic.gp".

Note that the intersection number methods have not yet been added to the GP and PARI pdf guides, but you can access the help by calling ?bqf_int from GP.
