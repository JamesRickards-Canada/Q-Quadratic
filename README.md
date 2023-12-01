# Q Quadratic
A PARI/GP package for integral binary quadratic forms and quaternion algebras over Q, with an emphasis on indefinite quadratic forms and indefinite quaternion algebras.

At the moment, this package is compatible with versions 2.13.1 and later. If you are on Windows, then you must use Windows Subsystem for Linux (WSL).

This package is not regularly maintained, and with updates to PARI/GP some things may break. If they do, please let me know, and I will try to fix them.

## Installation
1. git clone this repository

2. You need to know where the version of PARI/GP you want to use is installed. The default location is inside /usr/local, but this may change based on your Linux distro, or if you want to use it through SageMath. If you think it is installed in the default location, you can simply call "make".

3. Otherwise, call "make setup" to search for the correct files. By default the program searches in "/usr", but there is a chance it is not installed there (this sometimes happens on a server). If this is the case, you can supply an alternate location.

4. If the program finds potential matches, it will ask you to confirm which files are correct, and saves them to "pari_loc.txt". Once this step completes, a call to "make" will compile the project! Modifying the program (e.g. via git pull) won't require redoing this setup, unless the version of PARI/GP or Sage you use changes.

5. Call "gp isogeny" to start gp and load the methods. ?qq accesses the help.

6. Call "make clean" to clean up the object files created.

## Documentation

The file "QQuadratic_GP_guide" is a guide to the GP-accessible methods. It may be slightly out of date.

The file "QQuadratic_PARI_guide" is a guide to using this package in library mode. It too may be slightly out of date.

## Papers

All relevant computations from my [thesis](https://math.colorado.edu/~jari2770/PDFs/thesis.pdf) were made with this package. The same can be said for the 3 papers that emanated from this thesis:
* [Computing intersections of closed geodesics on the modular curve](https://doi.org/10.1016/j.jnt.2020.11.024)
* [Counting intersection numbers of closed geodesics on Shimura curves](https://rdcu.be/c7DBo)
* [Hecke operators acting on optimal embeddings in indefinite quaternion algebras](https://doi.org/10.4064/aa210723-11-7)

If you want to test the examples of intersection series from "Hecke operators acting on optimal embeddings in indefinite quaternion algebras", use the file "intersectionseries.gp".

## Warnings

This was the first package I wrote in PARI, and there are some poor syntax/decision choices, code that is somewhat redundant with more up-to-date PARI/GP versions (when I began coding in gp, the version was 2.7.6), and code that can be significantly sped up. I do not plan to come back here and do those modifications. Instead, they are slowly working their way into other projects. For example, if you are interested in fundamental domains, then you should look at the infinitely improved (Fundamental-domains-for-Shimura-curves)[https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves]. If most of this package ends up improved elsewhere, then I might come back and update this repository (or make a new one).

