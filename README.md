# Q Quadratic
A PARI/GP package for integral binary quadratic forms and quaternion algebras over Q, with an emphasis on indefinite quadratic forms and indefinite quaternion algebras.

This package is not regularly maintained, and with updates to PARI/GP some things may break. If they do, please let me know, and I will try to fix them.

## Installation Instructions

### Downloading the code
Call ```git clone https://github.com/JamesRickards-Canada/Q-Quadratic.git```. If you are on Windows, be sure to git clone from WSL (see below), as Windows line endings (carriage returns) may be added to files, causing issues.

### Prerequisites
* PARI/GP, but _not_ the downloaded ready-to-go binary. The PARI/GP website has binaries for Windows and Mac avaliable, but these will not work with the package. See below for OS specific instructions.
* Version at least 2.13.1, though the more up-to-date the better.
* You should have a guess as to the location of the ```pari.cfg``` file for the version of PARI/GP you are running. Suggestions on how to do this can be found below.

### Operating systems
* **Linux** - No further requirements
* **Windows** - You need to use Windows Subsytem for Linux. See the [guide](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf) I wrote for additional instructions.
* **Mac** - You need to have [Homebrew](https://brew.sh/) installed. This is also an easy way to install PARI/GP: ```brew install pari```

### Where is pari.cfg?
* The configuration file will search for this, but it is preferrable to not search your entire hard drive (as this can be very slow). So, you should at least supply a guess as to the location of ```pari.cfg```. Often only the top-level folder (e.g. ```/usr``` or ```/opt```) suffices.
* On Linux or WSL, if you build PARI/GP from source, it should be located in ```/usr/local/lib/pari/pari.cfg```, or at least somewhere in the ```/usr``` folder.
* On a Mac, if you install PARI/GP with Homebrew, it may be found in a folder like ```/opt/homebrew/Cellar/pari/VERSION/lib/pari```. Searching ```/opt/homebrew``` should be fine.
* If you are obtaining it through SageMath, it might be found where the library files of SageMath are
* Assuming you open PARI/GP with the command ```gp```, try ```type -a gp```, which will display where this command lives. The corresponding file(s) are likely symbolic links, and you can call ```readlink -f LOCATION``` on each of them to see where it lives. This can provide a clue as to the place to search for ```pari.cfg```.
* Another clue comes from gp itself. Open gp, and type ```default()```. Look for the entries ```datadir``` and ```help```. It is _often_ the case that ```datadir``` is in ```X/share/pari```, ```help``` is in ```X/bin/gphelp```, and ```pari.cfg``` is in ```X/lib/pari/pari.cfg```.

### Configuring and building the package
* From inside the project folder, call ```./configure``` to initialize the project. This helps you search for ```pari.cfg```, and stores the location to a file. You should supply it with a folder to search in!
* The script displays the corresponding versions of the found files, so if you have multiple versions, you can choose the correct one. This can be useful if you keep multiple copies of PARI/GP around.
* If the location of the installation of PARI/GP does not change, you do not need to reconfigure. If when you update PARI/GP there is a new location (e.g. if the version number is in the file path of ```pari.cfg```), you should call ```./configure``` again.
* Call ```make``` to build the project, and ```make clean``` to remove all .o object files. If you update to a new version of PARI/GP, you must remake the project.
* Once this is done, a call to ```gp qquadratic``` starts gp with the package installed!
* Call ```?qq``` or consult the [User's Manual](QQuadratic_GP_guide.pdf) to access further help for this package.

## Documentation

The file [QQuadratic_GP_guide.pdf](QQuadratic_GP_guide.pdf) is a guide to the GP-accessible methods. It may be slightly out of date.

The file [QQuadratic_PARI_guide.pdf](QQuadratic_PARI_guide.pdf) is a guide to using this package in library mode. It too may be slightly out of date.

## Papers

All relevant computations from my [thesis](https://jamesrickards-canada.github.io/PDFs/thesis.pdf) were made with this package. The same can be said for the 3 papers that emanated from this thesis:
* [Computing intersections of closed geodesics on the modular curve](https://doi.org/10.1016/j.jnt.2020.11.024)
* [Counting intersection numbers of closed geodesics on Shimura curves](https://rdcu.be/c7DBo)
* [Hecke operators acting on optimal embeddings in indefinite quaternion algebras](https://doi.org/10.4064/aa210723-11-7)

If you want to test the examples of intersection series from "Hecke operators acting on optimal embeddings in indefinite quaternion algebras", use the file "intersectionseries.gp".

## Warnings

This was the first package I wrote in PARI, and there are some poor syntax/decision choices, code that is somewhat redundant with more up-to-date PARI/GP versions (when I began coding in gp, the version was 2.7.6), and code that can be significantly sped up. I do not plan to come back here and do those modifications. Instead, they are slowly working their way into other projects. For example, if you are interested in fundamental domains, then you should look at the infinitely improved [Fundamental-domains-for-Shimura-curves](https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves). If most of this package ends up improved elsewhere, then I might come back and update this repository (or make a new one).

