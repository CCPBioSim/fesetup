============
Installation
============

Linux, Intel 32bit and 64bit

    Run the installer package from www.ccpbiosim.ac.uk/software in a convenient location (N will stand for the release you have downloaded).   You can run the installer with "--help" to see further options.  Here we describe interactive installation i.e. when run without any command line options.  The installer requires the xz compression tool to be installed on your system.

    	> cd /where/I/want/it		# replace the path to whatever you like
    	> ./FESetupN_Linux.sh           # extract all files into FESetupN/

    The installer will automatically detect which version to extract (either 32 or 64 bit). You will be asked to provide paths to AMBER, GROMACS, NAMD and DL_POLY.  It is strongly recommended to choose 'N' (the default, so just press Enter) when the first question suggests to use your existing $AMBERHOME.  Choose the internal path as suggested in the following question to avoid modifications to the original AmberTools installation.  Press Enter to accept defaults or to set an empty path if you do not have a certain MD package.  You will also be asked for a Python 2.7 interpreter (the default is "python", assumed to be located in the $PATH but read the note below).  You will find the FESetup script in ./FESetupN/FesetupMM/bin after successful installation (MM is either 32 or 64 depending on your hardware).
    Check that FESetup is working.  You can do this by running the test set from our first tutorial or just run the FESetup script without any command line parameters.  This will write the default input parameters to your terminal.
    You can copy/link the FESetup script to your PATH e.g. /usr/local/bin if you like, or create an alias to point to the script.

Our packages are self-contained and come with all relevant tools from AmberTools 16 including sander (pmemd still requires a full AMBER license).  To carry out standard MD simulations, in particular equilibration of your system, the abstract MD engine supports AMBER (both sander and pmemd), NAMD, GROMACS and DL_POLY.  Please note that currently we do not directly support NAMD's alchemical free energy methods though there is support for dual topology runs with AMBER inputs in NAMD (an additional PDB file is required to mark appearing/vanishing atoms, see NAMD manual).   Standard MD is supported for NAMD though.

Also note that you should use the standard Python interpreter (e.g. the one that comes as a package with your OS distribution or you download and compile from python.org).  Python versions that come with a package management systems of their own may break the assumptions that our installer makes with regards to shared libraries.  Specifially, anaconda appears to mess with the library search path and seems to disregard the setup in the FESetup script.

------------
Dependencies
------------

FESetup depends on various third-party software. All of these are included in the installer package. Here a list of dependencies for those who want to compile everything
themselves. Not listed are some dependencies which can be installed through the operating sytem's package management software. Some secondary dependencies are listed too.
Debian based systems have most libraries, toolkits and tools pre-compiled and ready to install through their package managment system.

Python 2.7

Sire/corelib 0.0.1, Sire/Python2 0.0.1: Qt4, Boost, GSL, BLAS/LAPACK, pcre3

Ambertools 16

OpenBabel 2.3.x: eigen, swig, xml2

RDKit 2016: numpy, Boost
