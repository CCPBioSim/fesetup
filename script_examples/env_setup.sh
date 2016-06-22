#  Copyright (C) 2012-2013,2016  Hannes H Loeffler
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  For full details of the license please see the COPYING file
#  that should have come with this distribution.

#
# required environment to run FESetup scripts
#
# dependencies: Python 2.7, Sire, Sire/Python2, OpenBabel + Python2, RDKit,
#               boost, gsl, qt4, eigen3, libz, xml2, numpy, py++0.95,
#               pygccxml0.95
#

FES_HOME=/home/hhl/projects/ccpbiosim/FESetup
#RDBASE=/usr/local/RDKit_2015_03_1
#DYLD_LIBRARY_PATH=$RDBASE/lib
#LD_LIBRARY_PATH=$RDBASE/lib

NAMDHOME=/usr/local/NAMD_2.11_Linux-x86-multicore
GMXHOME=/usr/local/gromacs-5.1.2
DLPOLYHOME=/usr/local/dl_poly_4.06
#AMBERHOME=

# first dir includes FESetup
PYTHONPATH=/home/hhl/projects/ccpbiosim/FESetup:$AMBERHOME/lib/python2.7/site-packages

export DYLD_LIBRARY_PATH LD_LIBRARY_PATH PYTHONPATH NAMDHOME GMXHOME DLPOLYHOME
export FES_HOME
