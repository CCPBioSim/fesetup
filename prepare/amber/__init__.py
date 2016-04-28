#  Copyright (C) 2012-2016  Hannes H Loeffler
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

"""
The amber package handles the AMBER type of force fields.  Currently, only
modern force fields (see Amber14 manual) are supported.
"""


__revision__ = "$Id$"



import os

from FESetup.prepare.amber.utils import *
from FESetup.prepare.amber.common import *
from FESetup.prepare.amber.ligand import *
from FESetup.prepare.amber.protein import *
from FESetup.prepare.amber.complex import *


AMBER_FF_TYPES = frozenset( ('protein.ff14SB',
                             'ff14SB', 'ff14ipq', 'ff12SB', 'ff10', 'ff99SB',
                             'ff99SBildn', 'ff99SBnmr', 'ff03.r1') )
AMBER_SOLVENT_TYPES = {
    'tip3p' : ['loadAmberParams frcmod.ionsjc_tip3p\n'
               'loadAmberParams frcmod.ionslrcm_%s_tip3p\n',
               'TIP3PBOX'],
    'tip4pew' : ['loadAmberParams frcmod.tip4pew\n'
                 'WAT = T4E\n'
                 'loadAmberParams frcmod.ionsjc_tip4pew\n'
                 'loadAmberParams frcmod.ionslrcm_%s_tip4pew\n',
                 'TIP4PEWBOX'],
    'spce' : ['loadAmberParams frcmod.spce\n'
              'WAT = SPC\n'
              'loadAmberParams frcmod.spce\n'
              'loadAmberParams frcmod.ionsjc_spce\n'
              'loadAmberParams frcmod.ionslrcm_%s_spce\n',
              'SPCBOX'],
    }

# NOTE: 12-6-4 type not supported yet because requires manipulation of
#       parmtop with parmed
AMBER_DIVALENT_TYPES = frozenset( ('hfe', 'cm', 'iod') )



def init(ff_type, solvent, div_ions, add_ons, mdengine, parmchk_version=2,
         gaff='gaff'):
    """
    Set force field and solvent types for *all* classes in the amber
    hierarchy.  Leap commands and names are directly written into the Common
    class *before* initialisation of the class.

    :param ff_type: the AMBER force field
    :type ff_type: string
    :param solvent: solvent force field
    :type solvent: string
    :param div_ions: divalent ion set (Li, Roberts, Chakravorty and Merz)
    :type div_ions: string
    :param add_ons: additional force fields
    :type add_ons: list
    :param mdengine: MD engine
    :type mdengine: string
    :param parmchk_version: version of parmchk, either 1 or 2
    :type parmchk_version: int
    :param gaff: GAFF version, either 'gaff' or 'gaff2'
    :type gaff: string
    """

    if ff_type in AMBER_FF_TYPES:
        ff_cmd = 'source leaprc.%s\n' % ff_type
        force_fields = [ff_type]

        for ff in add_ons:
            ff_cmd += 'source leaprc.%s\n' % ff
            force_fields.append(ff)
    else:
        print >> sys.stderr, 'Unsupported force field: %s' % ff_type
        print >> sys.stderr, 'Known types are: %s' % \
              ', '.join(AMBER_FF_TYPES)
        sys.exit(1)

    try:
        solvent_type = AMBER_SOLVENT_TYPES[solvent]
 
        # special treatment for namd (kludgy?)
        if '.namd' in mdengine.__module__ and solvent == 'tip4pew':
            solvent_type[0] += 'set default FlexibleWater on '

        solvent_load = solvent_type[0] % div_ions
        solvent_box = solvent_type[1]
    except KeyError, IndexError:
        print >> sys.stderr, 'Unsupported solvent type: %s' % solvent
        print >> sys.stderr, 'Known types are: %s' % \
              ', '.join(t for t in AMBER_SOLVENT_TYPES)
        sys.exit(1)

    # these are meant to be constants for the whole of the class hierarchy
    Common.ff_cmd = ff_cmd
    Common.ff_addons = add_ons
    Common.solvent = solvent
    Common.solvent_load = solvent_load
    Common.solvent_box = solvent_box
    Common.MDEngine = mdengine
    Common.parmchk_version = parmchk_version
    Common.gaff = gaff

    Common.force_fields = force_fields
