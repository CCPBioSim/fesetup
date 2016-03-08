#  Copyright (C) 2012-2014  Hannes H Loeffler
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

r"""
The ForceField class acts sort of a Strategy/Proxy: the actual force field and
morph clasess are swapped in dynamically and the actual Ligand, Protein and
Complex classes provided to the caller through the ForceField interface.
"""

__revision__ = "$Id$"


import sys, warnings

from FESetup import logger



class ForceField(object):
    """
    Swap in dynamically the right classes for a chosen force field.

    Expected to be present:
      Modules: forcefield, prepare, utils
      Classes: Ligand, Protein, Complex, Morph
      Functions: ff.utils.self_check(), ff.init()
    """

    def __init__(self, forcefield = 'amber', type = 'ff14SB',
                 solvent = 'tip3p', div_ions = 'hfe', add_ons = [],
                 mdengine = 'amber', parmchk_version = 2, gaff='gaff'):
        """
        Pull in the right classes for the chosen forcefield.

        :param forcefield: family of force fields
        :type forcefield: string
        :param type: the sub type within a force field family
        :type type: string
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

        try:
            ff = __import__(forcefield, globals = {'__name__': __name__})
            mde = __import__('mdengines', globals = {'__name__': __name__})
        except ImportError as detail:
            sys.exit('Error: %s\nEnsure all modules are properly installed' %
                     detail)
        except AttributeError as detail:
            sys.exit('Error: %s\nFailed to properly initialize %s' %
                     (detail, ff) )

        err_msg = ff.utils.self_check()

        if err_msg:
            sys.exit(err_msg)

        mdeng = mde.get_MDEngine(mdengine)
        ff.init(type, solvent, div_ions, add_ons, mdeng, parmchk_version, gaff)

        self.Ligand = ff.Ligand
        self.Protein = ff.Protein
        self.Complex = ff.Complex

        self.forcefield = forcefield
        self.type = type
        self.solvent = solvent
        self.div_ions = div_ions
        self.add_ons = add_ons
        self.mdengine = mdengine
        self.parmchk_version = parmchk_version


    def __repr__(self):
        return ('[forcefield=%s, type=%s, solvent=%s, div_ions=%s, '
                'add_ons=[%s], mdengine=%s, parmchk_version=%i]' %
                (self.forcefield, self.type, self.solvent, self.div_ions,
                 ','.join(self.add_ons), self.mdengine, self.parmchk_version) )
