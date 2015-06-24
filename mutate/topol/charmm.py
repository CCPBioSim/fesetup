#  Copyright (C) 2015  Hannes H Loeffler
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
Create perturbed topologies for CHARMM/PERT.
"""


__revision__ = "$Id$"


import os, sys
from collections import defaultdict

import sander
from FESetup import const, errors, logger
from FESetup.prepare.amber import charmm



VAC_INP_FILE = 'vac.inp'
SOL_INP_FILE = 'sol.inp'
STATE_INT = 'state_int'

TR_TABLE = {'pert': 'dummy', 'pert2': 'dummy2'}


class PertTopology(object):

    def __init__(self, FE_sub_type, sc_type, ff, con_morph, atoms_initial,
                 atoms_final, lig_initial, lig_final, atom_map,
                 reverse_atom_map, zz_atoms):

        if FE_sub_type[:4] != 'pert':
            raise errors.SetupError('Morph code only supports the '
                                    'CHARMM/PERT module')

        self.sc_type = sc_type
        self.ff = ff
        self.con_morph = con_morph
        self.atoms_initial = atoms_initial
        self.atoms_final = atoms_final
        self.lig_initial = lig_initial
        self.lig_final = lig_final
        self.atom_map = atom_map
        self.reverse_atom_map = reverse_atom_map
        self.zz_atoms = zz_atoms
        self.stype = FE_sub_type

        self.topol = None
        

    def setup(self, curr_dir, lig_morph, cmd1, cmd2):

        topol = sander.PertTopology(TR_TABLE[self.stype], self.sc_type, self.ff,
                                    self.con_morph, self.atoms_initial,
                                    self.atoms_final, self.lig_initial,
                                    self.lig_final, self.atom_map,
                                    self.reverse_atom_map, self.zz_atoms)

        topol.setup(curr_dir, lig_morph, cmd1, cmd2)

        lig0 = topol.lig0._parm_overwrite
        lig1 = topol.lig1._parm_overwrite
        top = topol.lig0.TOP_EXT
        rst = topol.lig0.RST_EXT

        # top0 and top1 written mainly for debugging purposes
        top0 = charmm.CharmmTop()
        top0.readParm(lig0 + top, lig0 + rst)
        top0.writePsf('state0.psf')
        top0.writeCrd('state0.cor')

        if self.stype == 'pert2':
            top_int = charmm.CharmmTop()
            top_int.readParm(STATE_INT + top, STATE_INT + rst)
            top_int.writePsf(STATE_INT + '.psf')
            top_int.writeCrd(STATE_INT + '.cor')

        top1 = charmm.CharmmTop()
        top1.readParm(lig1 + top, lig1 + rst)
        top1.writePsf('state1.psf')
        top1.writeCrd('state1.cor')

        top0.combine(top1)              # adds parms from top1 to top0!
        top0.writeRtfPrm('combined.rtf', 'combined.prm')

        self.topol = topol

        with open(VAC_INP_FILE, 'w') as inp:
            inp.write(VAC_INP)


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2):
        """
        Create only topology file, not coordinates.
        """

        # FIXME: support FE_sub_type
        self.topol.create_coords(curr_dir, dir_name, lig_morph, pdb_file,
                                 system, cmd1, cmd2)

        lig0 = self.topol.lig0._parm_overwrite
        lig1 = self.topol.lig1._parm_overwrite
        top = self.topol.lig0.TOP_EXT
        rst = self.topol.lig0.RST_EXT

        # top0 and top1 written mainly for debugging purposes
        top0 = charmm.CharmmTop()
        top0.readParm(lig0 + top, lig0 + rst)
        top0.writePsf('state0.psf')
        top0.writeCrd('state0.cor')

        if self.stype == 'pert2':
            top_int = charmm.CharmmTop()
            top_int.readParm(STATE_INT + top, STATE_INT + rst)
            top_int.writePsf(STATE_INT + '.psf')
            top_int.writeCrd(STATE_INT + '.cor')

        top1 = charmm.CharmmTop()
        top1.readParm(lig1 + top, lig1 + rst)
        top1.writePsf('state1.psf')
        top1.writeCrd('state1.cor')

        top0.combine(top1)              # adds parms from top1 to top0!
        top0.writeRtfPrm('combined.rtf', 'combined.prm')

        with open(SOL_INP_FILE, 'w') as inp:
            inp.write(SOL_INP)



VAC_INP = '''* vacuum
* FESetup
*

'''

SOL_INP = '''* solution
* FESEtup
*

'''
