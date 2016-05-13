#  Copyright (C) 2013-2016  Hannes H Loeffler
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
Create perturbed topologies.
"""


__revision__ = "$Id$"


import os

from FESetup import const, errors, logger
from FESetup.mutate import util

import amber


class PertTopology(object):

    def __init__(self, FE_sub_type, separate, ff, con_morph, atoms_initial,
                 atoms_final, lig_initial, lig_final, atom_map,
                 reverse_atom_map, zz_atoms, gaff):

        self.separate = separate
        self.ff = ff
        self.gaff = gaff
        self.con_morph = con_morph
        self.atoms_initial = atoms_initial
        self.atoms_final = atoms_final
        self.lig_initial = lig_initial
        self.lig_final = lig_final
        self.atom_map = atom_map
        self.reverse_atom_map = reverse_atom_map
        self.zz_atoms = zz_atoms

        self.files_created = []

        self.lig0 = None                # for GROMACS
        self.lig1 = None                # for GROMACS
        self.int_state = None           # for GROMACS
        self.parmtop = None             # for GROMACS
        self.inpcrd = None              # for GROMACS
        self.parmtop0 = None            # for GROMACS
        self.inpcrd0 = None             # for GROMACS
        self.parmtop1 = None            # for GROMACS
        self.inpcrd1 = None             # for GROMACS
        self.frcmod0 = None
        self.frcmod1 = None

        self.dummies0 = not all([a.atom for a in atom_map.keys()])
        self.dummies1 = not all([a.atom for a in atom_map.values()])

        self.mdin = True

        # special case: setup for Gromacs or CHARMM
        if FE_sub_type[0] == '_':
            self.FE_sub_type = FE_sub_type[1:]
            self.mdin = False
        else:
            want_softcore = FE_sub_type[:8] == 'softcore'
            self.FE_sub_type = ''

            # overwrite user choice, also for backward compatibility
            if self.separate and self.dummies0 and self.dummies1:
                if want_softcore:
                    self.FE_sub_type = 'softcore3'
                else:
                    self.FE_sub_type = 'dummy3'
            elif self.separate and (self.dummies0 or self.dummies1):
                if want_softcore:
                    self.FE_sub_type = 'softcore2'
                else:
                    self.FE_sub_type = 'dummy2'
            else:
                if want_softcore:
                    self.FE_sub_type = 'softcore'
                else:
                    self.FE_sub_type = 'dummy'

                if self.separate:
                    logger.write('Warning: linear transformation, not separating '
                                 'into vdw and electrostatic step\n')


    def setup(self, curr_dir, lig_morph, cmd1, cmd2):
        if self.FE_sub_type[:8] == 'softcore':
            state0, state1 = \
                    amber.softcore(lig_morph, self.lig_final,
                                   self.atom_map)

            pert0_info, pert1_info = None, None
        elif self.FE_sub_type[:5] == 'dummy':
            # also for Gromacs and CHARMM
            state0 = lig_morph
            state1, pert0_info, pert1_info = \
                    amber.dummy(lig_morph, self.con_morph,
                                self.lig_final, self.atom_map)
        else:
            raise NotImplementedError

        if self.mdin:
            amber.write_mdin(self.atoms_initial, self.atoms_final,
                             self.atom_map, 'sander', self.FE_sub_type, True)

        mol2_0 = os.path.join(curr_dir, const.MORPH_NAME + '0' +
                              const.MOL2_EXT)
        util.write_mol2(state0, mol2_0)

        frcmod0 = os.path.join(curr_dir, const.MORPH_NAME + '0.frcmod')

        lig0 = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_0,
                              start_fmt='mol2', frcmod=frcmod0,
                              gaff=self.gaff)

        lig0.set_atomtype(self.gaff)
        lig0._parmchk(mol2_0, 'mol2', frcmod0)
        lig0._parm_overwrite = 'state0'

        if pert0_info:
            lig0.prepare_top(pert=pert0_info)
        else:
            lig0.prepare_top()

        lig0.create_top(boxtype='', addcmd=cmd1 + cmd2)

        mol2_1 = os.path.join(curr_dir, const.MORPH_NAME + '1' +
                              const.MOL2_EXT)
        util.write_mol2(state1, mol2_1)

        frcmod1 = os.path.join(curr_dir, const.MORPH_NAME + '1.frcmod')

        lig1 = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_1,
                              start_fmt='mol2', frcmod=frcmod1,
                              gaff=self.gaff)

        lig1.set_atomtype(self.gaff)
        lig1._parmchk(mol2_1, 'mol2', frcmod1)
        lig1._parm_overwrite = 'state1'

        if pert1_info:
            lig1.prepare_top(pert=pert1_info)
        else:
            lig1.prepare_top()

        lig1.create_top(boxtype='', addcmd=cmd1 + cmd2)

        self.lig0 = lig0
        self.lig1 = lig1
        
        self.frcmod0 = frcmod0
        self.frcmod1 = frcmod1

        # NOTE: intermediate state assumed to not have dummies, so no
        #       missing parameters fixed through patching!
        if self.FE_sub_type == 'softcore2' or self.FE_sub_type == 'dummy2':
            ow_add = '_int'

            if self.FE_sub_type == 'dummy2':
                f = True
            else:
                f = False

            int_state = util.transfer_charges(state0, state1, self.atom_map, f)

            mol2_int = os.path.join(curr_dir, const.MORPH_NAME + ow_add +
                                    const.MOL2_EXT)
            util.write_mol2(int_state, mol2_int, resname = const.INT_NAME)

            lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_int,
                                 start_fmt='mol2', frcmod=frcmod0,
                                 gaff=self.gaff)
            lig.set_atomtype(self.gaff)
            lig._parm_overwrite = 'state_int'

            lig.prepare_top(add_frcmods=[frcmod1])
            lig.create_top(boxtype='', addcmd=cmd1 + cmd2)

        if self.FE_sub_type == 'dummy' or self.FE_sub_type == 'dummy2':
            top0 = lig0._parm_overwrite + lig0.TOP_EXT
            top1 = lig1._parm_overwrite + lig1.TOP_EXT

            util.patch_parmtop(top0, top1, ':%s' % const.LIGAND_NAME, '')

        if self.FE_sub_type == 'dummy3':
            ow_add = '_int'

            int_mol = util.zero_charges(state1, self.atom_map)

            mol2_int = os.path.join(curr_dir, const.MORPH_NAME + ow_add +
                                    const.MOL2_EXT)
            util.write_mol2(int_mol, mol2_int, resname = const.LIGAND_NAME)

            lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_int,
                                 start_fmt='mol2', frcmod=frcmod1,
                                 gaff=self.gaff)
            lig.set_atomtype(self.gaff)
            lig._parm_overwrite = 'state_int'

            lig.prepare_top(pert=pert1_info)
            lig.create_top(boxtype='')

            top0 = lig0._parm_overwrite + lig0.TOP_EXT
            int_name = lig._parm_overwrite + lig.TOP_EXT
            top1 = lig1._parm_overwrite + lig1.TOP_EXT

            util.patch_parmtop(top0, int_name, ':%s' % const.LIGAND_NAME, '')
            util.patch_parmtop(int_name, top1, ':%s' % const.LIGAND_NAME, '')

            self.int_state = lig


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2, boxdims):

        if self.FE_sub_type[:8] == 'softcore':
            state0, state1 = \
                    amber.softcore(lig_morph, self.lig_final,
                                   self.atom_map)

            pert0_info, pert1_info = None, None
        elif self.FE_sub_type[:5] == 'dummy':
            # also for Gromacs and CHARMM
            state0 = lig_morph
            state1, pert0_info, pert1_info = \
                    amber.dummy(lig_morph, self.con_morph,
                                self.lig_final, self.atom_map)

            leap_extra = 'source "%s"\n'
        else:
            raise NotImplementedError

        if self.mdin:
            amber.write_mdin(self.atoms_initial, self.atoms_final,
                             self.atom_map, 'sander', self.FE_sub_type, False)

        mol2_0 = os.path.join(curr_dir, const.MORPH_NAME + '0' +
                              const.MOL2_EXT)
        util.write_mol2(state0, mol2_0)

        com0 = self.ff.Complex(pdb_file, mol2_0)
        com0.box_dims = boxdims
        com0.ligand_fmt = 'mol2'
        com0.frcmod = self.frcmod0

        com0._parm_overwrite = 'state0'
        com0.ligand_fmt = 'mol2'

        if pert0_info:
            com0.prepare_top(gaff=self.gaff, pert=pert0_info)
        else:
            com0.prepare_top(gaff=self.gaff)

        com0.create_top(boxtype='set', addcmd=cmd1 + cmd2)

        mol2_1 = os.path.join(curr_dir, const.MORPH_NAME + '1' +
                              const.MOL2_EXT)
        util.write_mol2(state1, mol2_1)

        com1 = self.ff.Complex(pdb_file, mol2_1)
        com1.box_dims = boxdims
        com1.ligand_fmt = 'mol2'
        com1.frcmod = self.frcmod1

        com1._parm_overwrite = 'state1'

        if pert1_info:
            com1.prepare_top(gaff=self.gaff, pert=pert1_info)
        else:
            com1.prepare_top(gaff=self.gaff)

        com1.create_top(boxtype='set', addcmd=cmd1 + cmd2)

        if self.FE_sub_type == 'softcore2' or self.FE_sub_type == 'dummy2':
            ow_add = '_int'
            
            if self.FE_sub_type == 'dummy2':
                f = True
            else:
                f = False

            int_state = util.transfer_charges(state0, state1, self.atom_map, f)

            mol2_int = os.path.join(curr_dir, const.MORPH_NAME + ow_add +
                                    const.MOL2_EXT)
            util.write_mol2(int_state, mol2_int, resname = const.INT_NAME)

            com = self.ff.Complex(pdb_file, mol2_int)
            com.box_dims = boxdims
            com.ligand_fmt = 'mol2'
            com.frcmod = self.frcmod1
            com._parm_overwrite = 'state_int'
            com.prepare_top(add_frcmods=[self.frcmod0])
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)

        if self.FE_sub_type == 'dummy' or self.FE_sub_type == 'dummy2':
            top0 = com0._parm_overwrite + com0.TOP_EXT
            top1 = com1._parm_overwrite + com1.TOP_EXT

            util.patch_parmtop(top0, top1, ':%s' % const.LIGAND_NAME, '')
            
            self.parmtop = top0
            self.inpcrd = com0._parm_overwrite + com0.RST_EXT

        if self.FE_sub_type == 'dummy3':  # for GROMACS and CHARMM
            ow_add = '_int'

            int_mol = util.zero_charges(state1, self.atom_map)

            mol2_int = os.path.join(curr_dir, const.MORPH_NAME + ow_add +
                                    const.MOL2_EXT)
            util.write_mol2(int_mol, mol2_int, resname = const.LIGAND_NAME)

            com = self.ff.Complex(pdb_file, mol2_int)
            com.box_dims = boxdims
            com.ligand_fmt = 'mol2'
            com.frcmod = self.frcmod1
            com._parm_overwrite = 'state_int'

            if pert1_info:
                com.prepare_top(gaff=self.gaff, pert=pert1_info)
            else:
                com.prepare_top(gaff=self.gaff)

            com.create_top(boxtype='set')

            top0 = com0._parm_overwrite + com0.TOP_EXT
            int_name = com._parm_overwrite + com.TOP_EXT
            top1 = com1._parm_overwrite + com1.TOP_EXT

            util.patch_parmtop(top0, int_name, ':%s' % const.LIGAND_NAME, '')
            util.patch_parmtop(int_name, top1, ':%s' % const.LIGAND_NAME, '')

            self.parmtop0 = top0
            self.inpcrd0 = com0._parm_overwrite + com0.RST_EXT
            self.parmtop1 = top1
            self.inpcrd1 = com1._parm_overwrite + com1.RST_EXT
            self.int_state = com
