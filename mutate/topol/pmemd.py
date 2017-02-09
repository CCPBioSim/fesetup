#  Copyright (C) 2014-2016  Hannes H Loeffler
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
Create perturbed topologies for pmemd14 and later.
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

        self.frcmod0 = None
        self.frcmod1 = None

        self.dummies0 = not all([a.atom for a in atom_map.keys()])
        self.dummies1 = not all([a.atom for a in atom_map.values()])

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

        patch_parms = []

        if self.FE_sub_type[:8] == 'softcore':
            state0, state1 = amber.softcore(lig_morph, self.lig_final,
                                            self.atom_map)
            pert0_info, pert1_info = None, None
        elif self.FE_sub_type[:5] == 'dummy':
            state0 = lig_morph
            state1, pert0_info, pert1_info = \
                    amber.dummy(lig_morph, self.con_morph,
                                self.lig_final, self.atom_map)
        else:
            raise NotImplementedError

        amber.write_mdin(self.atoms_initial, self.atoms_final,
                         self.atom_map, 'pmemd', self.FE_sub_type, True)

        mol2_0 = os.path.join(curr_dir, const.MORPH_NAME + '0' +
                              const.MOL2_EXT)
        util.write_mol2(state0, mol2_0, resname = const.LIGAND0_NAME)

        frcmod0 = os.path.join(curr_dir, const.MORPH_NAME + '0.frcmod')

        mol2_1 = os.path.join(curr_dir, const.MORPH_NAME + '1' +
                              const.MOL2_EXT)
        util.write_mol2(state1, mol2_1, resname = const.LIGAND1_NAME)

        frcmod1 = os.path.join(curr_dir, const.MORPH_NAME + '1.frcmod')


        lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_0,
                             start_fmt='mol2', frcmod=frcmod0,
                             gaff=self.gaff)

        lig.set_atomtype(self.gaff)

        lig._parmchk(mol2_0, 'mol2', frcmod0)
        lig._parmchk(mol2_1, 'mol2', frcmod1)

        if self.FE_sub_type == 'softcore' or self.FE_sub_type == 'dummy':
            lig._parm_overwrite = 'onestep'
            self.files_created.extend(('onestep.parm7', 'onestep.rst7'))

        if self.FE_sub_type == 'softcore3' or self.FE_sub_type == 'dummy3':
            lig._parm_overwrite = 'vdw'
            self.files_created.extend(('vdw.parm7', 'vdw.rst7'))

        if self.FE_sub_type == 'softcore' or self.FE_sub_type == 'softcore3' \
               or self.FE_sub_type == 'dummy3':
            if pert0_info:
                pert0=pert0_info
            else:
                pert0=''

            if pert1_info:
                pert1=pert1_info
            else:
                pert1=''

            lig.prepare_top(pert=pert0)
            lig.leap.add_mol(mol2_1, 'mol2', [frcmod1], pert=pert1)
            lig.create_top(boxtype = '', addcmd = cmd1 + cmd2)

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

            lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_0,
                                 start_fmt='mol2', frcmod=frcmod0,
                                 gaff=self.gaff)
            lig.set_atomtype(self.gaff)

            if self.dummies0:
                lig._parm_overwrite = 'vdw'

                patch_parms.append( (lig._parm_overwrite,
                                     ':%s' % const.LIGAND0_NAME,
                                     ':%s' % const.INT_NAME) )
            else:
                lig._parm_overwrite = 'charge'

            lig.prepare_top(pert=pert0_info)
            # intermediate state does never have dummies
            lig.leap.add_mol(mol2_int, 'mol2', [frcmod1])
            lig.create_top(boxtype = '', addcmd = cmd1 + cmd2)

            lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_int,
                                 start_fmt='mol2', frcmod=frcmod1,
                                 gaff=self.gaff)
            lig.set_atomtype(self.gaff)

            if self.dummies1:
                lig._parm_overwrite = 'vdw'

                patch_parms.append( (lig._parm_overwrite,
                                     ':%s' % const.INT_NAME,
                                     ':%s' % const.LIGAND1_NAME) )
            else:
                lig._parm_overwrite = 'charge'

            # intermediate state does never have dummies
            lig.prepare_top()
            lig.leap.add_mol(mol2_1, 'mol2', [frcmod0], pert=pert1_info)
            lig.create_top(boxtype = '', addcmd = cmd1 + cmd2)

            self.files_created.extend(('charge.parm7', 'charge.rst7',
                                       'vdw.par7', 'vdw.rst7'))
        # FIXME: residue name will be both the same
        elif self.FE_sub_type == 'softcore3' or self.FE_sub_type == 'dummy3':
            lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_0,
                                 start_fmt='mol2', frcmod=frcmod0,
                                 gaff=self.gaff)
            lig.set_atomtype(self.gaff)
            lig._parm_overwrite = 'decharge'

            patch_parms.append( (lig._parm_overwrite,
                                 ':%s' % const.LIGAND0_NAME,
                                 ':%s' % const.LIGAND1_NAME) )

            lig.prepare_top(pert=pert0)
            lig.leap.add_mol(mol2_0, 'mol2', [frcmod0], pert=pert0)
            lig.create_top(boxtype = '', addcmd = cmd1 + cmd2)

            lig = self.ff.Ligand(const.MORPH_NAME, start_file=mol2_1,
                                 start_fmt='mol2', frcmod=frcmod1,
                                 gaff=self.gaff)
            lig.set_atomtype(self.gaff)
            lig._parm_overwrite = 'recharge'

            if pert1_info:
                pert=pert1_info
            else:
                pert=''

            lig.prepare_top(pert=pert1)
            lig.leap.add_mol(mol2_1, 'mol2', [frcmod1], pert=pert1)
            lig.create_top(boxtype = '', addcmd = cmd1 + cmd2)

            self.files_created.extend(('recharge.parm7', 'recharge.rst7',
                                       'decharge.parm7', 'decharge.rst7'))
        elif self.FE_sub_type == 'dummy':
            lig.prepare_top(pert=pert0_info)
            lig.leap.add_mol(mol2_1, 'mol2', [frcmod1], pert=pert1_info)
            lig.create_top(boxtype = '', addcmd = cmd1 + cmd2)

        self.frcmod0 = frcmod0
        self.frcmod1 = frcmod1

        if self.FE_sub_type[:5] == 'dummy':
            for prm in patch_parms:
                util.patch_parmtop(prm[0] + lig.TOP_EXT, "", prm[1], prm[2])


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2, boxdims):

        patch_parms = []

        if self.FE_sub_type[:8] == 'softcore':
            state0, state1 = \
                    amber.softcore(lig_morph, self.lig_final,
                                   self.atom_map)
            pert0_info, pert1_info = None, None
        elif self.FE_sub_type[:5] == 'dummy':
            state0 = lig_morph
            state1, pert0_info, pert1_info = \
                    amber.dummy(lig_morph, self.con_morph,
                                self.lig_final, self.atom_map)
        else:
            raise NotImplementedError

        amber.write_mdin(self.atoms_initial, self.atoms_final,
                         self.atom_map, 'pmemd', self.FE_sub_type, False)

        mol2_0 = os.path.join(curr_dir, const.MORPH_NAME + '0' +
                              const.MOL2_EXT)
        util.write_mol2(state0, mol2_0, resname = const.LIGAND0_NAME)

        mol2_1 = os.path.join(curr_dir, const.MORPH_NAME + '1' +
                              const.MOL2_EXT)
        util.write_mol2(state1, mol2_1, resname = const.LIGAND1_NAME)

        com = self.ff.Complex(pdb_file, mol2_0)
        com.__class__.SSBONDS_OFFSET = 2 # FIXME: kludge
        com.box_dims = boxdims
        com.ligand_fmt = 'mol2'
        com.frcmod = self.frcmod0

        if self.FE_sub_type == 'softcore' or self.FE_sub_type == 'dummy':
            com._parm_overwrite = 'onestep'

        if self.FE_sub_type == 'softcore3' or self.FE_sub_type == 'dummy3':
            com._parm_overwrite = 'vdw'

        if self.FE_sub_type == 'softcore' or self.FE_sub_type == 'softcore3' \
               or self.FE_sub_type == 'dummy3':
            if pert0_info:
                pert0=pert0_info
            else:
                pert0=''

            if pert1_info:
                pert1=pert1_info
            else:
                pert1=''

            com.prepare_top(gaff=self.gaff, pert=pert0)
            com.leap.add_mol(mol2_1, 'mol2', [self.frcmod1], pert=pert1)
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)

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

            com = self.ff.Complex(pdb_file, mol2_0)
            com.__class__.SSBONDS_OFFSET = 2 # FIXME: kludge
            com.box_dims = boxdims
            com.ligand_fmt = 'mol2'
            com.frcmod = self.frcmod0

            if self.dummies0:
                com._parm_overwrite = 'vdw'
                patch_parms.append( (com._parm_overwrite,
                                     ':%s' % const.LIGAND0_NAME,
                                     ':%s' % const.INT_NAME) )
            else:
                com._parm_overwrite = 'charge'

            com.prepare_top(gaff=self.gaff, pert=pert0_info)
            # intermediate state does never have dummies
            com.leap.add_mol(mol2_int, 'mol2', [self.frcmod1])
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)

            com = self.ff.Complex(pdb_file, mol2_int)
            com.__class__.SSBONDS_OFFSET = 2 # FIXME: kludge
            com.box_dims = boxdims
            com.ligand_fmt = 'mol2'
            com.frcmod = self.frcmod1

            if self.dummies1:
                com._parm_overwrite = 'vdw'
                patch_parms.append( (com._parm_overwrite,
                                     ':%s' % const.INT_NAME,
                                     ':%s' % const.LIGAND1_NAME) )
            else:
                com._parm_overwrite = 'charge'

            # intermediate state does never have dummies
            com.prepare_top(gaff=self.gaff)
            com.leap.add_mol(mol2_1, 'mol2', [self.frcmod0], pert=pert1_info)
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)

        # FIXME: residue name will be both the same
        elif self.FE_sub_type == 'softcore3' or self.FE_sub_type == 'dummy3':
            com = self.ff.Complex(pdb_file, mol2_0)
            com.__class__.SSBONDS_OFFSET = 2 # FIXME: kludge
            com.box_dims = boxdims
            com.ligand_fmt = 'mol2'
            com.frcmod = self.frcmod0
            com._parm_overwrite = 'decharge'

            com.prepare_top(gaff=self.gaff, pert=pert0)
            com.leap.add_mol(mol2_0, 'mol2', [self.frcmod0], pert=pert0)
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)

            com = self.ff.Complex(pdb_file, mol2_1)
            com.__class__.SSBONDS_OFFSET = 2 # FIXME: kludge
            com.box_dims = boxdims
            com.ligand_fmt = 'mol2'
            com.frcmod = self.frcmod1
            com._parm_overwrite = 'recharge'

            com.prepare_top(gaff=self.gaff, pert=pert1)
            com.leap.add_mol(mol2_1, 'mol2', [self.frcmod1], pert=pert1)
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)
        elif self.FE_sub_type == 'dummy':
            com.prepare_top(gaff=self.gaff, pert=pert0_info)
            com.leap.add_mol(mol2_1, 'mol2', [self.frcmod1], pert=pert1_info)
            com.create_top(boxtype='set', addcmd=cmd1 + cmd2)

        if self.FE_sub_type[:5] == 'dummy':
            for prm in patch_parms:
                util.patch_parmtop(prm[0] + com.TOP_EXT, "", prm[1], prm[2])
