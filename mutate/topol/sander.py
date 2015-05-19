#  Copyright (C) 2013-2014  Hannes H Loeffler
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

import Sire.Mol
import Sire.Units

from FESetup import const, errors, logger
from FESetup.mutate import util



class PertTopology(object):

    def __init__(self, FE_sub_type, sc_type, ff, con_morph, atoms_initial,
                 atoms_final, lig_initial, lig_final, atom_map,
                 reverse_atom_map, zz_atoms):

        self.FE_sub_type = FE_sub_type
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

        self.lig0 = None                # for GROMACS
        self.lig1 = None                # for GROMACS
        self.parmtop = None             # for GROMACS
        self.inpcrd = None              # for GROMACS
        self.frcmod0 = None
        self.frcmod1 = None


    def setup(self, curr_dir, lig_morph, cmd1, cmd2):

        if self.FE_sub_type == 'softcore':
            util.amber_input(self.atoms_initial, self.atoms_final,
                             self.atom_map, self.sc_type, '', True)

            state0, state1 = \
                    util.amber_softcore(lig_morph, self.lig_final,
                                        self.atom_map)

            leap_extra = ''
            ow_add = '_sc'
        elif self.FE_sub_type == 'dummy':
            state0 = lig_morph
            state1 = amber_dummy(lig_morph, self.con_morph,
                                 self.lig_final, self.atom_map)

            leap_extra = 'source "%s"\n'
            ow_add = '_dummy'
        else:
            raise NotImplementedError

        mol2_0 = os.path.join(curr_dir, const.MORPH_NAME + ow_add + '0' +
                              const.MOL2_EXT)
        util.write_mol2(state0, mol2_0)

        frcmod0 = os.path.join(curr_dir, const.MORPH_NAME + ow_add +
                               '0.frcmod')

        lig0 = self.ff.Ligand(const.MORPH_NAME, '', start_file = mol2_0,
                              start_fmt = 'mol2', frcmod = frcmod0)

        lig0.set_atomtype('gaff')

        lig0._parmchk(mol2_0, 'mol2', frcmod0)

        lig0._parm_overwrite = 'state0' + ow_add

        if leap_extra:
            lig0.create_top(boxtype = '', addcmd = cmd1 + cmd2,
                            addcmd2 = leap_extra %
                            os.path.join(curr_dir, const.LEAP_PERT0_FILE) )
        else:
            lig0.create_top(boxtype = '', addcmd = cmd1 + cmd2)

        mol2_1 = os.path.join(curr_dir, const.MORPH_NAME + ow_add + '1' +
                              const.MOL2_EXT)
        util.write_mol2(state1, mol2_1)

        frcmod1 = os.path.join(curr_dir, const.MORPH_NAME + ow_add +
                               '1.frcmod')

        lig1 = self.ff.Ligand(const.MORPH_NAME, '', start_file = mol2_1,
                              start_fmt = 'mol2', frcmod = frcmod1)

        lig1.set_atomtype('gaff')

        lig1._parmchk(mol2_1, 'mol2', frcmod1)

        lig1._parm_overwrite = 'state1' + ow_add

        if leap_extra:
            lig1.create_top(boxtype = '', addcmd = cmd1 + cmd2,
                            addcmd2 = leap_extra %
                            os.path.join(curr_dir, const.LEAP_PERT1_FILE) )
        else:
            lig1.create_top(boxtype = '', addcmd = cmd1 + cmd2)

        self.lig0 = lig0
        self.lig1 = lig1
        
        self.frcmod0 = frcmod0
        self.frcmod1 = frcmod1

        if self.FE_sub_type == 'dummy':
            top0 = lig0._parm_overwrite + lig0.TOP_EXT
            top1 = lig1._parm_overwrite + lig1.TOP_EXT

            util.patch_parmtop(top0, top1)


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2):

        if self.FE_sub_type == 'softcore':
            util.amber_input(self.atoms_initial, self.atoms_final,
                             self.atom_map, self.sc_type, '', False)

            state0, state1 = \
                    util.amber_softcore(lig_morph, self.lig_final,
                                        self.atom_map)

            leap_extra = ''
            ow_add = '_sc'
        elif self.FE_sub_type == 'dummy':
            state0 = lig_morph
            state1 = amber_dummy(lig_morph, self.con_morph,
                                 self.lig_final, self.atom_map)

            leap_extra = 'source "%s"\n'
            ow_add = '_dummy'
        else:
            raise NotImplementedError

        mol2_0 = os.path.join(curr_dir, const.MORPH_NAME + ow_add + '0' +
                              const.MOL2_EXT)
        util.write_mol2(state0, mol2_0)

        com0 = self.ff.Complex(pdb_file, mol2_0)
        com0.ligand_fmt = 'mol2'
        com0.frcmod = self.frcmod0

        com0._parm_overwrite = 'state0' + ow_add
        com0.ligand_fmt = 'mol2'

        if leap_extra:
            com0.create_top(boxtype = 'set',
                            addcmd = cmd1 + cmd2,
                            addcmd2 = leap_extra %
                            os.path.join(curr_dir, const.LEAP_PERT0_FILE) )
        else:
            com0.create_top(boxtype = 'set',
                            addcmd = cmd1 + cmd2)

        mol2_1 = os.path.join(curr_dir, const.MORPH_NAME + ow_add + '1' +
                              const.MOL2_EXT)
        util.write_mol2(state1, mol2_1)

        com1 = self.ff.Complex(pdb_file, mol2_1)
        com1.ligand_fmt = 'mol2'
        com1.frcmod = self.frcmod1

        com1._parm_overwrite = 'state1' + ow_add

        if leap_extra:
            com1.create_top(boxtype = 'set',
                            addcmd = cmd1 + cmd2,
                            addcmd2 = leap_extra %
                            os.path.join(curr_dir, const.LEAP_PERT1_FILE) )
        else:
            com1.create_top(boxtype = 'set',
                            addcmd = cmd1 + cmd2)

        if self.FE_sub_type == 'dummy':
            top0 = com0._parm_overwrite + com0.TOP_EXT
            top1 = com1._parm_overwrite + com1.TOP_EXT

            util.patch_parmtop(top0, top1)
            
            self.parmtop = top0
            self.inpcrd = com0._parm_overwrite + com0.RST_EXT


def amber_dummy(lig_morph, con_morph, lig_final, atom_map):
    """
    Support for AMBER/dummy which creates leap commands to deal with its
    valency check.  Also, create the final state from the morph.

    :param lig_morph: the morph molecule
    :type lig_morph: Sire.Mol.CutGroup
    :param con_morph: the connectivity of the morph
    :type con_morph: Sire.Mol.Connectivity
    :param lig_final: the final state molecule
    :type lig_final: Sire.Mol.Molecule
    :param atom_map: the forward atom map
    :type atom_map: dict of _AtomInfo to _AtomInfo
    :returns: initial state molecule, final state molecule
    :rtype: Sire.Mol.Molecule, Sire.Mol.Molecule
    """

    state1_m = Sire.Mol.Molecule(lig_morph)
    state1 = state1_m.edit()    # MolEditor

    leap_pert0 = []
    leap_pert1 = []


    for iinfo, finfo in atom_map.items():
        istr = iinfo.name.value()
        fstr = finfo.name.value()

        new = state1.atom(iinfo.index)  # AtomEditor

        # Dummy may be bonded to H which exceeds leap's valency limit.  Leap
        # will only keep the first bond encountered.  Solution:
        # "pert=true" and explicit bonding to DU, recreate all bonds to not
        # rely on which one leap keeps
        if not iinfo.atom:
            for idx in con_morph.connectionsTo(iinfo.index):
                atom = lig_morph.select(idx)
                #resnum = atom.residue().number().value()
                name = '%s' % atom.name().value()
                ambertype = '%s' % atom.property('ambertype')

                # FIXME: hard-coded unit name
                if ambertype.upper().startswith('H') and \
                       name.startswith('H'):
                    ur = 's.1' #% resnum
                    leap_pert0.append('set %s.%s pert true' % (ur, name) )
                    leap_pert0.append('deletebond %s.%s %s.%s' %
                                      (ur, istr, ur, name) ) # not necessary
                    leap_pert0.append('bond %s.%s %s.%s' % (ur, istr, ur, name) )

        if fstr.startsWith('H') and con_morph.nConnections(iinfo.index) > 1:
            ur = 's.1' #% new.residue().number().value()
            leap_pert1.append('set %s.%s pert true' % (ur, fstr) )

            for bond_index in con_morph.connectionsTo(iinfo.index):
                atom1 = lig_morph.select(bond_index)
                #resnum = atom1.residue().number().value()

                rname = util.search_atominfo(atom1.index(), atom_map)
                name = '%s' % rname.name.value()
                ur = 's.1' #% resnum
                
                # FIXME: hard-coded unit name
                leap_pert1.append('deletebond %s.%s %s.%s' % (ur, fstr, ur, name) )
                leap_pert1.append('bond %s.%s %s.%s' % (ur, fstr, ur, name) )


        if not finfo.atom:
            charge = 0.0 * Sire.Units.mod_electron
            ambertype = const.DUMMY_TYPE
        else:
            base = lig_final.atoms().select(finfo.index)

            charge = base.property('charge')
            ambertype = '%s' % base.property('ambertype')

            new.setProperty('element', Sire.Mol.Element(istr) )

        # tag atom name because Sire throws exception when duplicate atom names
        new.rename(Sire.Mol.AtomName('\x00%s\x00' % fstr ) )

        new.setProperty('ambertype', ambertype)
        new.setProperty('charge', charge)

        # NOTE: state1 is already MolEditor but if we do not call molecule() on
        #       it, Sire segfaults...
        state1 = new.molecule()

    # second pass to remove atom name tags
    for atom in state1.molecule().atoms():
        new = state1.atom(atom.index() )

        name = atom.name().value().replace('\x00', '')

        new.rename(Sire.Mol.AtomName(name) )
        state1 = new.molecule()

    with open(const.LEAP_PERT0_FILE, 'w') as lp:
        lp.write('\n'.join(leap_pert0) )
        lp.write('\n')

    with open(const.LEAP_PERT1_FILE, 'w') as lp:
        lp.write('\n'.join(leap_pert1) )
        lp.write('\n')

    state1 = state1.commit()

    return state1
