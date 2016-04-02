#  Copyright (C) 2016  Hannes H Loeffler, Julien Michel
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
Leap related functionality
"""

__revision__ = "$Id$"



class Leap(object):
    """
    Hold data/commands for leap input.
    """

    def __init__(self, ff, sload):
        """
        :param ff: force fields
        :type ff: list
        :param sload: solvent leap commands
        :type sload: string
        """

        self.force_fields = set(ff)
        self.solvents = sload
        self.user_params = []
        self.mols = []


    def add_mol(self, mol_file, ftype, mods=[], pert=''):
        """
        Add molecule info.
        
        :param mol_file: filename of the input structure
        :type mol_file: string
        :param ftype: file type of the input structure, either PDB or mol2
        :type ftype: str
        :param mods: file name of frcmods
        :type mods: list
        :param pert: instructions for set pert true
        :type pert: str
        """

        self.mols.append( (mol_file, ftype, mods, pert) )


    def add_force_field(self, ff):
        self.force_fields.add(ff)


    def generate_init(self):
        leap_cmds = []

        for ff in self.force_fields:
            leap_cmds.append('source "leaprc.%s"' % ff)

        leap_cmds.append(self.solvents)

        for up in self.user_params:
            leap_cmds(up)

        load_cmd = {'pdb': 'loadPDB', 'mol2': 'loadmol2'}

        mnames = []
        mol_cnt = 0
        mod_cnt = 0

        for mol in self.mols:
            mol_file, ftype, mods, pert = mol

            if mods:
                for mod in mods:
                    leap_cmds.append('mod%i = loadAmberParams "%s"'
                                     % (mod_cnt, mod) )
                    mod_cnt += 1

            mname = 'cmp' + str(mol_cnt)
            mnames.append(mname)
            leap_cmds.append('%s = %s "%s"' % (mname, load_cmd[ftype], mol_file))

            if pert:
                set_pert = ''
                # FIXME: hard-coded residue number 1
                pert_cmds = ('set {0}.1.{1} pert true\n'
                             'deletebond {0}.1.{1} {0}.1.{2}\n'
                             'bond {0}.1.{1} {0}.1.{2}\n')

                for a, b in pert:
                    if a.startswith('H'):
                        set_pert += pert_cmds.format(mname, a, b)
                    elif b.startswith('H'):
                        set_pert += pert_cmds.format(mname, b, a)
                    else:
                        set_pert += '# unknown: %s, %s\n' % (a, b)

                leap_cmds.append(set_pert)

            mol_cnt += 1

        leap_cmds.append('s = combine {%s}\n' % (' '.join(mnames)))

        return '\n'.join(leap_cmds)
