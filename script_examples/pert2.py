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

import FESetup as setup
import FESetup.prepare as prep
import FESetup.mutate as mutate


setup.create_logger('pert2.log')

base = 'ligands'

# force field, sub type, water model, divalent ions, MD engine
amber = prep.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')

protein = amber.Protein('data', './thrombin/protein',
                             start_file = 'min.pdb')
ligand1 = amber.Ligand('3A', base)
ligand2 = amber.Ligand('3B', base)

complex = amber.Complex(protein, ligand1)

with mutate.Morph(ligand1, ligand2) as mol:
  mol.map_atoms()
  mol.create_coords()
  mol.make_pert_file()
  mol.create_coords(ligand1)
  mol.create_coords(complex)
