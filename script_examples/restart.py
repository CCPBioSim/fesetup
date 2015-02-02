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

#
# simple example to demonstrate how to start from previously and successfully
# constructed ligand and receptor
#

import sys
import FESetup as setup
import FESetup.prepare as prep

# force field, sub type, water model, divalent ions, MD engine
amber = prep.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')

setup.create_logger('test.log')

protein = amber.Protein('data', './thrombin/protein',
                             start_file = 'min.pdb')
ligand = amber.Ligand('3K', './thrombin/poses')

with amber.Complex(protein, ligand) as complex:
  complex.create_top(boxtype = 'rectangular', neutralize = True)
  complex.setup_MDEngine('sander')
  complex.minimize(nsteps = 200, restraint = 'backbone')
  complex.prot_flex()
  complex.md(nsteps = 200, restraint = 'notsolvent')
  complex.flatten_rings()
