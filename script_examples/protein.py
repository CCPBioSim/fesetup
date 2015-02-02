#  Copyright (C) 2014  Hannes H Loeffler
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
#
#  Setup of a protein with FESetup
#

import sys
from FESetup import create_logger, prepare, errors


create_logger('glnbp.log')

# prepare protein
base = './bench'
protein_name = 'GlnBP'

# force field, sub type, water model, divalent ions, add on ff, MD engine
amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')
#amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'namd')
#amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'gromacs')

with amber.Protein(protein_name, base) as protein:
  protein.protonate_propka(pH = 7.0)
  protein.create_top(boxtype = 'rectangular', boxlength = 10.0)

  # binary, prefix, postfix
  protein.setup_MDEngine('pmemd.MPI', 'mpirun -np 4')
  #protein.setup_MDEngine('namd2', '', '+p4 +isomalloc_sync')
  #protein.setup_MDEngine('mdrun', '', '-nt 4')

  protein.minimize('%ALL', restraint = 'heavy', nsteps = 300, ncyc = 100)
  #protein.minimize('%STD', restraint = 'heavy', nsteps = 300, ncyc = 100)

  protein.md('%RANDOMV', restraint = 'backbone', restr_force = 5.0,
             nsteps = 2000, T = 300.0)
  protein.md('%PRESS', restraint = 'backbone', restr_force = 2.0,
             nsteps = 5000, T = 300.0, p = 1.0)
