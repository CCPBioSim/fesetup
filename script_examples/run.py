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

import sys
from FESetup import create_logger, prepare, errors


create_logger('test.log')

# prepare protein
base = './thrombin/protein'
protein_name = 'data'

# force field, sub type, water model, divalent ions, add on ff, MD engine
amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')

with amber.Protein(protein_name, base, start_file = '2ZC9sims.pdb') as protein:
  protein.protonate_propka(pH = 7.0)
  protein.create_top(boxtype = '')

  protein.setup_MDEngine('sander')

  protein.minimize(namelist = '%FIXH', restraint = 'heavy',
                          nsteps = 50)
  protein.minimize(namelist = '%ALL', restraint = 'backbone',
                          nsteps = 450)


# prepare ligands
base = './thrombin/poses'
ligand_names = ['3A','3B','3C','3D','3E','3F','3G','3H','3I','3J','3K','5A','5B','5C','5D','5E','5F','5G','5H','5I','5J','5K']

failed = []

for name in ligand_names:
  with amber.Ligand(name, base) as ligand:
    try:
      ligand.prepare()
      ligand.param()

      ligand.create_top(boxtype = '')
      ligand.flex()

      ligand.conf_search(numconf = 100, geomsteps = 10)
      ligand.align()

      ligand.create_top(boxtype = 'rectangular', neutralize = True)

      ligand.setup_MDEngine('sander')

      ligand.minimize(restraint = 'notsolvent')
      ligand.md(nsteps = 200, restraint = 'notsolvent')

# prepare complexes
      with amber.Complex(protein, ligand) as complex:
        complex.create_top(boxtype = 'rectangular', neutralize = True)

        complex.setup_MDEngine('pmemd.MPI', 'mpirun -np 4')

        complex.minimize(namelist = '%FIXH', restraint = 'heavy', nsteps = 50)
        complex.minimize(namelist = '%ALL', restraint = 'backbone',
                         nsteps = 450)

        complex.prot_flex()
        complex.flatten_rings()

        # FIXME: does not work yet
        #complex.align(filt = True)
 
        complex.md('%HEAT', restraint = 'backbone', nsteps = 5000, T = 300.0)
        complex.md('%PRESS', restraint = 'backbone', nsteps = 5000,
                   T = 300.0, p = 1.0)

        for rforce in range(10, -1, -2):
          complex.md(namelist = '%PRESS', restraint = 'backbone',
                            nsteps = 2000, restr_force = rforce)


    except errors.SetupError, why:
      failed.append(name)
      print '%s failed: %s' % (name, why)

if failed:
  print >>sys.stderr, 'The following ligands have failed:'

  for name in failed:
    print >>sys.stderr, '  %s' % name
