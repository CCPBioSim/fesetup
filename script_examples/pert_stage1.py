#  Copyright (C) 2012-2014,2016  Hannes H Loeffler
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



import os
import sys
from FESetup import create_logger, prepare, errors, DirManager


create_logger('pert.log')

top = os.getcwd()

# force field, sub type, water model, divalent ions, MD engine
amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')

ligand_names = ['3A', '3B']

failed = []

for name in ligand_names:
  ligand_file = os.path.join(top, 'thrombin/poses/%s/ligand.pdb' % name)
  ligand_wd = os.path.join(top, '_ligand', name)

  print 'Making ligand %s...' % name
  ligand = amber.Ligand(name, ligand_file)

  with DirManager(ligand_wd):
    try:
      ligand.prepare()
      ligand.param()

      ligand.create_top(boxtype='')
      ligand.flex()

      ligand.create_top(boxtype='rectangular', neutralize=True)

      ligand.setup_MDEngine('sander')

      ligand.minimize(namelist='%ALL', restraint='notsolvent')
      ligand.md(namelist='%HEAT', nsteps=500, restraint='notsolvent', T=300.0)
      ligand.md(namelist='%CONSTT', nsteps=500, restraint='notsolvent', T=300.0)

    except errors.SetupError, why:
      failed.append(name)
      print '%s failed: %s' % (name, why)

if failed:
  print >>sys.stderr, 'The following ligands have failed:'

  for name in failed:
    print >>sys.stderr, '  %s' % name
