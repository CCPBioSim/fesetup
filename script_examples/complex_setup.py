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


# This script is a simple demonstration on how to set up several complexces
# from a single protein with several ligands



import os
import sys

from FESetup import create_logger, prepare, errors, DirManager


create_logger('complex_setup.log')

top = os.getcwd()

# force field, sub type, water model, divalent ions, add on ff, MD engine
amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')

# make protein
protein_file = os.path.join(top, 'thrombin/protein/2ZC9/protein.pdb')
protein_wd = os.path.join(top, '_protein', '2ZC9')

print 'Making protein 2ZC9...'
protein = amber.Protein('2ZC9', protein_file)

with DirManager(protein_wd):
    protein.protonate_propka(pH=7.0)
    protein.get_charge()                  # must be done explicitly!
    protein.prepare_top()                 # must be done explicitly!
    protein.create_top(boxtype = '')

    # will look in AMBERHOME for this
    protein.setup_MDEngine('sander.MPI', 'mpirun -np 2')
    protein.minimize(namelist='%FIXH', restraint='heavy', nsteps=50)
    protein.minimize(namelist='%ALL', restraint='backbone', nsteps=150)

# make ligands
ligand_names = ['3A','3B','3C','3D','3E','3F','3G','3H','3I','3J','3K',
  '5A','5B','5C','5D','5E','5F','5G','5H','5I','5J','5K']

failed = []

for name in ligand_names:
    ligand_file = os.path.join(top, 'thrombin/poses/%s/ligand.pdb' % name)
    ligand_wd = os.path.join(top, '_ligands', name)

    print 'Making ligand %s...' % name
    ligand = amber.Ligand(name, ligand_file)

    with DirManager(ligand_wd):
        try:
            ligand.prepare()
            ligand.param()

            ligand.prepare_top()
            ligand.create_top(boxtype = '')

            #ligand.conf_search(numconf=100, geomsteps=10)
            #ligand.align()

            ligand.create_top(boxtype='rectangular', neutralize=True)

            ligand.setup_MDEngine('sander')
            ligand.minimize(namelist='%ALL', restraint='notsolvent')
            ligand.md(namelist='%HEAT', nsteps=200, restraint='notsolvent',
                      T=300.0)
            ligand.md(namelist='%CONSTT', nsteps=200, restraint='notsolvent',
                      T=300.0)

            print 'Making complex with %s...' % name
            complex = amber.Complex(protein, ligand)

            # make complexes
            with DirManager(os.path.join(top, '_complex', '2ZC9-%s' % name)):
                complex.copy_files((ligand_wd, protein_wd),
                                   (ligand.orig_file, ligand.frcmod,
                                    protein.orig_file, 'ligand.ac'))

                complex.ligand_fmt = ligand.mol_fmt
                complex.prepare_top()
                complex.create_top(boxtype = 'rectangular', neutralize = True)

                complex.setup_MDEngine('pmemd.MPI', 'mpirun -np 4')

                complex.minimize(namelist='%FIXH', restraint='heavy', nsteps=50)
                complex.minimize(namelist='%ALL', restraint='backbone', nsteps=150)

                #complex.prot_flex()
                #complex.flatten_rings()

                # these MD steps will take a long time...
                complex.md('%HEAT', restraint='backbone', nsteps=2000, T=300.0)
                complex.md('%PRESS', restraint='backbone', nsteps=5000, T=300.0,
                           p=1.0)

                for rforce in range(4, -1, -2):
                    complex.md(namelist='%PRESS', restraint='backbone', nsteps=500,
                               restr_force=rforce)


        except errors.SetupError, why:
            failed.append(name)
            print '%s failed: %s' % (name, why)

if failed:
    print >>sys.stderr, 'The following ligands have failed:'

    for name in failed:
        print >>sys.stderr, '  %s' % name
