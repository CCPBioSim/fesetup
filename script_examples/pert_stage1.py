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
from FESetup.modelconf import ModelConfig

import modelio


create_logger('pert_stage1.log')

top = os.getcwd()

# force field, sub type, water model, divalent ions, MD engine
amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')

# make protein
protein_file = os.path.join(top, 'thrombin/protein/2ZC9/protein.pdb')
protein_wd = os.path.join(top, '_protein', '2ZC9')

name = '2ZC9'
print 'Making protein %s...' % name
protein = amber.Protein(name, protein_file)
model = ModelConfig(name)
model_filename = name + os.extsep + 'model'

with DirManager(protein_wd):
    protein.protonate_propka(pH=7.0)
    protein.get_charge()                  # must be done explicitly!
    protein.prepare_top()                 # must be done explicitly!
    protein.create_top(boxtype = '')

    # will look in AMBERHOME for this
    protein.setup_MDEngine('sander.MPI', 'mpirun -np 2')
    protein.minimize(namelist='%ALL', restraint='backbone', nsteps=100)

    model['charge.total'] = protein.charge
    model['forcefield'] = 'AMBER/ff14SB'
    model['molecule.type'] = 'biomolecule'

    model['crd.original'] = protein.mol_file
    model.add_file(protein.mol_file)

    modelio.save_model(model, protein, model_filename, '..')


# make ligand
ligand_names = ['3A', '3B']

failed = []

for name in ligand_names:
    ligand_file = os.path.join(top, 'thrombin/poses/%s/ligand.pdb' % name)
    ligand_wd = os.path.join(top, '_ligands', name)

    print 'Making ligand %s...' % name
    ligand = amber.Ligand(name, ligand_file)
    model = ModelConfig(name)
    model_filename = name + os.extsep + 'model'

    with DirManager(ligand_wd):
        try:
            ligand.prepare()
            ligand.param()

            ligand.prepare_top()
            ligand.create_top(boxtype='')
            ligand.create_top(boxtype='rectangular', neutralize=True)

            ligand.setup_MDEngine('sander')

            ligand.minimize(namelist='%ALL', restraint='notsolvent')
            ligand.md(namelist='%HEAT', nsteps=500, restraint='notsolvent',
                      T=300.0)
            ligand.md(namelist='%CONSTT', nsteps=500, restraint='notsolvent',
                      T=300.0)

            model['charge.total'] = ligand.charge
            model['charge.method'] = 'AM1-BCC'
            model['forcefield'] = ligand.gaff
            model['molecule.type'] = 'ligand'

            model['crd.original'] = ligand.mol_file
            model.add_file(ligand.mol_file)

            modelio.save_model(model, ligand, model_filename, '..')

        except errors.SetupError, why:
            failed.append(name)
            print '%s failed: %s' % (name, why)

if failed:
    print >>sys.stderr, 'The following ligands have failed:'

    for name in failed:
        print >>sys.stderr, '  %s' % name
