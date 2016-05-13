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
import warnings

from FESetup import create_logger, prepare, errors
from FESetup.modelconf import ModelConfig

# FIXME: That's here solely to suppress a warning over a fmcs/Sire double
# data type registration collision.  Impact limited as much as possible but
# still potentially dangerous.  Fix actual problem instead!
with warnings.catch_warnings():
    warnings.filterwarnings('ignore',
                            'to-Python converter for.*already registered')

    from FESetup import mutate

import modelio


create_logger('pert_stage2.log')

# force field, sub type, water model, divalent ions, MD engine
amber = prepare.ForceField('amber', 'ff14SB', 'tip3p', 'cm', [], 'amber')
model = modelio.read_model('_protein/2ZC9.model')
protein = amber.Protein(model['name'])

protein.charge = float(model['charge.total'])
protein.amber_top = model['top.filename']
protein.amber_crd = model['crd.filename']
protein.orig_file = model['crd.original']

if 'box.dimensions' in model:
    protein.box_dims = model['box.dimensions']


ligand_pair = []
workdirs = []

for name in ('3A', '3B'):
    model = modelio.read_model('_ligands/%s.model' % name)
    ligand =  amber.Ligand(name)
  
    ligand.charge = float(model['charge.total'])
    ligand.amber_top = model['top.filename']
    ligand.amber_crd = model['crd.filename']
    ligand.orig_file = model['crd.original']

    if 'box.dimensions' in model:
        ligand.box_dims = model['box.dimensions']

    ligand_pair.append(ligand)
    workdirs.append(os.path.join(os.getcwd(), '_ligands', name))

with mutate.Morph(ligand_pair[0], ligand_pair[1], workdirs[0], workdirs[1], amber,
                  FE_type='pmemd/softcore', separate=True) as mol:
    mol.setup('', '', '')
    #mol.create_coords(ligand1, 'morph', workdirs[0])
