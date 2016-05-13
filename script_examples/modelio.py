#  Copyright (C) 2016  Hannes H Loeffler
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
import shutil

from FESetup import errors
from FESetup.modelconf import ModelConfig



def read_model(filename):

    model = ModelConfig()
    model.read(filename)

    try:
        model.check_data(model['data.hash'], model['data.hash_type'])

        if not model['is.valid']:
            # FIXME: change exception type
            raise errors.SetupError('invalid model %s' % model['name'])
    except KeyError:
        raise errors.SetupError('invalid model %s' % model['name'])

    model.check_keys()

    return model


def save_model(model, mol, filename, dest_dir):

    # add only latest info
    model['crd.filename'] = mol.amber_crd
    model['crd.filetype'] = 'amber-rst7'
    model['top.filename'] = mol.amber_top
    model['top.filetype'] = 'amber-parm7'

    model.add_file(mol.amber_crd)
    model.add_file(mol.amber_top)

    # FIXME: check for minimisation/equilibration
    if mol.box_dims:
        box_dims = []

        for v in mol.box_dims:
            box_dims.append(float(v))

        model['box.dimensions'] = box_dims
        model['box.format'] = 'bla'  # FIXME: boxlengths-angle

    # the final blessing
    model['is.valid'] = 1

    model.write(filename)

    mname = os.path.join(dest_dir, filename)

    if os.access(mname, os.F_OK):
        os.remove(mname)

    shutil.move(filename, dest_dir)
