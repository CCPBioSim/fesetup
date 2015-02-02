#  Copyright (C) 2013  Hannes H Loeffler
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
Version 1 of the ModelConfig class which is a DataDict with certain required
dictionary keys.
"""

__revision__ = "$Id$"


import os, time
from FESetup import const
from FESetup.datadict import DataDict, DataDictError



class ModelConfig(DataDict):

    # None is used to signal unset values
    mandatory = {
        'version': 1,
        'name': None,
        'forcefield': None,
        'type': None,
        'crd.filename': None,
        'crd.filetype': None,
        'top.filename': None,
        'top.filetype': None,
        'charge.filename': None,
        'charge.filetype': None,
        'charge.type': None,
        'isvalid': None,
        'ismorph': None,
        'supports_md': None,
        'supports_mc': None,
        'data.checksum': None ,
        'data.checksum_type': 'sha1',
        'data.compression_type': 'bz2'
        }

    optional = {
        'header': None,
        'time': None,
        'mc_type': None,
        'state0': None,
        'state1': None
        }


    def __init__(self, name = 'unnamed'):
        super(ModelConfig, self).__init__(self.mandatory)

        #self.update(self.optional)

        self['name'] = name
        self.files = set()
        self.filename = self['name'] + const.MODEL_EXT


    def add_file(self, file):
        self.files.add(file)


    def remove_files(self):
        self.files = set()


    def write(self):
        self['data.checksum'] = self.add_files(self.files,
                                               self['data.checksum_type'],
                                               self['data.compression_type'])
        for elem in self.mandatory:
            if elem not in self:
                raise DataDictError('mandatory field "%s" not present' %
                                    self[elem])
            else:
                if self[elem] == None:
                    raise DataDictError('mandatory field "%s" not set' %
                                        elem)

        self['isvalid'] = True
        self['time'] = time.ctime()

        super(ModelConfig, self).write(self.filename)
