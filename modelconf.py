#  Copyright (C) 2013,2016  Hannes H Loeffler
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
Version 0 of the ModelConfig class which is a DataDict with certain required
dictionary keys describing associated files in a packaged archive.
"""

__revision__ = "$Id$"


import time

from datadict import DataDict, DataDictError



class ModelConfig(DataDict):
    """
    A class to hold a dictionary plus associated files.  This is a particular
    implementation of the more generic DataDict.
    """

    # None is used to signal unset values
    mandatory_fields = {
        'version': 0,
        'name': None,
        'charge.total': None,
        'forcefield': None,
        'molecule.type': None,
        'crd.filename': None,
        'crd.filetype': None,
        'top.filename': None,
        'top.filetype': None,
        'data.hash': None ,
        'data.hash_type': 'sha1',
        'data.compression_type': 'bz2',
        'is.valid': None,
        'timestamp': None
        }


    def __init__(self, name=''):
        """
        :param name: the model name
        :type name: string
        """

        super(ModelConfig, self).__init__(self.__class__.mandatory_fields)

        self['name'] = name
        self.files = set()


    def add_file(self, filename):
        """
        Add a file name to this class.

        :param filename: the file name to be added
        :type filename: string
        """

        self.files.add(filename)


    def remove_all_files(self):
        """
        Remove all files from the class.
        """

        self.files = set()


    def write(self, filename):
        """
        Write out this class into a dictionary plus compressed file archive.

        :param filename: the file name
        :type filename: string
        """

        if self.files:
            self['data.hash'] = self.add_files(self.files,
                                               self['data.hash_type'],
                                               self['data.compression_type'])
        self['timestamp'] = time.ctime()

        self.check_keys()

        super(ModelConfig, self).write(filename)


    def check_keys(self):
        """
        Check if the mandatory keys are present.
        """

        for elem in self.__class__.mandatory_fields:
            if elem not in self:
                raise DataDictError('mandatory field "%s" not present' %
                                    self[elem])
            else:
                if self[elem] == None:
                    raise DataDictError('mandatory field "%s" not set' %
                                        elem)
