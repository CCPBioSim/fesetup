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
DataDict class.
"""

__revision__ = "$Id$"


import os
import time
import tarfile
import cStringIO
import hashlib



class DataDictError(Exception):
    pass


class DataDict(dict):
    """Simple class for a standard dictionary plus an added data part.

    Arbitrary key-value pairs can be created just like with the built-in
    dict.  A data part consisting of an arbitrary set of file is a mandatory
    part of the class.  This set of file is internally handled via a tar
    (pax) archive and compressed (bz2 or gz).  A checksum is calculated for
    the data.

    The class can be written to a file which is a mixture of ASCII key-value
    pairs and binary data.  The format is each key=value pair on a line by
    itself followed by a END_TAG.  The data part is then written to the final
    line.
    """


    END_TAG = '_END'


    # NOTE: check if subclassing from dict is really a good idea
    def __init__(self, *args, **kwargs):
        super(DataDict, self).__init__(*args, **kwargs)

        self.data = None


    # maybe also overwrite __str__?


    def write(self, filename):
        """Write dictionary header (key/value pairs) and data to filename."""

        if not self.data:
            raise DataDictError('no data files')

        with open(filename, 'wb') as outfile:
            for key, val in sorted(self.iteritems() ):
                outfile.write('%s = %s\n' % (key, val) )

            outfile.write('%s\n' % self.__class__.END_TAG)

            outfile.write(self.data)


    def read(self, filename):
        """Read a datadict file."""

        lineno = 0

        with open(filename, 'rb') as infile:
            while True:
                line = infile.readline()

                if line[:-1] == self.__class__.END_TAG or not line:
                    break

                lineno += 1
                line = line.strip()

                if line == '' or line[0] == '#':
                    continue

                equal_found = False
                key = []
                val = []

                for pos, char in enumerate(line):
                    if char == '=':
                        equal_found = True
                        continue

                    # comments will be lost
                    if char == ' ':
                        if line[pos+1] == '#':
                            break

                    if equal_found:
                        val.append(char)
                    else:
                        key.append(char)

                if not equal_found:
                    raise DataDictError('input line %i does not contain a '
                                        'key/value pair' % lineno)

                key = ''.join(key).strip()
                val = ''.join(val).strip()

                self[key] = val

            data = infile.read()

        if not data:
            raise DataDictError('%s does not contain data files' % filename)

        self.data = data


    # NOTE: merge with write()?
    def add_files(self, files, hash_type = 'sha1',
                       compression_type = 'bz2'):
        """Add a list of files to an interal tar(pax) archive.  A manifest is
        automatically created.  The tar(pax) archive will be compressed.  A
        checksum will be computed for the final archive."""
        
        memtar = cStringIO.StringIO()
        manifest = []

        with tarfile.open(mode = 'w:%s' % compression_type,
                          format = tarfile.PAX_FORMAT,
                          fileobj = memtar) as tar:

            for name in files:
                tar.add(name)

                hash_val = hashlib.new(hash_type)

                with open(name, 'rb') as member:
                    hash_val.update(member.read())

                manifest.append('%s  %s\n' % (hash_val.hexdigest(), name) )

            manifest = ''.join(manifest)

            tinfo = tarfile.TarInfo('manifest')
            tinfo.size = len(manifest)
            tinfo.mtime = time.time()
            tar.addfile(tinfo, cStringIO.StringIO(manifest) )

        data = memtar.getvalue()
        hash_val = hashlib.new(hash_type)
        hash_val.update(data)

        self.data = data

        return hash_val.hexdigest()


    def check_data(self, hash, hash_type):
        """Check checksum of data."""

        if not self.data:
            raise DataDictError('no data files')

        hash_val = hashlib.new(hash_type)
        hash_val.update(self.data)

        if hash != hash_val.hexdigest():
            raise DataDictError('data corruption: checksum doesn\'t match')


    def list(self, compression_type = '*'):
        """List contents of internal tar(pax) archive.

        Compression type by default is transparent.  Other possible types
        are 'gz' and 'bz2'.
        """

        if not self.data:
            raise DataDictError('not data files')

        with tarfile.open(mode = 'r:%s' % compression_type,
                          format = tarfile.PAX_FORMAT,
                          fileobj = cStringIO.StringIO(self.data) ) as tar:
            tar.list()


    def extract(self, compression_type = '*', direc = '.'):
        """Extract contents of internal tar(pax) archive.

        Compression type by default is transparent.  Other possible types
        are 'gz' and 'bz2'.
        """

        if not self.data:
            raise DataDictError('not data files')

        with tarfile.open(mode = 'r:%s' % compression_type,
                          format = tarfile.PAX_FORMAT,
                          fileobj = cStringIO.StringIO(self.data) ) as tar:
            tar.extractall(direc)


if __name__ == '__main__':
    import sys
    import time

    defaults = {
        'version' : 1,
        'nonsense' : 'still nonsense',
        'name' : 'Super Drug Molecule(TM)',
        'type' : 'ligand',
        'forcefield' : 'gaff',
        'supports_md' : False,
        'supports_mc' : False,
        'is_valid' : False
        }

    opts = DataDict(defaults)
    opts['time'] = time.ctime()

    # ModelConfig does not have any default entries!
    opts['data.checksum_type'] = 'sha1'
    opts['data.compression_type'] = 'bz2'
    opts['data.checksum'] = \
        opts.add_files(sys.argv[1:],
        opts['data.checksum_type'], opts['data.compression_type'])

    del opts['nonsense']

    opts['is_valid'] = True
    opts.write('test.model')
    opts.list()

    new_opts = DataDict()

    new_opts.read('test.model')

    new_opts.check_data(new_opts['data.checksum'], new_opts['data.checksum_type'])
    new_opts.list()
    #new_opts.extract(direc = '/tmp')

    new_opts.write('test2.model')
