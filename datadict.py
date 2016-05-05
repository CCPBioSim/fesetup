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
The DataDict class to store and retrieve a dictionary plus associated
data files.

Note: this holds the data package in memory so may need to be revised later.
"""

__revision__ = "$Id$"


import os
import time
import tarfile
import cStringIO
import hashlib



def strip_eol_comment(s, comment_chars='#;', space_chars=' \t'):
    """
    Strip end-of-line comments which are a comment character preceded by a
    space character.

    :param comment_chars: all recognized comment characters
    :type comment_chars: str
    :param space_chars: all recognized white-space character
    :type space_chars: str
    """

    if len(s) < 3:
        return s

    for pos in range(1, len(s)):
        if s[pos] in comment_chars and s[pos-1] in space_chars:
            return s[:pos-1]

    return s


class DataDictError(Exception):
    pass


class DataDict(dict):
    """
    Simple class for a standard dictionary plus an added data package

    Arbitrary key-value pairs can be created just like with the built-in
    dict.  The data part consisting of an arbitrary set of files is a mandatory
    part of the class.  This set of file is internally handled via a tar
    (pax) archive and compressed (bz2 or gz).  A hash is calculated for
    the data and a manifest added to the archive.  The manifest contains the
    individual hashes of the archive data.  The hashes are also stored for
    each file in the comment field of the pax header.

    The class can be written to a file as a mixture of ASCII key-value pairs
    and the binary data archive.  The format is each key=value pair on a line
    by itself followed by the END_TAG.  The data part is then written to the
    final (logical) line.
    """


    _END_TAG = '_END'


    # NOTE: check if subclassing from dict is really a good idea
    def __init__(self, *args, **kwargs):
        super(DataDict, self).__init__(*args, **kwargs)

        self.data = None


    def write(self, filename):
        """
        Write dictionary header (key/value pairs) and data to filename.

        :param filename: the file name to be write to
        :type filename: string
        """

        if not self.data:
            raise DataDictError('no data files')

        with open(filename, 'wb') as outfile:
            for key, val in sorted(self.iteritems() ):
                outfile.write('%s = %s\n' % (key, val) )

            outfile.write('%s\n' % self.__class__._END_TAG)
            outfile.write(self.data)


    def read(self, filename):
        """
        Read a datadict file.  Strips all comments.

        :param filename: the file name to be read from
        :type filename: string
        """

        lineno = 0

        with open(filename, 'rb') as infile:
            while True:
                line = infile.readline()

                if line[:-1] == self.__class__._END_TAG or not line:
                    break

                lineno += 1
                line = line.strip()

                if not line or line[0] == '#':
                    continue

                try:
                    key, val = line.split('=', 1)
                except ValueError:
                    raise DataDictError('input line %i does not contain a '
                                        'key = value pair' % lineno)

                key = key.strip()
                val = strip_eol_comment(val).strip()

                self[key] = val

            data = infile.read()        # FIXME: maybe we just keep a pointer
                                        #        to the archive location?

        if not data:
            raise DataDictError('%s does not contain data files' % filename)

        self.data = data


    def add_files(self, files, hash_type='sha1', compression_type='bz2'):
        """
        Add a list of files to an internal tar(pax) archive.  A manifest is
        automatically created and contains the hashes of every file.  Hashes
        are also stored in the comment field of the PAX header of every file.
        The tar(pax) archive will be compressed.  A hash will be computed for
        the final archive.

        :param files: the file names to be added
        :type files: set of strings
        :param hash_type: the type of hash (see hashlib.algorithms)
        :type hash_type: string
        :param compression_type: the compression type (gz or bz2)
        :type compression_type: string
        """
        
        memtar = cStringIO.StringIO()
        manifest = []

        with tarfile.open(mode = 'w:%s' % compression_type,
                          format = tarfile.PAX_FORMAT,
                          fileobj = memtar) as tar:

            for name in files:
                tinfo = tar.gettarinfo(name)
                hash_val = hashlib.new(hash_type)

                with open(name, 'rb') as member:
                    hash_val.update(member.read())

                hexdig = hash_val.hexdigest()

                tinfo.pax_headers = {'comment': hexdig}
                tar.addfile(tinfo, open(name, 'rb'))

                manifest.append('%s  %s\n' % (hexdig, name) )

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


    def check_data(self, hashv, hash_type):
        """
        Check hash of data against the hash value.

        :param hashv: the hash to be checked against
        :type hashv: string
        :param hash_type: the type of hash as supported by hashlib
        :type hash_type: string
        """

        if not self.data:
            raise DataDictError('no data files')

        hash_val = hashlib.new(hash_type)
        hash_val.update(self.data)

        if hashv != hash_val.hexdigest():
            raise DataDictError('data corruption: hash doesn\'t match')


    def list(self, compression_type='*'):
        """
        List contents of internal tar(pax) archive.

        :param compression_type: the compression type (*=transparent)
        :type compression_type: string
        """

        if not self.data:
            raise DataDictError('no data files')

        with tarfile.open(mode='r:%s' % compression_type,
                          format=tarfile.PAX_FORMAT,
                          fileobj=cStringIO.StringIO(self.data) ) as tar:
            tar.list()


    def extract(self, compression_type='*', direc='.'):
        """
        Extract contents of internal tar(pax) archive.

        :param compression_type: the compression type (*=transparent)
        :type compression_type: string
        :param direc: path to extract to
        :type direc: string
        """

        if not self.data:
            raise DataDictError('no data files')

        with tarfile.open(mode='r:%s' % compression_type,
                          format=tarfile.PAX_FORMAT,
                          fileobj=cStringIO.StringIO(self.data) ) as tar:
            tar.extractall(direc)


    def __str__(self):
        return '\n'.join('{} = {}'.format(key, value)
                         for key, value in self.iteritems()) \
                         + '\n(data package size: {})'.format(len(self.data))



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
    opts['data.hash_type'] = 'sha1'
    opts['data.compression_type'] = 'bz2'
    opts['data.hash'] = \
        opts.add_files(sys.argv[1:],
        opts['data.hash_type'], opts['data.compression_type'])

    del opts['nonsense']

    opts['is_valid'] = True
    opts.write('test.model')
    opts.list()

    new_opts = DataDict()

    new_opts.read('test.model')

    new_opts.check_data(new_opts['data.hash'], new_opts['data.hash_type'])
    new_opts.list()
    #new_opts.extract(direc = '/tmp')

    new_opts.write('test2.model')
    print new_opts
