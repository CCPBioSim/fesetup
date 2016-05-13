#  Copyright (C) 2014  Hannes H Loeffler
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
MD engine base
"""

__revision__ = "$Id$"



import os, sys

from parmed.amber.mask import AmberMask
from parmed.amber.readparm import AmberParm

from FESetup import const



MIN_PREFIX = 'min'
MD_PREFIX = 'md'
RST_EXT = os.extsep + 'rst7'

_rs = ' ntr = 1, restraint_wt = %.2f,\n restraintmask="%s",'

# FIXME: more specific, may not be just protein+water
#        make FESetup aware what it is dealing with and do this dynamically
_restraint_table = {
    'backbone': '(@CA,C,N,O & !:WAT) | @O3\',C3\',C4\',C5\',O5\',OP1,OP1',
    'bb_lig': '(@CA,C,N,O & !:WAT) | @O3\',C3\',C4\',C5\',O5\',OP1,OP1 | '
              '(:LIG & !@H=)',
    'heavy': '!:WAT & !@H=',
    'protein': ':' + ','.join(const.AMBER_PROTEIN_RESIDUES),
    'nucleic': ':' + ','.join(const.AMBER_NUCLEIC_RESIDUES),
    'notligand': '!:' + const.LIGAND_NAME,
    'notsolvent': '!:WAT,HOH,T3P,T4P'
    }


class MDEBase(object):
    """
    MD engine base.
    """

    def __init__(self):
        self.run_no = 1


    def mask_indexes(self, parmtop, mask):
        """
        Create AMBER mask index generator from parmtop file.

        :param parmtop: parmtop filename
        :type parmtop: string
        :param mask: AMBER mask
        :type mask: integer
        :returns: mask index generator
        """
                
        p = AmberParm(parmtop)
        m = AmberMask(p, mask)

        return m.Selected()


    def _write_rst7(self, natoms, xx, yy, zz, coords, vels, center = 'False'):
        """
        Write AMBER .rst7 file

        Format:
        20A4 title
        I5,5E15.7 #atoms time(ps)    why 5E15.7?
        6F12.7 coords
        ...
        6F12.7 vels (if MD)
        ...
        6F12.7 box size

        :returns: file name of created rst7 file
        """

        delx = dely = delz = 0.0

        # FIXME: only cuboid box
        if center:
            minc = [ sys.float_info[0],  sys.float_info[0],  sys.float_info[0] ]
            maxc = [-sys.float_info[0], -sys.float_info[0], -sys.float_info[0] ]

            for i in range(0, len(coords), 3):
                if coords[i]   < minc[0]:
                    minc[0] = coords[i]
                if coords[i+1] < minc[1]:
                    minc[1] = coords[i+1]
                if coords[i+2] < minc[2]:
                    minc[2] = coords[i+2]

                if coords[i]   > maxc[0]:
                    maxc[0] = coords[i]
                if coords[i+1] > maxc[1]:
                    maxc[1] = coords[i+1]
                if coords[i+2] > maxc[2]:
                    maxc[2] = coords[i+2]

            # FIXME: do we have do consider vdW radii?
            delx = minc[0] - (xx - maxc[0] + minc[0]) / 2
            dely = minc[1] - (yy - maxc[1] + minc[1]) / 2
            delz = minc[2] - (zz - maxc[2] + minc[2]) / 2


        with open(self.prev + RST_EXT, 'w') as rst7:
            rst7.write('converted with FESetup\n')
            rst7.write('%5i%15.7f\n' % (natoms, 0.0) )

            cnt = 0
            nl_done = False

            for i in range(0, len(coords), 3):
                cnt += 3
                nl_done = False

                x = coords[i]   - delx
                y = coords[i+1] - dely
                z = coords[i+2] - delz

                rst7.write('%12.7f%12.7f%12.7f' % (x, y, z) )

                if cnt >= 6:
                    rst7.write('\n')
                    cnt = 0
                    nl_done = True

            if cnt < 6 and not nl_done:
                rst7.write('\n')

            cnt = 0
            nl_done = False

            for coord in vels:
                cnt += 1
                nl_done = False

                rst7.write('%12.7f' % coord)

                if cnt >= 6:
                    rst7.write('\n')
                    cnt = 0
                    nl_done = True

            if cnt < 6 and not nl_done:
                rst7.write('\n')

            # FIXME: only cuboid box
            rst7.write('%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n' %
                       (xx, yy, zz, 90.0, 90.0, 90.0) )

        return self.prev + RST_EXT
