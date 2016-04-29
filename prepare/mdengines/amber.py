#  Copyright (C) 2014,2016  Hannes H Loeffler
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
Amber MD engine: sander or pmemd
"""

__revision__ = "$Id$"



import os

import mdebase
from FESetup import errors, logger
from FESetup.prepare.amber import utils



def is_periodic(topfile):
    """
    Check if AMBER topology file was made for a periodic system.

    :param topfile: the topology file to be checked for periodicity
    :type topfile: string
    :rtype: bool
    """

    with open(topfile, 'r') as top:
        for line in top:
            if line.startswith("%FLAG BOX_DIMENSIONS"):
                return True

    return False


class MDEngine(mdebase.MDEBase):
    """
    Amber MD engine.
    """

    # FIXME: files are specific to AMBER but needed for conversion for
    #        other MD packages
    def __init__(self, amber_top, amber_crd, sander_crd, sander_rst,
                 amber_pdb, box_dims=None, solvent=None, mdprog='sander',
                 mdpref='', mdpost=''):

        super(MDEngine, self).__init__()

        # FIXME: do not assume that top/crd are in AMBER format
        self.update_files(amber_top, amber_crd, sander_crd, sander_rst,
                          amber_pdb)

        self.mdpref = mdpref

        self.mdprog = ''
        self._self_check(mdprog)


    def update_files(self, amber_top, amber_crd, sander_crd, sander_rst,
                     amber_pdb):
        self.amber_top = amber_top
        self.amber_crd = amber_crd
        self.sander_crd = sander_crd
        self.sander_rst = sander_rst
        self.amber_pdb = amber_pdb

        if is_periodic(self.amber_top):
            self.min_periodic = ' ntb = 1,\n'
            self.md_periodic = ''           # must be in pre-defined namelist
        else: # FIXME: this needs a closer look
            self.min_periodic = ' ntb = 0, cut = 99999.0,\n'
            self.md_periodic = ' ntb = 0, igb = 2, cut = 16.0, nrespa = 2\n'


    def minimize(self, namelist='%ALL', nsteps=100, ncyc=10, restr_str='',
                 restr_force=5.0):
        """
        Use the AMBER/sander module to minimize a system.

        :param config: MD paramers, if start with '%' respective default is
           chosen, otherwise a full sander namelist as string
        :type config: string
        :param nsteps: maximum number of conjugated gradient steps
        :type nsteps: integer
        :param ncyc: number of initial steepest decent steps
        :type ncyc: float
        :param restr_str: pre-defined restraints or AMBER mask
        :type restr_str: string
        :param restr_force: force constant for the restraints in kcal/mol/AA
        :type restr_force: float
        :raises: SetupError
        """

        if namelist[0] == '%':
            pname = namelist[1:]
            logger.write('Running minimisation with protocol %s' % pname)

            try:
                namelist = PROTOCOLS['MIN_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such min protocol predefined: '
                                        '%s' % namelist)

            restraint = ''

            if restr_str:
                try:
                    mask = mdebase._restraint_table[restr_str]
                except KeyError:
                    mask = restr_str

                restraint = (mdebase._rs % (restr_force, mask) )

            namelist = namelist.format(nsteps, ncyc, nsteps / 5,
                                       self.min_periodic + restraint)

        self._run_mdprog(mdebase.MIN_PREFIX, namelist, mask, False)


    def md(self, namelist='', nsteps=1000, T=300.0, p=1.0,
           restr_str='', restr_force=5.0, nrel=1, wrap=True, dt=0.002):
        """
        Use the sander module from AMBER to run molecular dynamics on a system.

        :param namelist: MD paramers, if start with '%' respective default is
           chosen, otherwise a full sander namelist as string
        :type namelist: string
        :param nsteps: maximum number of MD steps
        :type nsteps: integer
        :param T: temperature
        :type T: float
        :param p: pressure
        :type p: float
        :param restr_str: pre-defined restraints or AMBER mask
        :type restr_str: string
        :param restr_force: force constant for the restraints in kcal/mol/AA
        :type restr_force: float
        :param nrel: number of restraint relaxation steps
        :type nrel: integer
        :param wrap: wrap coordinates in a periodic system
        :type wrap: bool
        :param dt: timestep in ps
        :type dt: float
        :raises: SetupError
        """

        constp = False

        if namelist[0] == '%':
            pname = namelist[1:]
            logger.write('Running MD with protocol %s' % pname)

            if pname == 'PRESS' or pname == 'SHRINK':
                constp = True

            try:
                namelist = PROTOCOLS['MD_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such MD protocol predefined: '
                                        '%s' % namelist)

            restraint = ''

            if restr_str:
                try:
                    mask = mdebase._restraint_table[restr_str]
                except KeyError:
                    mask = restr_str

                restraint = (mdebase._rs % (restr_force, mask) )

            namelist = namelist.format(nsteps, nsteps / 5, nsteps / 10,
                                       self.md_periodic + restraint,
                                       '{:d}'.format(wrap), T, p, dt)


        self._run_mdprog(mdebase.MD_PREFIX, namelist, mask, constp)


    def get_box_dims(self):
        """
        Extract box information from rst7 file.

        :returns: box dimensions
        """

        with open(self.sander_rst, 'r') as rst:
            for line in rst:
                box_dims = line

        # FIXME: rectangular box only
        return [float(d) for d in box_dims.split()]


    def _run_mdprog(self, prefix, namelist, mask, constp):
        """
        Run sander/pmemd from the Amber package.
        """

        prefix += '%05i'

        prefix = prefix % self.run_no
        self.sander_rst = prefix + mdebase.RST_EXT

        with open(prefix + os.extsep + 'in', 'w') as mdin:
            mdin.writelines(namelist)

        # NOTE: we assume trajectory will be written in NetCDF
        #       format (ioutfm = 1)
        flags = ('-i {0}.in -p {1} -c {2} '
                 '-O -o {0}.out -e {0}.en -x {0}.nc -inf {0}.info -r {3}')

        if mask:
            # Constant pressure with positional restraints shifts coordinates.
            if constp:
                flags += ' -ref %s' % self.sander_crd
            else:
                flags += ' -ref %s' % self.amber_crd

        err = utils.run_amber(self.mdpref + ' ' + self.mdprog,
                              flags.format(prefix, self.amber_top,
                                           self.sander_crd, self.sander_rst) )

        if err:
            logger.write('sander/pmemd failed with message %s' % err[1])
            raise errors.SetupError('error in sander run %s: %s' %
                                    (prefix, err[1]) )

        self.sander_crd = self.sander_rst
        self.run_no += 1


    def to_rst7(self):
        pass


    def _self_check(self, mdprog):
        """
        Check NAMD installation.

        :param mdprog: the namd executable provided to the class
        :type mdprog: str
        """

        if not 'AMBERHOME' in os.environ:
            raise errors.SetupError('AMBERHOME not set')

        self.mdprog = os.path.join(os.environ['AMBERHOME'], 'bin', mdprog)

        if not os.access(self.mdprog, os.X_OK):
            raise errors.SetupError('AMBERHOME does not have a %s binary' %
                                    mdprog)


PROTOCOLS = dict(
    MIN_FIXH = '''Fix hydrogens
 &cntrl
   imin = 1, ntmin = 1,
   maxcyc = {0}, ncyc = {1},
   ntpr = {2}, ntwe = {2},
   dx0 = 1.0D-7,
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   ntc = 2, noshakemask = '!:WAT,HOH,T3P,T4P,T4E',
   {3}
 /
''',

    MIN_ALL = '''Minimise whole system
 &cntrl
   imin = 1, ntmin = 1,
   maxcyc = {0}, ncyc = {1},
   ntpr = {2}, ntwe = {2},
   dx0 = 1.0D-7,
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   {3}
 /
''',

    MD_SHRINK = '''shrink box
 &cntrl
   imin = 0, nstlim = {0}, irest = 0, ntx = 1, dt = {7},
   ntt = 1, temp0 = {5}, tautp = 0.5,
   ntp = 1, pres0 = {6}, taup = 0.1,
   ntb = 2,
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = {4},
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   ntwe = {1}, ntwx = {1}, ntpr = {2},
   {3}
 /
''',
    
    MD_RANDOMV = '''Maxwell velocities
 &cntrl
   imin = 0, nstlim = {0}, irest = 0, ntx = 1, dt = {7},
   ntt = 1, temp0 = {5}, tempi = {5}, tautp = 1.0,
   ntb = 1,
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = {4},
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   ntwe = {1}, ntwx = {1}, ntpr = {2},
   {3}
 /
''',
    
# we may later want to add ntxo = 2, i.e. restart file in NetCDF format
    MD_HEAT = '''heat the system
 &cntrl
   imin = 0, nstlim = {0}, irest = 0, ntx = 1, dt = {7},
   nmropt = 1,
   ntt = 1, temp0 = {5}, tempi = 5.0, tautp = 1.0,
   ntb = 1, pres0 = {6},
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = {4},
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   ntwe = {1}, ntwx = {1}, ntpr = {2},
   {3}
 /

 &wt
   type = 'TEMP0',
   istep1 = 0, istep2 = {0},                                      
   value1 = 5.0, value2 = 300.0
 /

 &wt type = 'END'
 /
''',

    MD_CONSTT = '''constant temperature
 &cntrl
   imin = 0, nstlim = {0}, irest = 1, ntx = 5, dt = {7},
   ntt = 1, temp0 = {5}, tautp = 1.0,
   ntb = 1,
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = {4},
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   ntwe = {1}, ntwx = {1}, ntpr = {2},
   {3}
 /
''',

# FIXME: isotropic pressure scaling only
    MD_PRESS = '''constant pressure
 &cntrl
   imin = 0, nstlim = {0}, irest = 1, ntx = 5, dt = {7},
   ntt = 1, temp0 = {5}, tautp = 1.0,
   ntp = 1, pres0 = {6}, taup = 0.5,
   ntb = 2,
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = {4},
   ntwe = {1}, ntwx = {1}, ntpr = {2},
   ntxo = 1,  ! 2 is default in AMBER16 for NetCDF file
   {3}
 /
'''
)
