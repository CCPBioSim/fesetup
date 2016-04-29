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
DL_POLY MD engine
"""

__revision__ = "$Id$"



import os, sys, shutil, math

import mdebase
from FESetup import const, errors, logger
from FESetup.prepare.amber import dlpoly, utils


FIELD_FILENAME = 'FIELD'
CONFIG_FILENAME = 'CONFIG'
CONTROL_FILENAME = 'CONTROL'
STATIS_FILENAME = 'STATIS'
OUTPUT_FILENAME = 'OUTPUT'
REVCON_FILENAME = 'REVCON'
HISTORY_FILENAME = 'HISTORY'

MOVE_LIST = (FIELD_FILENAME, CONFIG_FILENAME, STATIS_FILENAME,
             CONTROL_FILENAME, OUTPUT_FILENAME, HISTORY_FILENAME)

NO_STEPS_PER_100DEG = 2
START_T = 5.0


class MDEngine(mdebase.MDEBase):
    """
    Gromacs MD engine.
    """

    def __init__(self, amber_top, amber_crd, sander_crd, sander_rst,
                 amber_pdb, box_dims = None, solvent = None,
                 mdprog = 'DLPOLY.Z', mdpref = '', mdpost = ''):

        super(MDEngine, self).__init__()

        self.update_files(amber_top, amber_crd, sander_crd, sander_rst,
                          amber_pdb)

        self.mdpref = mdpref
        self.mdpost = mdpost

        self.prev = ''

        self.mdprog = ''
        self._self_check(mdprog)


    def update_files(self, amber_top, amber_crd, sander_crd, sander_rst,
                     amber_pdb):
        self.amber_top = amber_top
        self.amber_crd = amber_crd
        self.sander_crd = sander_crd
        self.sander_rst = sander_rst
        self.amber_pdb = amber_pdb

        self.dlpoly = dlpoly.DLPolyField()
        self.dlpoly.readParm(amber_top, amber_crd)
        self.dlpoly.writeConfig(CONFIG_FILENAME)


    def minimize(self, config = '%STD', nsteps = 100, ncyc = 100,
                 mask = '', restr_force = 5.0):
        """
        Minimize a system.

        :param config: MD paramers, if start with '%' respective default is
           chosen, otherwise a full mdrun input as string
        :type config: string
        :param nsteps: maximum number of conjugated gradient steps
        :type nsteps: integer
        :param ncyc: number of initial steepest decent steps
        :type ncyc: float
        :param mask: pre-defined restraints or AMBER mask
        :type mask: string
        :param restr_force: force constant for the restraints in kcal/mol/AA
        :type restr_force: float
        :raises: SetupError
        """

        suffix = '%05i' % self.run_no

        if config[0] == '%':
            pname = config[1:]
            logger.write('Running minimisation with protocol %s' % pname)

            try:
                config = PROTOCOLS['MIN_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such min protocol predefined: '
                                        '%s' % config)

        self._run_mdprog(suffix, config, mask, restr_force)


    def md(self, config = '', nsteps = 1000, T = 300.0, p = 1.0,
           mask = '', restr_force = 5.0, nrel = 1, wrap = True, dt = 0.002):
        """
        Run molecular dynamics on a system.

        :param config: MD paramers, if start with '%' respective default is
           chosen, otherwise a full mdrun input as string
        :type config: string
        :param nsteps: maximum number of MD steps
        :type nsteps: integer
        :param T: temperature
        :type T: float
        :param p: pressure
        :type p: float
        :param mask: pre-defined restraints or AMBER mask
        :type mask: string
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

        suffix = '%05i' % self.run_no
        cfg = config

        if config[0] == '%':
            pname = config[1:]
            logger.write('Running MD (%s) with protocol %s' %
                         (suffix, pname) )

            try:
                config = PROTOCOLS['MD_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such MD protocol predefined: '
                                        '%s' % config)

        p /= const.ATM2BAR * 1000       # convert to kbar

        # DL_POLY does not have a built-in heating protocol
        if cfg == '%HEAT':
            ninc = int(NO_STEPS_PER_100DEG * T / 100.0)
            tinc = T / (ninc - 1)
            nst = nsteps / ninc

            for i in range(0, ninc):
                if i == 0:
                    T_curr = START_T
                else:
                    T_curr = i * tinc

                ctrl = config.format(T_curr, nst, nst / 5, nst / 10,
                                     nst / 5, p, dt)

                self._run_mdprog(suffix, ctrl, mask, restr_force)

                suffix = '%05i' % self.run_no
        else:
            ctrl = config.format(T, nsteps, nsteps / 5, nsteps / 10,
                                   nsteps / 5, p, dt)

            self._run_mdprog(suffix, ctrl, mask, restr_force)


    def get_box_dims(self):
        """
        Extract box information from rst7 file.

        :returns: box dimensions
        """

        config_file = CONFIG_FILENAME
        line_no = 0
        box_dims = []

        with open(config_file, 'r') as cfg:
            for line in cfg:
                line_no += 1

                if line_no == 1:
                    continue

                if line_no == 2:
                    try:
                        imcon = int(line.split()[1])
                    except ValueError:
                        raise errors.SetupError('invalid line %i in file %s'
                                                % (line_no, config_file) )
                    continue

                if line_no < 5 and imcon > 0:
                    # FIXME: rectangular box only
                    try:
                        x, y, z = line.split()
                        box_dims.append(float(x))

                        x, y, z = cfg.next().split()
                        box_dims.append(float(y))

                        x, y, z = cfg.next().split()
                        box_dims.append(float(z))
                    except ValueError:
                        raise errors.SetupError('invalid cell in file %s'
                                                % config_file)
                    break

        box_dims.extend( (90.0, 90.0, 90.0) )
        return box_dims


    def _run_mdprog(self, suffix, config, mask, restr_force):
        """
        Run DL_POLY executable.
        """

        if mask:
             self._make_restraints(mask, restr_force)

        self.dlpoly.writeField(FIELD_FILENAME)

        with open(CONTROL_FILENAME, 'w') as mdin:
            mdin.writelines(config)

        retc, out, err = utils.run_exe(' '.join((self.mdpref, self.mdprog,
                                                 self.mdpost)))

        if retc:
            logger.write(err)
            raise errors.SetupError('%s has failed (see logfile)' %
                                    self.mdprog)

        for mfile in MOVE_LIST:
            try:
                shutil.move(mfile, mfile + os.extsep + suffix)
            except IOError:             # some files may not be created
                continue

        try:
            shutil.copy2(REVCON_FILENAME, REVCON_FILENAME + os.extsep + suffix)
            shutil.move(REVCON_FILENAME, CONFIG_FILENAME)
        except IOError as why:
            raise errors.SetupError(why)

        self.run_no += 1


    def _make_restraints(self, restr, k):

        self.dlpoly.posres = []

        if k < sys.float_info.epsilon:
            return

        try:
            mask = mdebase._restraint_table[restr]
        except KeyError:
            mask = restr

        mask_idx = list(self.mask_indexes(self.amber_top, mask) )

        for idx in mask_idx:
            self.dlpoly.posres.append( (idx + 1, k) )


    def to_rst7(self):
        """
        Extract coordinates, velocities and box dimensions from REVCON file.
        Units of velocities = Angstrom/ps

        :raises: SetupError
        """

        # CONFIG is actually the most recent REVCON, see _run_mdprog()
        config_file = CONFIG_FILENAME
        self.prev = config_file + os.extsep + '%05i' % (self.run_no - 1)

        line_no = 0
        natoms = 0
        cell = []
        coords = []
        vels = []

        with open(config_file, 'r') as crdvel:
            for line in crdvel:
                line_no += 1

                if line_no == 1: continue

                if line_no == 2:
                    try:
                        levcfg, imcon = line.split()[0:2]
                        lvgcfg = int(levcfg)
                        imcon = int(imcon)
                    except ValueError:
                        raise errors.SetupError('invalid line %i in file %s'
                                                % (line_no, config_file) )

                    continue

                if line_no < 5 and imcon > 0:
                    try:
                        x, y, z = line.split()
                        cell.extend( (float(x), float(y), float(z) ) )

                        x, y, z = crdvel.next().split()
                        cell.extend( (float(x), float(y), float(z) ) )

                        x, y, z = crdvel.next().split()
                        cell.extend( (float(x), float(y), float(z) ) )
                    except ValueError:
                        raise errors.SetupError('invalid cell in file %s'
                                                % config_file)

                    line_no += 3

                    continue

                natoms += 1

                if levcfg >= 0:
                    try:
                        x, y, z = crdvel.next().split()
                        coords.extend( (float(x), float(y), float(z) ) )
                    except ValueError:
                        raise errors.SetupError('invalid coords in file %s, '
                                                'line %s' % (crdvel, line_no) )
                    line_no += 1

                if levcfg >= 1:
                    try:
                        x, y, z = crdvel.next().split()
                        vels.extend( (float(x), float(y), float(z) ) )
                    except ValueError:
                        raise errors.SetupError('invalid vels in file %s, '
                                                'line %s' % (crdvel, line_no) )
                    line_no += 1

                if levcfg >= 2:
                    crdvel.next().split()
                    line_no += 1

        if vels:
            for i in range(natoms * 3):
                vels[i] = float(vels[i]) / const.AMBER_VELCONV
        else:
            vels = [0.0] * natoms * 3

        # FIXME: non-orthorombic box
        la = math.sqrt(cell[0]**2 + cell[1]**2 + cell[2]**2)
        lb = math.sqrt(cell[3]**2 + cell[4]**2 + cell[5]**2)
        lc = math.sqrt(cell[6]**2 + cell[7]**2 + cell[8]**2)

        # FIXME: may fail
        self.dlpoly.unwrap(coords, la, lb, lc)
        self.sander_crd = self._write_rst7(natoms, la, lb, lc, coords, vels,
                                           True)

    def _self_check(self, mdprog):
        """
        Check DL_POLY installation.

        :param mdprog: the dl_poly executable provided to the class
        :type mdprog: str
        """

        if not 'DLPOLYHOME' in os.environ:
            raise errors.SetupError('DLPOLYHOME not set')

        self.mdprog = os.path.join(os.environ['DLPOLYHOME'], 'execute',
                                   mdprog)

        if not os.access(self.mdprog, os.X_OK):
            raise errors.SetupError('DLPOLYHOME does not have a %s binary' %
                                    mdprog)


PROTOCOLS = dict(
    # there appears to be no control of the CG algorithm, run as 1 step
    MIN_STD = '''=== minimization ===
no topology
no strict
optimise distance 0.01
temperature 1.0
steps 1
equil 1
print 1
stats 1
timestep 0.0005
rcut 8.0
delr 1.0
ewald precision 1e-6
mxshak 250
shake 1.0e-5
job time 9999.0
close time 10.0
finish
''',

    MD_SHRINK = '''=== shrink box ===
no topology
no strict

temperature {0}
press {5}
ensemble npt berendsen 0.5 0.1

steps {1}
equil {1}
print {2}
stats {3}
traj 0 {4} 0

timestep {6}
rcut 8.0
delr 1.0

ewald precision 1e-6

mxshak 500
shake 1.0e-5

job time 999999.0
close time 10.0

finish
''',

    # does pseudo do something similar?
    MD_RANDOMV = '''=== Pseudo ensemble ===
no topology
no strict

temperature {0}
ensemble nvt berendsen 1.0
pseudo langevin 2.0 {0}

steps {1}
equil {1}
print {2}
stats {3}
traj 0 {4} 0

timestep {6}
rcut 8.0
delr 1.0

ewald precision 1e-6

mxshak 500
shake 1.0e-5

job time 999999.0
close time 10.0

finish
''',

    MD_HEAT = '''=== heating ===
no topology
no strict

temperature {0}
ensemble nvt berendsen 1.0

steps {1}
equil {1}
print {2}
stats {3}
traj 0 {4} 0

timestep {6}
rcut 8.0
delr 1.0

ewald precision 1e-6

mxshak 500
shake 1.0e-5

job time 999999.0
close time 10.0

finish
''',

    MD_CONSTT = '''=== constant temperature ===
no topology
no strict

temperature {0}
ensemble nvt berendsen 1.0

steps {1}
equil {1}
print {2}
stats {3}
traj 0 {4} 0

timestep {6}
rcut 8.0
delr 1.0

ewald precision 1e-6

mxshak 500
shake 1.0e-5

job time 999999.0
close time 10.0

finish
''',

    MD_PRESS = '''=== constant pressure and temperature ===
no topology
no strict

temperature {0}
press {5}
ensemble npt berendsen 1.0 0.1

steps {1}
equil {1}
print {2}
stats {3}
traj 0 {4} 0

timestep {6}
rcut 8.0
delr 1.0

ewald precision 1e-6

mxshak 500
shake 1.0e-5

job time 999999.0
close time 10.0

finish
'''
)
