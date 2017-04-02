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
Gromacs MD engine: sander or pmemd
"""

__revision__ = "$Id$"



import os
import glob
import re

import mdebase
from FESetup import const, errors, logger
from FESetup.prepare.amber import gromacs, utils


# assume standard GROMACS file name conventions
GROMACS5_EXE_NAMES = ['gmx_d', 'gmx', 'gmx_mpi_d', 'gmx_mpi']
GROMACS_SUFFIXES = ['_d', '', '_mpi_d', '_mpi']


def _check_exe(bindir, exe_name):
    full_path = os.path.join(bindir, exe_name)

    if os.access(full_path, os.X_OK):
        return full_path

    for suffix in GROMACS_SUFFIXES:
        if os.access(full_path, os.X_OK):
            return full_path

    return None

def _get_suffix(name, suffix):
    for i, char in enumerate(name):
        if char == suffix:
            return name[i:]

    return ''


class MDEngine(mdebase.MDEBase):
    """
    Gromacs MD engine.
    """

    def __init__(self, amber_top, amber_crd, sander_crd, sander_rst,
                 amber_pdb, box_dims=None, solvent=None, mdprog='mdrun',
                 mdpref='', mdpost=''):

        super(MDEngine, self).__init__()

        self.update_files(amber_top, amber_crd, sander_crd, sander_rst,
                          amber_pdb)

        self.mdpref = mdpref
        self.mdpost = mdpost

        self.prev = ''
        self.prefix = ''

        self.mdprog = ''
        self.grompp = ''
        self.gmxdump = ''

        self._self_check(mdprog)


    def update_files(self, amber_top, amber_crd, sander_crd, sander_rst,
                     amber_pdb):
        self.amber_top = amber_top
        self.amber_crd = amber_crd
        self.sander_crd = sander_crd
        self.sander_rst = sander_rst
        self.amber_pdb = amber_pdb

        self.top = os.path.splitext(amber_top)[0] + os.extsep + 'top'
        self.gro = os.path.splitext(amber_crd)[0] + os.extsep + 'gro'

        gtop = gromacs.GromacsTop()

        gtop.readParm(amber_top, amber_crd)
        self.molidx = gtop.molidx

        gtop.writeTop(self.top, '', '', False)
        gtop.writeGro(self.gro)

        self.gtop = gtop


    def minimize(self, config='%STD', nsteps=100, ncyc=100, mask='',
                 restr_force=5.0):
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

        prefix = mdebase.MIN_PREFIX + '%05i' % self.run_no

        if config[0] == '%':
            pname = config[1:]
            logger.write('Running minimisation with protocol %s' % pname)

            try:
                config = PROTOCOLS['MIN_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such min protocol predefined: '
                                        '%s' % config)

            if mask:
                cons = 'define = -DPOSRES'
            else:
                cons = ''

            # read coordinates from self.amber_pdb because origin is at center
            # of cell, it is on one edge in the parm file!
            #
            # FIXME: water model, box shape
            config = config.format(cons, nsteps, ncyc, nsteps / 10, nsteps / 5)

        self._run_mdprog(prefix, config, mask, restr_force)
        self.prefix = prefix


    def md(self, config='', nsteps=1000, T=300.0, p=1.0, mask='',
           restr_force=5.0, nrel=1, wrap=True, dt=0.002):
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

        prefix = mdebase.MD_PREFIX + '%05i' % self.run_no

        if self.run_no == 1:
            gen_vel = 'yes'
        else:
            gen_vel = 'no'

        if config[0] == '%':
            pname = config[1:]
            logger.write('Running MD with protocol %s' % pname)

            try:
                config = PROTOCOLS['MD_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such MD protocol predefined: '
                                        '%s' % config)

            if mask:
                cons = 'define = -DPOSRES'
            else:
                cons = ''

            # read coordinates from self.amber_pdb because origin is at center
            # of cell, it is on one edge in the parm file!
            #
            # FIXME: water model, box shape
            config = config.format(cons, nsteps, nsteps / 10, nsteps / 5,
                                   dt, nsteps * dt, gen_vel, T, p)

        self._run_mdprog(prefix, config, mask, restr_force)
        self.prefix = prefix


    def get_box_dims(self):
        """
        Extract box information from rst7 file.

        :returns: box dimensions
        """

        gro_file = self.prefix + os.extsep + 'gro'

        with open(gro_file, 'r') as rst:
            for line in rst:
                last_line = line

        # FIXME: rectangular box only
        box_dims = [float(d) / const.A2NM for d in last_line.split()]
        box_dims.extend( (90.0, 90.0, 90.0) )

        return box_dims


    def _run_mdprog(self, prefix, config, mask, restr_force):
        """
        Run mdrun from the Gromacs package.
        """

        filename = prefix + os.extsep
        config_filename = '_' + filename + 'mdp'

        with open(config_filename, 'w') as mdin:
            mdin.writelines(config)

        if self.run_no == 1:
            params = ('-f %s -c %s -r %s -p %s -o %s -po %s' %
                      (config_filename, self.gro, self.gro, self.top,
                       filename + 'tpr', filename + 'mdp') )
        else:
            edr = self.prev + os.extsep + 'edr'
            trr = self.prev + os.extsep + 'trr'

            params = ('-f %s -c %s -r %s -e %s -t %s -p %s -o %s -po %s' %
                      (config_filename, self.prev + os.extsep + 'gro',
                       self.gro, edr, trr, self.top, filename + 'tpr',
                       filename + 'mdp'))

        if mask:
             self._make_restraints(mask, restr_force)

        retc, out, err = utils.run_exe(' '.join((self.grompp, params)))

        if retc:
            logger.write(err)
            raise errors.SetupError('%s has failed (see logfile)' %
                                    self.grompp)

        params = '-deffnm %s' % prefix

        retc, out, err = utils.run_exe(' '.join((self.mdpref, self.mdprog,
                                                 self.mdpost, params)))

        if retc:
            logger.write(err)
            raise errors.SetupError('%s has failed (see logfile) has failed' %
                                    self.mdprog)

        self.run_no += 1
        self.prev = prefix


    def _make_restraints(self, restr, k):

        k *= const.CALA2JNM2

        try:
            mask = mdebase._restraint_table[restr]
        except KeyError:
            mask = restr

        mask_idx = list(self.mask_indexes(self.amber_top, mask) )

        # FIXME: each molecule type only recorded once, so this means that
        #        the mask is assumed to be the same for each type
        for name, data in self.molidx.iteritems():
            idx_list = data[0]
            natoms = data[1]
            molt_idx = []

            # FIXME: only checks first occurence of molecule in mask index
            for i in idx_list[:natoms]:
                if i in mask_idx:
                    molt_idx.append(i)

            posres_file = (const.GROMACS_POSRES_PREFIX + name +
                           const.GROMACS_ITP_EXT)

            if molt_idx:
                mi = min(molt_idx)

                with open(posres_file, 'w') as posres:
                    posres.write('[ position_restraints ]\n')

                    for i in molt_idx:
                        posres.write('%i 1 %.2f %.2f %.2f\n' %
                                     (i + 1 - mi, k, k, k) )
            else:
                open(posres_file, 'w').close()


    def to_rst7(self):
        """
        Use gmxdump to extract coordinates, velocities and box dimensions
        from either the trajectory and convert to AMBER ASCII .rst7

        :raises: SetupError
        """

        params = '-f %s' % (self.prev + os.extsep + 'trr')
        retc, out, err = utils.run_exe(' '.join((self.gmxdump, params)))

        if retc:
            logger.write(err)
            raise errors.SetupError('%s has failed (see logfile) has failed' %
                                    self.gmxdump)

        lidx = out.rfind('frame ')
        box = []
        coords = []
        vels = []

        # FIXME: check for consistency
        for line in out[lidx:].rstrip().split('\n'):
            if 'natoms= ' in line:
                natoms = line.split()[1]

                continue

            m = re.match('      box\[.*\]={(.*)}', line)

            if m:
                t = m.group(1).split(',')
                box.extend(t)

            m = re.match('      x\[.*\]={(.*)}', line)

            if m:
                t = m.group(1).split(',')
                coords.extend(t)

            m = re.match('      v\[.*\]={(.*)}', line)

            if m:
                t = m.group(1).split(',')
                vels.extend(t)

        natoms = int(natoms)

        for i in range(natoms * 3):
            coords[i] = float(coords[i]) / const.A2NM

        # Gromacs stores velocities in nm/ps, Amber is A/time unit where
        # time unit is 1/20.455 ps
        if vels:
            vel_conv = const.AMBER_VELCONV / 10.0

            for i in range(natoms * 3):
                vels[i] = float(vels[i]) / vel_conv
        else:
            vels = [0.0] * natoms * 3

        xx = float(box[0]) / const.A2NM
        yy = float(box[4]) / const.A2NM
        zz = float(box[8]) / const.A2NM

        # FIXME: may fail
        self.gtop.unwrap(coords, xx, yy, zz)
        self.sander_crd = self._write_rst7(natoms, xx, yy, zz, coords, vels,
                                           False)

    def _self_check(self, mdprog):
        """
        Check which Gromacs version is in GMXHOME and determine what
        executables to use.

        :param mdprog: the mdrun executable provided to the class
        :type mdprog: str
        """

        if not 'GMXHOME' in os.environ:
            raise errors.SetupError('GMXHOME not set')

        bindir = os.path.join(os.environ['GMXHOME'], 'bin')

        if not os.path.isdir(bindir):
            raise errors.SetupError('GMXHOME does not have a bin directory')

        version = 0
        suffixes = []

        for exe_name in GROMACS5_EXE_NAMES:
            full_path = os.path.join(bindir, exe_name)

            if os.access(full_path, os.X_OK):
                version = 5
                suffixes.append(_get_suffix(exe_name, '_'))

        if version == 5:
            prefix = mdprog.split()[0]

            if prefix[:3] != 'gmx':
                suffix = _get_suffix(mdprog, '_')

                # FIXME: proper exception...
                if suffix not in GROMACS_SUFFIXES:
                    raise errors.SetupError('only the following GROMACS '
                                            'suffixes are supported: ' +
                                            ', '.join(GROMACS_SUFFIXES))

                if suffix in suffixes:
                    prefix = 'gmx' + suffix
                else:
                    raise errors.SetupError('gmx%s does not exist in GMXHOME' %
                                            suffix)

                self.mdprog = os.path.join(bindir, prefix, 'mdrun')

            full_path = os.path.join(bindir, prefix)

            self.mdprog = ' '.join((full_path, 'mdrun'))
            self.grompp = ' '.join((full_path, 'grompp'))
            self.gmxdump = ' '.join((full_path, 'dump'))

            return
        else:
            suffix = _get_suffix(mdprog, '_')

            if suffix not in GROMACS_SUFFIXES:
                raise errors.SetupError('only the following GROMACS '
                                        'suffixes are supported: ' +
                                        ', '.join(GROMACS_SUFFIXES))

            # FIXME: check if suffixed version actually exists
            #        if not use the executables that are found
            self.grompp = _check_exe(bindir, 'grompp' + suffix)
            self.gmxdump = _check_exe(bindir, 'gmxdump' + suffix)

            full_path = os.path.join(bindir, mdprog)

            # special case of user supplied GROMACS4 mdrun_suffix
            if os.access(full_path, os.X_OK):
                self.mdprog = full_path

            if self.grompp and self.mdprog and self.gmxdump:
                return

        raise errors.SetupError('GMXHOME does not have any useful executables')



NB_PARAMS = '''\
cutoff-scheme        = Verlet
nstcalclr            = 1
rlist                = 0.8
nstlist              = 10
ns_type              = grid
pbc                  = xyz

coulombtype          = PME
coulomb-modifier     = None
rcoulomb             = 0.8
fourierspacing       = 0.12
pme_order            = 4
ewald_rtol           = 1.0E-5
optimize_fft         = yes

vdwtype              = cutoff
vdw-modifier         = None
rvdw                 = 0.8
DispCorr             = AllEnerPres'''

CONS_PARAMS = '''\
constraints          = h-bonds
constraint_algorithm = Lincs
continuation         = no
lincs_order          = 4
lincs_warnangle      = 30'''

PROTOCOLS = dict(
    # Gromacs appears to be sensitive to this and minimisation should be run
    # first, it also looks like this should be done through the
    # double precision executable
    MIN_STD = '''; minimization
{0}
integrator            = steep
nsteps                = {1}
nstcgsteep            = {2}
nstlog                = {3}
nstxout               = {4}
emtol                 = 10.0
emstep                = 0.01
%s

%s
''' % (NB_PARAMS, CONS_PARAMS),

    MD_SHRINK = '''; shrink box
{0}
integrator           = md
tinit                = 0.0
dt                   = {4}
nsteps               = {1}
nstcomm              = 1000

nstxout              = {2}
nstvout              = {2}
nstfout              = 0
nstlog               = {3}
nstenergy            = {2}
nstxtcout            = 0

%s

tcoupl               = Berendsen
tc_grps              = System
tau_t                = 0.5
ref_t                = {7}

pcoupl               = Berendsen
pcoupltype           = isotropic
tau_p                = 0.2
compressibility      = 4.5e-5
ref_p                = {8}
refcoord-scaling     = com

gen_vel              = {6}
gen_temp             = {7}
gen_seed             = 173529

%s
''' % (NB_PARAMS, CONS_PARAMS),

    MD_RANDOMV = '''; Maxwell velocities
{0}
integrator            = md
tinit                 = 0.0
dt                    = {4}
nsteps                = {1}
nstcomm               = 1000

nstxout               = {2}
nstvout               = {2}
nstfout               = 0
nstlog                = {3}
nstenergy             = {2}
nstxtcout             = 0

%s

tcoupl                = Berendsen
tc_grps               = System
tau_t                 = 1.0
ref_t                 = {7}
annealing             = single
annealing-npoints     = 2
annealing-time        = 0.0 {5}
annealing-temp        = 5.0 300.0

pcoupl                = no

gen_vel               = yes
gen_temp              = {7}
gen_seed              = 173529

%s
''' % (NB_PARAMS, CONS_PARAMS),

    MD_HEAT = '''; heating
{0}
integrator           = md
tinit                = 0.0
dt                   = {4}
nsteps               = {1}
nstcomm              = 1000

nstxout              = {2}
nstvout              = {2}
nstfout              = 0
nstlog               = {3}
nstenergy            = {2}
nstxtcout            = 0

%s

tcoupl               = Berendsen
tc_grps              = System
tau_t                = 1.0
ref_t                = {7}
annealing            = single
annealing-npoints    = 2
annealing-time       = 0.0 {5}
annealing-temp       = 5.0 300.0

pcoupl               = no

gen_vel              = {6}
gen_temp             = 5.0
gen_seed             = 173529

%s
''' % (NB_PARAMS, CONS_PARAMS),

    MD_CONSTT = '''; constant temperature
{0}
integrator           = md
tinit                = 0.0
dt                   = {4}
nsteps               = {1}
nstcomm              = 1000

nstxout              = {2}
nstvout              = {2}
nstfout              = 0
nstlog               = {3}
nstenergy            = {2}
nstxtcout            = 0

%s

tcoupl               = Berendsen
tc_grps              = System
tau_t                = 1.0
ref_t                = {7}

pcoupl               = no

gen_vel              = {6}
gen_temp             = {7}
gen_seed             = 173529

%s
''' % (NB_PARAMS, CONS_PARAMS),

    # goes faster if temperature is higher...
    MD_PRESS = '''; constant pressure and temperature
{0}
integrator           = md
tinit                = 0.0
dt                   = {4}
nsteps               = {1}
nstcomm              = 1000

nstxout              = {2}
nstvout              = {2}
nstfout              = 0
nstlog               = {3}
nstenergy            = {2}
nstxtcout            = 0

%s

tcoupl               = Berendsen
tc_grps              = System
tau_t                = 1.0
ref_t                = {7}

pcoupl               = Berendsen
pcoupltype           = isotropic
tau_p                = 0.2
compressibility      = 4.5e-5
ref_p                = {8}
refcoord-scaling     = com

gen_vel              = {6}
gen_temp             = {7}
gen_seed             = 173529

%s
''' % (NB_PARAMS, CONS_PARAMS)
)
