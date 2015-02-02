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
Gromacs MD engine: sander or pmemd
"""

__revision__ = "$Id$"



import os, re

import mdebase
from FESetup import const, errors, logger
from FESetup.prepare.amber import gromacs, utils



class MDEngine(mdebase.MDEBase):
    """
    Gromacs MD engine.
    """

    def __init__(self, amber_top, amber_crd, sander_crd, sander_rst,
                 amber_pdb, box_dims = None, solvent = None, mdprog = 'mdrun',
                 mdpref = '', mdpost = ''):

        super(MDEngine, self).__init__()

        self.update_files(amber_top, amber_crd, sander_crd, sander_rst,
                          amber_pdb)

        self.mdprog = mdprog
        self.mdpref = mdpref
        self.mdpost = mdpost

        self.prev = ''


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
            params = ('-f %s -c %s -r %s -e %s -t %s -p %s -o %s -po %s' %
                      (config_filename, self.prev + os.extsep + 'gro', self.gro,
                       self.prev + os.extsep + 'edr',
                       self.prev + os.extsep + 'trr',
                       self.top, filename + 'tpr', filename + 'mdp') )

        if mask:
             self._make_restraints(mask, restr_force)

        out, err = utils.run_gromacs('grompp', '', '', params)

        if out == 1:
            logger.write(err)
            raise errors.SetupError('grompp has failed')

        params = '-deffnm %s' % prefix

        out, err = utils.run_gromacs(self.mdprog, self.mdpref, self.mdpost, params)

        if out == 1:
            logger.write(err)
            raise errors.SetupError('Gromacs MD program has failed')

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
        out, err = utils.run_gromacs('gmxdump', '', '', params)

        if out == 1:
            logger.write(err)
            raise errors.SetupError('Gromacs gmxdump has failed')

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
        self._write_rst7(natoms, xx, yy, zz, coords, vels, False)


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
    # Gromacs appears to be sensitive to this and minimisation should be run first
    MIN_STD = '''; minimization
{0}
integrator            = cg
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
