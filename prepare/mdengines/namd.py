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
NAMD MD engine
"""

__revision__ = "$Id$"



import os, sys, struct

import mdebase
from FESetup import const, errors, logger
from FESetup.prepare.amber import utils


SOLVENT_TRANS = {
    'tip3p': 'tip3',
    'tip4pew': 'tip4',
    'spce': 'tip3'
    }

STEPS_PER_CYCLE = 20


def namd_velcoor(filename):
    """
    Return number of atoms and coordinates or velocities from NAMD binary
    coor or vel file.

    :param filename: file name of NAMD coordinates or velocities
    :type filename: string
    :raises: SetupError
    :returns: number of atoms, coordinates or velocities
    :rtype: integer, list
    """

    natoms_expected = (os.stat(filename).st_size - 4) / 24

    endian = '<'

    with open(filename, 'rb') as coor:
        first4 = coor.read(4)
        natoms = struct.unpack('%si' % endian, first4)[0]

        if natoms != natoms_expected:
            endian = '>'
            natoms = struct.unpack('%si' % endian, first4)[0]

            if natoms != natoms_expected:
                raise errors.SetupError('BUG: unknown endianess in %s' %
                                        filename)

        n = 3 * natoms
        coords = struct.unpack('%s%id' % (endian, n), coor.read(n * 8) )

    return natoms, coords


# FIXME: write individual steps into TCL script and run only once
class MDEngine(mdebase.MDEBase):
    """
    NAM MD engine.
    """

    # FIXME: files are specific to AMBER but needed for conversion for
    #        other MD packages
    def __init__(self, amber_top, amber_crd, sander_crd, sander_rst,
                 amber_pdb, box_dims=[0.0, 0.0, 0.0], solvent='tip3',
                 mdprog='namd2', mdpref='', mdpost=''):

        super(MDEngine, self).__init__()

        # FIXME: do not assume that top/crd are in AMBER format
        self.update_files(amber_top, amber_crd, sander_crd, sander_rst,
                          amber_pdb, box_dims, solvent)

        self.mdpref = mdpref
        self.mdpost = mdpost

        self.prev = ''
        self.prefix = ''

        self.mdprog = ''
        self._self_check(mdprog)


    def update_files(self, amber_top, amber_crd, sander_crd, sander_rst,
                     amber_pdb, box_dims, solvent):
        self.amber_top = amber_top
        self.amber_crd = amber_crd
        self.sander_crd = sander_crd
        self.sander_rst = sander_rst
        self.amber_pdb = amber_pdb

        try:
            self.solvent = SOLVENT_TRANS[solvent]
        except KeyError:
            raise errors.SetupError('namd does not know about solvent model %s'
                                    % solvent)

        # FIXME: only cuboid box
        self.xx, self.yy, self.zz = box_dims[0:3]


    def minimize(self, config='%STD', nsteps=100, ncyc=10, mask='',
                 restr_force=5.0):
        """
        Minimize a system.

        :param config: MD paramers, if start with '%' respective default is
           chosen, otherwise a full sander namelist as string
        :type config: string
        :param nsteps: maximum number of conjugated gradient steps
        :type nsteps: integer
        :param ncyc: number of initial steepest decent steps (unused)
        :type ncyc: float
        :param mask: pre-defined restraints or AMBER mask
        :type mask: string
        :param restr_force: force constant for the restraints in kcal/mol/AA
        :type restr_force: float
        :param box_dims: box dimensions
        :type box_dims: list
        :raises: SetupError
        """

        prefix = mdebase.MIN_PREFIX + '%05i' % self.run_no
        restr_filen = prefix + '_restr' + const.PDB_EXT

        if mask:
            self._make_restraints(restr_filen, mask, restr_force)

        if config[0] == '%':
            pname = config[1:]
            logger.write('Running minimisation with protocol %s' % pname)

            try:
                config = PROTOCOLS['AMBER_MIN_%s' % pname]
            except KeyError:
                raise errors.SetupError('no such min protocol predefined: '
                                        '%s' % config)

            if mask:
                cons = 'on'
            else:
                cons = 'off'

            # read coordinates from self.amber_pdb because origin is at center
            # of cell, it is on one edge in the parm file!
            #
            # FIXME: box shape
            config = config.format(self.prev, nsteps, self.amber_top,
                                   self.amber_pdb, self.solvent, self.xx,
                                   self.yy, self.zz, prefix, cons, restr_filen,
                                   STEPS_PER_CYCLE)

        self._run_mdprog(prefix, config)
        self.prefix = prefix


    def md(self, config='%STD', nsteps=1000, T=300.0, p=1.0, mask='',
           restr_force=5.0, nrel=1, wrap=True, dt=0.002):
        """
        Run molecular dynamics on a system.

        :param namelist: MD paramers, if start with '%' respective default is
           chosen, otherwise a full sander namelist as string
        :type namelist: string
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

        dt = dt * 1000.0
        prefix = mdebase.MD_PREFIX + '%05i' % self.run_no
        
        restr_filen = prefix + '_restr' + const.PDB_EXT

        if mask:
            self._make_restraints(restr_filen, mask, restr_force)
            restb = 'T'
        else:
            restb = 'F'

        if nsteps % STEPS_PER_CYCLE:
            nsteps = (nsteps / STEPS_PER_CYCLE + 1) * STEPS_PER_CYCLE
            logger.write('Warning: nsteps has been changed to %i to ensure it '
                         'to be multiple of stepsPerCycle' % nsteps)

        if config[0] == '%':
            pname = config[1:]
            logger.write('Running MD with protocol %s' % pname)

            try:
                job_type = pname.lower()
                config = PROTOCOLS['AMBER_MD_STD']
            except KeyError:
                raise errors.SetupError('no such MD protocol predefined: '
                                        '%s' % namelist)

            if wrap:
                wrap = 'on'
            else:
                wrap = 'off'

            # FIXME: box shape
            config = config.format(self.prev, prefix, job_type, T, p, nsteps,
                                   self.amber_top, self.amber_pdb,
                                   self.solvent, self.xx, self.yy, self.zz,
                                   wrap, nrel, restr_filen, restb,
                                   STEPS_PER_CYCLE, dt)

        self._run_mdprog(prefix, config)
        self.prefix = prefix


    def get_box_dims(self):
        """
        Extract box information from rst7 file.

        :returns: box dimensions
        """

        xst_file = self.prefix + os.extsep + 'xst'

        with open(xst_file, 'r') as rst:
            for line in rst:
                last_line = line

        # FIXME: rectangular box only
        d = last_line.split()
        box_dims = [float(d[1]), float(d[5]), float(d[9]), 90.0, 90.0, 90.0]
                        
        return box_dims


    def to_rst7(self):
        """
        Convert binary NAMD coordinates and velocities plus extended system
        information (xsc) into AMBER ASCII .rst7

        :raises: SetupError
        """
        
        prefix = self.prev + os.extsep

        natoms, coords = namd_velcoor(prefix + 'coor')
        ncheck, vels = namd_velcoor(prefix + 'vel')

        if natoms != ncheck:
            raise errors.SetupError('different number of atoms in coor(%i) '
                                    'and vel(%i) files' % (natoms, ncheck) )

        with open(prefix + 'xsc') as xsc:
            for line in xsc:
                ext = line

        ext = ext.split()

        # FIXME: only cuboid box
        xx, yy, zz = float(ext[1]), float(ext[5]), float(ext[9])

        self.sander_crd = self._write_rst7(natoms, xx, yy, zz, coords, vels,
                                           True)


    def _run_mdprog(self, prefix, config):
        """
        Run namd.
        """

        filename = prefix + os.extsep
        config_filename = filename + 'in'
        
        with open(config_filename, 'w') as mdin:
            mdin.writelines(config)

        retc, out, err = utils.run_exe(' '.join((self.mdpref, self.mdprog,
                                                 self.mdpost,
                                                 config_filename)))

        with open(filename + 'out', 'w') as outfile:
            outfile.writelines(out)

        if retc:
            logger.write(err)
            raise errors.SetupError('%s has failed (see logfile)' %
                                    self.mdprog)

        self.run_no += 1
        self.prev = prefix


    def _make_restraints(self, ofilen, restr, k):

        try:
            mask = mdebase._restraint_table[restr]
        except KeyError:
            mask = restr

        # FIXME: assumes AMBER parmtop
        indexes = list(self.mask_indexes(self.amber_top, mask) )
        acnt = 0

        with open(ofilen, 'w') as opdb:
            with open(self.amber_pdb, 'r') as ipdb:
                for line in ipdb:
                    if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
                        if acnt in indexes:
                            tmp = line[:60] + '%6.2f' % k + line[66:]
                        else:
                            tmp = line

                        acnt += 1
                    else:
                        tmp = line

                    opdb.write(tmp)


    def _self_check(self, mdprog):
        """
        Check NAMD installation.

        :param mdprog: the namd executable provided to the class
        :type mdprog: str
        """

        if not 'NAMDHOME' in os.environ:
            raise errors.SetupError('NAMDHOME not set')

        self.mdprog = os.path.join(os.environ['NAMDHOME'], mdprog)

        if not os.access(self.mdprog, os.X_OK):
            raise errors.SetupError('NAMDHOME does not have a %s binary' %
                                    mdprog)
        
        

PROTOCOLS = dict(
    # NAMD also appears to need prior minimisation, possibly because of tight
    # pressure coupling
    AMBER_MIN_STD = '''
set prev_step     "{0}"

minimization      on
minTinyStep       1.0E-3
minBabyStep       1.0E-7

set nsteps        {1}
numsteps          $nsteps

amber             yes
parmfile          "{2}"
coordinates       "{3}"
readexclusions    yes
exclude           scaled1-4
1-4scaling        0.833333
scnb              2.0
zeromomentum      on
LJcorrection      on

watermodel        {4}
useSettle         on
rigidBonds        all
rigidIterations   300
rigidTolerance    1.0e-8
rigidDieOnError   off

PME               on
PMEGridSpacing    1.0
nonbondedFreq     1
fullElectFrequency 1

switching         off
cutoff            8.0
pairlistsPerCycle 1
stepspercycle     {11}

if {{$prev_step != ""}} {{
  binvelocities       "$prev_step.vel"
  bincoordinates      "$prev_step.coor"
  ExtendedSystem      "$prev_step.xsc"
  firsttimestep       0
}} else {{
  cellBasisVector1  {5}  0   0
  cellBasisVector2   0  {6}  0
  cellBasisVector3   0   0  {7}
}}

outputname        "{8}"
outputEnergies    [expr $nsteps / 5]
outputTiming      $nsteps
binaryoutput      yes

restartname       "{8}"
restartfreq       [expr $nsteps / 10]
binaryrestart     yes

XSTfile           "{8}.xst"
XSTfreq           [expr $nsteps / 5]

constraints       {9}
consexp           2
consref           "{10}"
conskfile         "{10}"
conskcol          B
''',

    AMBER_MD_STD = '''
set prev_step     "{0}"
set job_type      "{2}"

set restraint     {15}

set startT        5.0
set finalT        {3}
set p             {4}
set nsteps        {5}

amber             yes
parmfile          "{6}"
coordinates       "{7}"
readexclusions    yes
exclude           scaled1-4
1-4scaling        0.833333
scnb              2.0
zeromomentum      on
LJcorrection      on

watermodel        {8}
useSettle         on
rigidBonds        all
rigidIterations   300
rigidTolerance    1.0e-8
rigidDieOnError   off

PME                on
PMEGridSpacing     1.0
nonbondedFreq      1
fullElectFrequency 1

switching         off
cutoff            8.0
pairlistsPerCycle 1
stepspercycle     {16}

langevin          on
langevinHydrogen  on
langevinDamping   1

if {{$job_type == "press" || $job_type == "shrink" || $job_type == "relres"}} {{
  margin               2.5
  BerendsenPressure  on
  BerendsenPressureTarget $p
  BerendsenPressureCompressibility 4.5E-5
  BerendsenPressureRelaxationTime 50.0
}} else {{
  useGroupPressure  yes
  LangevinPiston    off
  BerendsenPressure off
}}

outputname       "{1}"
outputEnergies   [expr $nsteps / 10]
outputPressure   [expr $nsteps / 10]
outputTiming     $nsteps
binaryoutput     yes

restartname      "{1}"
restartfreq      [expr $nsteps / 5]
binaryrestart    yes

XSTfile          "{1}.xst"
XSTfreq          [expr $nsteps / 10]

wrapAll          {12}
wrapNearest      on
DCDfile          "{1}.dcd"
DCDUnitCell      yes
DCDfreq          [expr $nsteps / 5]

if {{$restraint == "T"}} {{
  constraints       on
  consexp           2
  consref           "{14}"
  conskfile         "{14}"
  conskcol          B
}}

if {{$prev_step != ""}} {{
  binvelocities       "$prev_step.vel"
  bincoordinates      "$prev_step.coor"
  ExtendedSystem      "$prev_step.xsc"
  firsttimestep       0
}} elseif {{$job_type == "randomv"}} {{
  cellBasisVector1    {9}  0  0
  cellBasisVector2    0  {10} 0
  cellBasisVector3    0  0  {11}
  temperature $finalT
}} else {{
  cellBasisVector1    {9}  0  0
  cellBasisVector2    0  {10} 0
  cellBasisVector3    0  0  {11}
  temperature $startT
}}

timestep {17}


proc check_multiple {{p1 p2}} {{
  if {{[expr $p1 % $p2] != 0}} {{
    set p1 [expr ($p1 / $p2 + 1) * $p2]
  }}

  return $p1
}}

if {{$job_type == "randomv"}} {{
  langevinTemp $finalT
  run $nsteps
}} elseif {{$job_type == "heat"}} {{
  set tinc [check_multiple [expr $nsteps / 10] {16}]

  for {{set step 0}} {{$step < [expr $nsteps - $tinc]}} {{incr step $tinc}} {{
    langevinTemp [expr {{($finalT - $startT) * $step / $nsteps + $startT}} ]
    run $tinc
  }}

  langevinTemp $finalT
  run $tinc
}} elseif {{$job_type == "press" || $job_type == "shrink" || $job_type == "constt"}} {{
  langevinTemp $finalT
  run $nsteps
}} elseif {{$job_type == "relres"}} {{
  langevinTemp $finalT
  set rinc [check_multiple [expr $nsteps / {13:d}] {16}]

  if {{ $restraint == "T"}} {{
    for {{set step 0}} {{$step < $nsteps}} {{incr step $rinc}} {{
      constraintScaling [ expr {{ -double($step) / $nsteps + 1.0 }} ]
      run $rinc
    }}
  }} else {{
    run $nsteps
  }}
}}
'''
)
