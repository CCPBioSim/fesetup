#  Copyright (C) 2012-2016  Hannes H Loeffler, Julien Michel
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
The Common parent class for the free energy setup of a ligand-protein complex
plus various helper functions.

The Common class provides common methods for the setup of a ligand-protein
complex.  Classes for ligand, protein receptor and complex subclasses inherit
from Common.  Common can not be instantiated directly.  Methods are
sequentially coupled, i.e. suitable wrapper code is needed.

Various helper functions to run external programs, etc. are defined too.

Common methods:
    minimize, md, flatten_ring, get_box_dims
"""

__revision__ = "$Id$"



import os, sys, re, shutil
import openbabel as ob
import pybel

import utils                            # relative import
from FESetup import const, errors, logger, report
from leap import Leap

import Sire.IO



def ssbonds(ss_file):
    """
    Read file with SS-bond information.

    :param ss_file: filename with disulfide bond information
    :type ss_file: string

    :returns: pair of SS-bond indices
    """

    lcnt = 0
    pairs = []

    with open(ss_file, 'r') as ssb:
        for line in ssb:
            lcnt += 1
            l = line.lstrip()

            if l.startswith('#') or l == '':
                continue

            try:
                a, b = line.split()
                a, b = int(a), int(b)

                if a == b:
                    raise errors.SetupError('malformed file %s in '
                            'line %i: residue numbers cannot be the '
                            'same' % (lcnt, ss_file) )

                pairs.append( (a, b) )
            except (ValueError, TypeError) as why:
                raise errors.SetupError('malformed file %s in '
                            'line %i: %s' % (ss_file, lcnt, why) )

    return pairs


def read_box_file(box_file):
    """
    Read box dimensions from file.

    :param box_file: filename with box dimensions
    :type box_file: string

    :returns: box dimensions
    """

    boxdata = None

    with open(box_file, 'r') as box:
        for line in box:
            l = line.lstrip()

            if l.startswith('#') or l == '':
                continue

            bb_len = line.split()

            if len(bb_len) == 1:
                try:
                    boxdata = '%s' % float(bb_len[0])
                except ValueError:
                    raise errors.SetupError('boxfile %s contains '
                                            'non-float data' % box_file)
            elif len(bb_len) >= 3:
                # tleap generates rst7 with origin in corner (but PDB
                # has it in the centre???), setbox takes _buffer_
                try:
                    boxdata = '{ %s %s %s }' % (float(bb_len[0]) / 2.0,
                                                float(bb_len[1]) / 2.0,
                                                float(bb_len[2]) / 2.0)
                except ValueError:
                    raise errors.SetupError('boxfile %s contains '
                                            'non-float data' % box_file)
            else:
                raise errors.SetupError('boxfile %s must contain either '
                                        '1 or more than 3 floats' %
                                        box_file)

    return boxdata


# FIXME: check very carefully the flow of method calls.  The current assumption
# is that they are called in a particular order.  Think about work flow: what is
# optional? do we have a well-defined order?  Make sure current work mol2 is the
# right one!

class Common(object):
    """
    Common methods for the setup protocol.  Used through subclasses for
    ligand, protein an complex setup.
    """

    PROT_MAP = {
        'HIS' : 'HIP',
        'ASP' : 'ASH',
        'GLU' : 'GLH',
        'LYS' : 'LYN',
        'CYS' : 'CYM',
        'ARG' : 'ARG',
        'TYR' : 'TYR'
        }

    TOP_EXT = os.extsep + 'parm7'
    RST_EXT = os.extsep + 'rst7'


    def __new__(cls, *args, **kwargs):
        """Prevent direct instantiation."""

        if cls is Common:
            raise TypeError('Common may not be instantiated directly.')

        return object.__new__(cls)


    def __init__(self, mol_name, basedir, workdir, overwrite):
        """
        :param mol_name: (file) name of molecules
        :type mol_name: string
        :param basedir: the base directory, i.e. location of initial data files
        :type basedir: string
        :param workdir: the working directory
        :type workdir: string
        :param overwrite: overwrite files in the working directory from basedir
        :type overwrite: string
        :raises: SetupError
        """

        self.mol_file = ''
        self.mol_fmt = ''

        self.sander_rst = ''
        self.sander_crd = ''
        self.namd_prefix = ''

        self.min_no = 0                 # number of current minimisation step
        self.md_no = 0                  # number of current MD step

        self.charge = 0.0
        self.box_dims = []

        self.amber_top = ''             # only set in _amber_top_common
        self.amber_crd = ''             # only set in _amber_top_common
        self.amber_pdb = ''

        # set in __init__.py before this __init__(): do not set here or anywhere
        # else as they are meant to be constants for the whole class hierarchy!
        #self.ff_cmd
        #self.solvent
        #self.solvent_load 
        #self.solvent_box 
        #self.MDEngine

        self.leap = Leap(self.force_fields, self.solvent_load)

        self.mdengine = None

        self.mol_name = mol_name

        self.topdir = os.getcwd()

        self.dst = os.path.join(self.topdir, workdir, self.mol_name)
        self.rel_dst = os.path.join(workdir, self.mol_name)

        self._parm_overwrite = None

        if not basedir:
            return

        self.basedir = basedir

        logger.write('*** Working on %s ***\n' % self.mol_name)

        src = os.path.join(basedir, self.mol_name)

        try:
            # FXIME: catch exceptions due to failed file actions
            if not os.access(self.dst, os.F_OK):
                os.makedirs(self.dst)

                logger.write('Copying directory contents of %s to %s' %
                             (src, self.dst) )

                for filename in os.listdir(src):
                    shutil.copy(os.path.join(src, filename), self.dst)
            elif overwrite:
                logger.write('Overwrite mode: copying directory contents of %s '
                             'to %s' % (src, self.dst) )

                for filename in os.listdir(src):
                    shutil.copy(os.path.join(src, filename), self.dst)


        except OSError as why:
            raise errors.SetupError(why)

        self.model = None


    # context manager used to keep track of directory changes
    def __enter__(self):
        '''Enter directory dst.'''
        
        logger.write('Entering %s' % self.dst)
        os.chdir(self.dst)

        return self


    def __exit__(self, typ, value, traceback):
        '''Leave directory dst and return to topdir.'''

        logger.write('Entering %s\n' % self.topdir)
        os.chdir(self.topdir)

        return


    def _amber_top_common(self, boxtype = '', boxlength = '10.0',
                          boxfile = None, neutralize = False,
                          align=None, remove_first = False):
        """Common scripting commands for leap.  Internal function only."""

        leapin = self.leap.generate_init()

        if os.access(const.SSBOND_FILE, os.R_OK):
            pairs = ssbonds(const.SSBOND_FILE)
            cmd = []

            for a, b in pairs:
                cmd.append('bond s.%i.SG s.%i.SG\n' % (a, b) )

            leapin += ''.join(cmd)

        boxdata = None

        if align:
            leapin += 'alignAxes s\n'

        if boxfile:
            boxdata = read_box_file(boxfile)

        # bug #976: input may contain ions and water
        #leapin += self.solvent_load

        if boxtype:
            self.amber_top = const.LEAP_SOLVATED + self.TOP_EXT
            self.amber_crd = const.LEAP_SOLVATED + self.RST_EXT
            self.amber_pdb = const.LEAP_SOLVATED + const.PDB_EXT

            if boxtype == 'octahedron':
                leapin += ('solvateOct s %s %s 0.75\n' %
                           (self.solvent_box, boxlength) )
            elif boxtype == 'rectangular':
                leapin += ('solvateBox s %s %s 0.75\n' %
                           (self.solvent_box, boxlength) )
            elif boxtype == 'set':      # assumes coords already centered
                if boxdata:
                    leapin += ('set default nocenter on\nsetBox s centers %s\n'
                               % boxdata)
                else:
                    leapin += 'setBox s vdw\n'
            else:
                raise errors.SetupError('Unknown box type: %s' % boxtype)
        else:
            self.amber_top = const.LEAP_VACUUM + self.TOP_EXT
            self.amber_crd = const.LEAP_VACUUM + self.RST_EXT
            self.amber_pdb = const.LEAP_VACUUM + const.PDB_EXT

        if self._parm_overwrite:        # how cute...
            self.amber_top = self._parm_overwrite + self.TOP_EXT
            self.amber_crd = self._parm_overwrite + self.RST_EXT
            self.amber_pdb = self._parm_overwrite + const.PDB_EXT

        if neutralize:
            nions = abs(round(self.charge) )

            # add explicit number of ions because leap just truncates the int!
            if self.charge < 0.0:
                leapin += 'addions s Na+ %i\n' % nions
            elif self.charge > 0.0:
                leapin += 'addions s Cl- %i\n' % nions

        leapin += ('saveAmberParm s "%s" "%s"\nsavepdb s "%s"\n' %
                   (self.amber_top, self.amber_crd, self.amber_pdb) )

        if remove_first:
            # NOTE: this only removes the first unit!
            leapin += ('remove s s.1\nsaveAmberParm s "%s" "%s"\n'
                       'savepdb s "%s"\n' % (const.NOT_FIRST_TOP,
                                             const.NOT_FIRST_CRD,
                                             const.NOT_FIRST_PDB) )

        leapin += 'quit\n'

        self.sander_crd = self.amber_crd

        return leapin


    def setup_MDEngine(self, mdprog = 'sander', mdpref = '', mdpost = ''):
        """
        Instantiate MD engine.
        """

        self.get_box_dims(self.sander_crd)

        self.mdengine = self.MDEngine(self.amber_top, self.amber_crd,
                                      self.sander_crd, self.sander_rst,
                                      self.amber_pdb, self.box_dims,
                                      self.solvent, mdprog, mdpref, mdpost)


    @report
    def minimize(self, namelist = '%ALL', nsteps = 100, ncyc = 10,
                 restraint = '', restr_force = 10.0):
        """
        Use the AMBER/sander module to minimize a system.

        :param config: MD paramers, if start with '%' respective default is
           chosen, otherwise a full sander namelist as string
        :type config: string
        :param nsteps: maximum number of conjugated gradient steps
        :type nsteps: integer
        :param ncyc: number of initial steepest decent steps
        :type ncyc: float
        :param restraint: pre-defined restraints
        :type restraint: string
        :param restr_force: force constant for the restraints in kcal/mol/AA
        :type restr_force: float
        :raises: SetupError
        """

        if not self.amber_top or not self.amber_crd:
            raise errors.SetupError('Topology and/or coordinates missing.  '
                                    'Please run amber_create_top first.')

        self.mdengine.minimize(namelist, nsteps, ncyc, restraint, restr_force)


    @report
    def md(self, namelist = '', nsteps = 1000, T = 300.0, p = 1.0,
           restraint = '', restr_force = 10.0, nrestr = 1, wrap = True):
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
        :param restraint: pre-defined restraints
            backbone, heavy, protein, notligand, notsolvent or empty string
            for no restraints.  Any other string is passed on as a mask.  Must
            be set to any true value if namelist is not a preset template and
            restraints are to be used in the custom namelist.
        :type restraint: string
        :param restr_force: force constant for the restraints in kcal/mol/AA
        :type restr_force: float
        :param nrel: number of restraint relaxation steps
        :type nrel: integer
        :param wrap: wrap coordinates in a periodic system
        :type wrap: bool
        :raises: SetupError
        """

        if not namelist:
            raise errors.SetupError('No namelist supplied.')

        # FIXME: clean up nrestr to be consistent between MD programs
        self.mdengine.md(namelist, nsteps, T, p, restraint, restr_force,
                         nrestr, wrap)


    def to_rst7(self):
        self.mdengine.to_rst7()


    @report
    def flatten_rings(self):
        """
        Make aromatic rings planar.  Useful for MC like Sire.

        :raises: Sire error
        """


        import Sire.IO
        import Sire.Mol

        import warnings

        # FIXME: suppress a warning over a fmcs/Sire double data type
        #        registration collision
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', 'to-Python converter '
                                    'for.*already registered')
            import Sire.Move

        import Sire.Units
        import Sire.Config              # NOTE: make sure paths are set
                                        #       correctly in __init__.py!


        logger.write('flattening rings of residues %s' %
                     ', '.join(const.AROMATICS) )

        amber = Sire.IO.Amber()
        molecules = amber.readCrdTop(self.amber_crd, self.amber_top)[0]

        zmat_maker = Sire.IO.ZmatrixMaker()
        protein_zmatrices = os.path.join(Sire.Config.parameter_directory,
                                         'amber.zmatrices')
        zmat_maker.loadTemplates(protein_zmatrices)

        newmols = Sire.Mol.Molecules()

        molnums = molecules.molNums()

        for idx in range(0, len(molnums)):
            molnum = molnums[idx]
            curr_mol = molecules.at(molnum).molecule()

            try:
                curr_mol = zmat_maker.applyTemplates(curr_mol)
            except UserWarning as error:
                error_type = re.search(const.RE_SIRE_ERROR_STR,
                                       str(error)).group(1)

                if (error_type == 'SireError::invalid_key' or 
                    error_type == 'SireBase::missing_property'):
                    newmols.add(curr_mol)
                    continue
                else:
                    raise error

            zmatrix = curr_mol.property('z-matrix')
            zmatrixcoords = Sire.Move.ZMatrixCoords(zmatrix, curr_mol)
            zmatlines = zmatrixcoords.lines()

            for zline in zmatlines:
                at0 = curr_mol.select(zline.atom() )
                at1 = curr_mol.select(zline.bond() )
                at2 = curr_mol.select(zline.angle() )
                at3 = curr_mol.select(zline.dihedral() )

                if (str(at0.residue().name().value()) in const.AROMATICS and 
                    str(at1.residue().name().value()) in const.AROMATICS and
                    str(at2.residue().name().value()) in const.AROMATICS and 
                    str(at3.residue().name().value()) in const.AROMATICS):

                    dihname = '%s-%s-%s-%s' % (at0.name().value(),
                                               at1.name().value(), 
                                               at2.name().value(),
                                               at3.name().value() )

                    if dihname in const.TORSIONS_TO_ZERO:
                        zmatrixcoords.setDihedral(zline.atom(),
                                                  0.0 * Sire.Units.degrees)
                    elif dihname in const.TORSIONS_TO_ONEHUNDREDEIGHTY:
                        zmatrixcoords.setDihedral(zline.atom(),
                                                  180.0 * Sire.Units.degrees)

            ncoor = zmatrixcoords.toCartesian()
            molec = curr_mol.edit().setProperty('coordinates', ncoor).commit()
            newmols.add(molec)

        Sire.IO.PDB().write(newmols, const.FLAT_RINGS_FILE)

        self.amber_pdb = const.FLAT_RINGS_FILE
        self.mol_file = self.amber_pdb
        self.mol_fmt = 'pdb'


    def get_box_dims(self, rst_filen = None):
        """
        Get box dimensions from sander .rst file.

        :param rst_filen: name of .rst file
        :type rst_filen: string
        """

        if not rst_filen:
            rst_filen = self.sander_rst

        with open(rst_filen, 'r') as rst:
            for line in rst:
                self.box_dims = line

        with open(const.BOX_DIMS, 'w') as box:
            box.write(self.box_dims)

        self.box_dims = self.box_dims.split()


    def get_box_info(self):
        """
        Use Sire to get information about the system.

        :returns: volume
        :returns: density
        :returns: box dimensions
        :rtype: float
        """

        # Sire.Mol.Molecules, Sire.Vol.PeriodicBox or Sire.Vol.Cartesian
        molecules, space = \
                   Sire.IO.Amber().readCrdTop(self.amber_crd, self.amber_top)

        volume, density = 0.0, 0.0
        x, y, z = 0.0, 0.0, 0.0
        
        if space.isPeriodic():
            volume = space.volume().value()

            # NOTE: currently rectangular box only
            x = space.dimensions().x()
            y = space.dimensions().y()
            z = space.dimensions().z()
 
            total_mass = 0.0
            
            for num in molecules.molNums():
                mol = molecules.at(num).molecule()

                for atom in mol.atoms():
                    total_mass += atom.property('mass').value()

            density = total_mass * const.AMU2GRAMS / volume

        return volume, (x, y, z), density
