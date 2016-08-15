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



def ssbonds(ss_file, offset=0):
    """
    Read file with SS-bond information.

    :param ss_file: filename with disulfide bond information
    :type ss_file: str
    :param offset: offset to work around leap indexing
    :type offset: int

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
                a, b = int(a) + offset, int(b) + offset

                if a == b:
                    raise errors.SetupError('malformed file %s in '
                            'line %i: residue numbers cannot be the '
                            'same' % (lcnt, ss_file) )

                pairs.append( (a, b) )
            except (ValueError, TypeError) as why:
                raise errors.SetupError('malformed file %s in '
                            'line %i: %s' % (ss_file, lcnt, why) )

    return pairs


class Common(object):
    """
    Common methods for the setup protocol.  Used through subclasses for
    ligand, protein an complex setup.  Methods are expected to be called
    in the right 'order' thus a second layer needs to take control of this.
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
    SSBONDS_OFFSET = 0


    def __new__(cls, *args, **kwargs):
        """Prevent direct instantiation."""

        if cls is Common:
            raise TypeError('Common may not be instantiated directly.')

        return object.__new__(cls)


    def __init__(self, mol_name):
        """
        :param mol_name: (file) name of molecules
        :type mol_name: string
        :raises: SetupError
        """

        self.mol_name = mol_name
        self.mol_file = ''
        self.mol_fmt = ''

        self.sander_rst = ''
        self.sander_crd = ''
        self.ssbond_file = ''
        self.namd_prefix = ''

        self.min_no = 0                 # number of current minimisation step
        self.md_no = 0                  # number of current MD step

        self.charge = 0.0
        self.box_dims = []
        self.volume = 0.0
        self.density = 0.0

        self.amber_top = ''             # only set in _amber_top_common
        self.amber_crd = ''             # only set in _amber_top_common
        self.amber_pdb = ''

        self._parm_overwrite = None
        self.mdengine = None

        self.leap = Leap(self.force_fields, self.solvent_load)

        # set in __init__.py before this __init__(): do not set here or anywhere
        # else as they are meant to be constants for the whole class hierarchy!
        #self.ff_cmd, self.ff_addon, self.solvent, self.solvent_load 
        #self.solvent_box, self.MDEngine, self.parmchk_version, self.gaff


    # FIXME: remove
    @staticmethod
    def copy_files(srcs, filenames=None, overwrite=False):
        """
        Copy files from basedir to the current directory.

        :param src: the source directories to copy from
        :type src: list of str
        :param filenames: the filenames to be copied
        :type filenames: list of str
        :param overwrite: overwrite the files in the current dir?
        :type overwrite: bool
        """

        try:
            for src in srcs:
                logger.write('%sopying directory contents of %s to %s' %
                             ('Overwrite mode: c' if overwrite else 'C', src,
                              os.getcwd()))

                if not filenames:
                    filenames = os.listdir(src)

                for filename in filenames:
                    if overwrite or not os.access(filename, os.F_OK):
                        src_file = os.path.join(src, filename)

                        # FIXME: only here to accommodate Complex
                        if os.access(src_file, os.F_OK):
                            shutil.copy(src_file, '.')

        except OSError as why:
            raise errors.SetupError(why)


    def _amber_top_common(self, boxtype='', boxlength='10.0', neutralize=0,
                          align=None, remove_first=False, conc=0.0, dens=1.0):
        """Common scripting commands for leap.  Internal function only.

        :param boxtype: rectangular, octahedron or set
        :type boxtype: string
        :param boxlength: box length in Angstrom
        :type boxlength: float
        :param neutralize: 1=minimum number of ions, 2=use conc
        :type neutralize: int
        :param align: align axes?
        :type align: bool
        :param remove_first: remove first molecule?
        :type remove_first: bool
        :param conc: ion concentration in mol/litres
        :type conc: float
        :param dens: expected target density
        :type dens: float
        :raises: SetupError
        """

        leapin = self.leap.generate_init()

        if os.access(const.SSBOND_FILE, os.R_OK):
            pairs = ssbonds(const.SSBOND_FILE, self.__class__.SSBONDS_OFFSET)
            cmd = []

            for a, b in pairs:
                cmd.append('bond s.%i.SG s.%i.SG\n' % (a, b) )

            leapin += ''.join(cmd)

            self.ssbond_file = const.SSBOND_FILE

        boxdata = None

        if align:
            leapin += 'alignAxes s\n'

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
                leapin += ('set default nocenter on\nsetBox s centers {%s}\n'
                           % ' '.join([str(float(b) / 2.0)
                                       for b in self.box_dims[:3]]))
            elif  boxtype == 'setvdw':
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

        if neutralize == 1:             # minimum number of ions
            nions = abs(round(self.charge) )

            # add explicit number of ions because leap just truncates the int!
            if self.charge < 0.0:
                leapin += 'addIons s Na+ %i\n' % nions
            elif self.charge > 0.0:
                leapin += 'addIons s Cl- %i\n' % nions
        elif neutralize == 2:           # neutralise to set concentration
            self.get_box_info()

            self.amber_top = const.LEAP_IONIZED + self.TOP_EXT
            self.amber_crd = const.LEAP_IONIZED + self.RST_EXT
            self.amber_pdb = const.LEAP_IONIZED + const.PDB_EXT

            nions = abs(round(self.charge) )

            # FIXME: check if this is correct
            volume = self.volume * self.density / dens

            # 1 mol/l = 6.022140857*10^23 particles/litre (NIST)
            # 1 A^3   = 10^-27 l
            npart = round(0.0006022141 * conc * volume)

            if self.charge < 0.0:
                npos = npart + nions
                nneg = npart
            elif self.charge > 0.0:
                npos = npart
                nneg = npart + nions
            else:
                npos = nneg = npart

            logger.write('box info: V = %f A^3, rho = %f g/cc\n'
                         'charge = %f; computed #pos = %i, #neg = %i\n' %
                         (self.volume , self.density, self.charge, npos,
                          nneg) )

            leapin += 'addIonsRand s Na+ %i Cl- %i 2.0\n' % (npos, nneg)


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
        self.amber_crd = self.mdengine.sander_crd  # FIXME: only for AMBER


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

        # FIXME: do we also want to density?
        self.box_dims = self.mdengine.get_box_dims()
        self.amber_crd = self.mdengine.sander_crd  # FIXME: only for AMBER


    def to_rst7(self):
        """
        Ask MD engine to convert coordinates to rst7 format.

        FIXME: must be called explicitly!
        """

        self.mdengine.to_rst7()
        self.amber_crd = self.mdengine.sander_crd


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


    # called in morph.py (1x), common.py/setup_MDEngine (1x)
    def get_box_dims(self, filename=None):
        """
        Get box dimensions from sander .rst file.

        :param filename: name of .rst file
        :type filename: string
        """

        if not filename:
            filename = self.sander_rst

        with open(filename, 'r') as rst:
            for line in rst:
                self.box_dims = line

        self.box_dims = self.box_dims.split()


    # called in common.py/_amber_top_common (1x)
    def get_box_info(self):
        """
        Use Sire to get information about the system: volume, density,
        box dimensions.

        NOTE: Sire prmtop reader is slow!
        """

        # Sire.Mol.Molecules, Sire.Vol.PeriodicBox or Sire.Vol.Cartesian
        molecules, space = \
                   Sire.IO.Amber().readCrdTop(self.amber_crd, self.amber_top)

        if space.isPeriodic():
            self.volume = space.volume().value()  # in A^3

            # NOTE: currently rectangular box only
            x = space.dimensions().x()
            y = space.dimensions().y()
            z = space.dimensions().z()
            self.box_dims = (x, y, z)   # in Angstrom
 
            total_mass = 0.0            # in amu
            
            for num in molecules.molNums():
                mol = molecules.at(num).molecule()

                for atom in mol.atoms():
                    total_mass += atom.property('mass').value()

            # in g/cc
            self.density = total_mass * const.AMU2GRAMS / self.volume
