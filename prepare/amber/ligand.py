#  Copyright (C) 2012-2014  Hannes H Loeffler
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
A class to build a ligand with FESetup.  Derives from Common.

The Ligand class does parameterization and coordinate manipulation
based on a single coordinate input file.  AMBER topology and coordinate
files can be created both for vacuum and solution.  Flexibility
information for Sire can be computed too.  A conformationl search tool
(Openbabel) may be used to determine unbound conformations.

Ligand methods:
    preminimize, amber_param, conf_search, flex, amber_create_top,
    transfer_charges
"""


__revision__ = "$Id$"



import os, math, shutil

import openbabel as ob

import FESetup
from FESetup import const, errors, logger
from . import dlfield
from FESetup.modelconf import ModelConfig
from common import *
import utils

import Sire.IO



SQM_OUT = 'sqm.out'

GAUSS_INP = 'esp.in'
GUS_INP = 'esp.inp'
GUS_HEADER = '''\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=OPTIMIZE COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. $END
 $STATPT NSTEP=100 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY $END
 $GUESS  GUESS=HUCKEL $END
'''


class Ligand(Common):
    """The ligand setup class."""

    # import force field independent functionality (to avoid mixins)
    from FESetup.prepare.ligutil import prepare, align, conf_search, flex, \
         preminimize


    def __init__(self, ligand_name, basedir, start_file = 'ligand.pdb',
                 start_fmt = 'pdb', workdir = const.LIGAND_WORKDIR,
                 frcmod = const.LIGAND_FRCMOD_FILE, overwrite = False):
        """
        :param ligand_name: name of the ligand, will be used as directory name
        :type ligand_name: string
        :param basedir: base directory containing start_file
        :type basedir: string
        :param start_file: the file name of the ligand
        :type start_file: string
        :param start_fmt: format of the ligand file
        :type start_fmt: string
        :param workdir: output work directory
        :type workdir: string
        :param frcmod: name of the leap frcmod file
        :type frcmod: string
        :param overwrite: overwrite files in the working directory from basedir
        :type overwrite: bool
        """

        super(Ligand, self).__init__(ligand_name, basedir, workdir, overwrite)

        self.mol_file = start_file
        self.mol_fmt = start_fmt
        self.frcmod = frcmod       # force field modifications e.g. via parmchk

        # reference to original coordinates, may be converted to other format
        # in convert()
        self.orig_file = start_file
        self.orig_fmt = start_fmt

        self.charge = 0.0

        self.ref_file = ''
        self.ref_fmt = ''

        # FIXME: we assume that we are working on the bound ligand!
        self.mol_atomtype = 'sybyl'

        self.workdir = workdir

        self.model = ModelConfig(ligand_name)


    @report
    def param(self, sqm_strategy = None):
        """
        Compute symmetrized AM1/BCC charges and generate missing forcefield
        parameters. Runs antechamber, respgen, parmchk. Finally generated MOL2
        file is in Sybyl format.  GAFF atom names are needed internally by AMBER.

        :param sqm_strategy: a strategy patter using preminimize() and setting
           the SCF convergence criterion for sqm
        :type sqm_strategy: list of 2-tuples
        :raises: SetupError
        """

        logger.write('Deriving AMBER/GAFF force field parameters')

        antechamber = utils.check_amber('antechamber')

        ac_cmd = [
            '-i %s' % self.mol_file,    # input file
            '-fi %s' % self.mol_fmt,    # input file format
            '-o %s' % const.LIGAND_AC_FILE, # output file
            '-fo ac',                   # output file format
            '-c bcc',                   # charge method
            '-nc %s' % str(self.charge), # net molecular charge
            '-m 1',                     # FIXME: spin multiplicity (sqm only 1)
            '-df 2',                    # 0 = mopac, 2 = sqm
            '-at gaff',                 # write GAFF types
            '-du y',                    # fix duplicate atom names
            '-an y',                    # adjust atom names
            '-j 4',                     # atom/bond type prediction = full
            '-s 2',                     # status information = verbose
            '-eq 2',                    # equalise atom charges (path+geometry)
            '-pf y',                    # clean up temporary files
            '-rn %s' % const.LIGAND_NAME
            ]

        # NOTE: The main problem is SCF convergence. If this happens MM
        #       minimisation is used to hope to obtain a better structure with a
        #       better wavefunction.  This obviously depends on a sensible
        #       assignment of force field parameters which may fail if the
        #       structure is "too" distorted and no bonding information, etc. are
        #       available a priori.
        if not sqm_strategy:
            sqm_strategy = (
                (0, '1.0d-10', 1, 1, 0, 1000),
                (50, '1.0d-10', 1, 1, 0, 1000),
                (0, '1.0d-9', 1, 1, 0, 1000),
                (50, '1.0d-9', 1, 1, 0, 1000),
                (50, '1.0d-9', 0, 1, 0, 1000)
                # If the molecule has not converged yet there may be some deeper
                # problem.  Thus no attempt at a more elaborate protocol, e.g.
                #   scfconv=1.0d-9,vshift=0.1,ndiis_attempts=200,ndiis_matrices=6
                # for ZINC03814826/28/31/38; ZINC03814832 can't be converged
                # Some of those "break": proton transfer, "decarboxylation"
                )

        logger.write('Optimizing structure and creating AM1/BCC charges')

        for premin, scfconv, tight, psdiag, ndiisa, maxcyc in sqm_strategy:
            converged = False

            if premin:
                self.preminimize(nsteps = premin)
            
            sqm_nlv = ("qm_theory='AM1',grms_tol=0.0002,tight_p_conv=%i,\n  "
                       "scfconv=%s,itrmax=500,pseudo_diag=%i,\n  "
                       "ndiis_attempts=%i,maxcyc=%i,\n" %
                       (tight, scfconv, psdiag, ndiisa, maxcyc) )
            ek = ['-ek "%s"' % sqm_nlv]  # sqm namelist variables

            # FIXME: Buffering messes with the stdout output order of
            #        antechamber (last line comes first).  Use stdbuf, pexpect
            #        or pty (probably Linux only)?
            err = utils.run_amber(antechamber, ' '.join(ac_cmd + ek) )

            if err:
                if 'the assigned bond types may be wrong' in err[0]:
                    logger.write('Error: antechamber failed to assign '
                                 'atom/bond types properly\n')
                    raise errors.SetupError('antechamber cannot assign atom '
                                            'and/or bond types, check input '
                                            'structure, e.g. with acdoctor')

                sce = False

                with open(SQM_OUT, 'r') as sqm:
                    for line in sqm:
                        if 'Unable to achieve self consistency' in line:
                            logger.write('Warning: SCF has not converged '
                                         'with %i %s\n' % (premin, scfconv) )
                            sce = True
                            break

                        if 'odd number of electrons' in line:
                            logger.write('Error: odd electron number\n')
                            raise errors.SetupError('wrong ligand charge, or '
                                                    'radical')

                if not sce:
                    raise errors.SetupError('unknown error see log file '
                                            'and %s file' %
                                            os.path.join(self.dst, SQM_OUT) )
            else:
                converged = True
                break

        if not converged:
            if sce:
                logger.write('Error: SCF has not converged\n')
                raise errors.SetupError('SCF has not converged')
            else:
                logger.write('Error: failed to produce atom charges\n')
                raise errors.SetupError('failed to produce atom charges')

        logger.write('SCF has converged with %i preminimisation steps and '
                     'scfconv = %s kcal/mol\n' % (premin, scfconv) )

        # antechamber does not write the optimised  coordinates from sqm.pdb
        # into const.LIGAND_AC_FILE
        # FIXME: do not read and write to same file?
        utils.run_amber(antechamber,
                        '-i %s -fi ac '
                        '-o %s -fo ac '
                        '-a %s -fa pdb -ao crd '
                        '-s 2 -pf y' %
                        (const.LIGAND_AC_FILE, const.LIGAND_AC_FILE,
                         const.SQM_PDB_FILE) )

        ngconv = 0
        H_form = 'unknown'
        grad = 'unknown'

        with open(SQM_OUT, 'r') as sqm:
            for line in sqm:
                if line.startswith('xmin'):
                    ngconv = int(line[4:10].strip() )
                    H_form = line[10:33].strip()
                    grad = line[33:].strip()

        if ngconv >= maxcyc:
            logger.write('Warning: maximum number of geometry optimisation '
                         'steps reached (%i), gradient = %s '
                         '(grms_tol=0.0002), check %s file\n'
                         % (maxcyc, grad, SQM_OUT) )
        else:
            logger.write('Geometry has converged after %i steps, heat of '
                         'formation: %s and gradient = %s\n' %
                         (ngconv, H_form, grad) )


        self._parmchk(const.LIGAND_AC_FILE, 'ac', self.frcmod)

        charges = []

        with open(const.LIGAND_AC_FILE, 'r') as acfile:
            for line in acfile:
                if line[:4] == 'ATOM':
                    charges.append(float(line[54:64]) )

        if filter(lambda ch: math.fabs(ch) > const.MAX_CHARGE, charges):
            logger.write('Warning: some atom charges > %.2f' %
                         const.MAX_CHARGE)

        total_charge = sum(charges)
        dec_frac = total_charge - round(total_charge)

        if abs(dec_frac) > const.MAX_CHARGE_DIFF:
            logger.write('Warning: total molecule charge (%f) is far from '
                         'being an integer' % total_charge)
        
        corr = dec_frac / len(charges)

        for idx, charge in enumerate(charges):
            charges[idx] = charge - corr

        with open(const.CORR_CH_FILE, 'w') as chfile:
           for charge in charges:
               chfile.write('%f\n' % charge)

        utils.run_amber(antechamber,
                        '-i %s -fi ac '
                        '-o %s -fo ac '
                        '-cf %s -c rc'
                        ' -s 2 -pf y' %
                        (const.LIGAND_AC_FILE, const.CORR_AC_FILE,
                         const.CORR_CH_FILE) )

        shutil.copyfile(const.LIGAND_AC_FILE,
                        const.LIGAND_AC_FILE + os.extsep + '0')
        shutil.move(const.CORR_AC_FILE, const.LIGAND_AC_FILE)

        self.charge = sum(charges)
        logger.write('Total molecule charge is %.2f\n' % self.charge)

        with open(const.CHARGE_FILE, 'w') as chf:
            chf.write('%s' % self.charge)

        self.ref_file = self.mol_file
        self.ref_fmt = self.mol_fmt

        #self.model.add_file(self.mol_file)
        #self.model.add_file(const.LIGAND_AC_FILE)
        #self.model['charge'] = self.charge
        #self.model['charge.filename'] = const.LIGAND_AC_FILE
        #self.model['charge.filetype'] = 'ac'
        #self.model['charge.type'] = 'AM1-BCC'

        #self.model['forcefield'] = 'gaff'
        #self.model['type'] = 'ligand'
        #self.model['isvalid'] = False
        #self.model['ismorph'] = False
        #self.model['supports_md'] = True
        #self.model['supports_mc'] = False


    def _parmchk(self, infile, informat, outfile):
        if self.parmchk_version > 1:
            parmchk = utils.check_amber('parmchk%s' %
                                        str(self.parmchk_version) )
        else:
            parmchk = utils.check_amber('parmchk')

        logger.write('Creating frcmod file')

        params = '-i %s -f %s -o %s -a N ' % (infile, informat, outfile)

        # FIXME: parmchk only reads on parmfile, possible solution: write
        #        temporary frcmod files and paste together?
        if self.ff_addons:
            addon = self.ff_addons[0]

            # FIXME: Can it get any uglier? Consistent file naming, ey...
            if addon.startswith('GLYCAM_06'):
                addon = addon[:10]
            
            params += ' -p %s' % (os.path.join(os.environ['AMBERHOME'], 'dat',
                                              'leap', 'parm', addon) +
                                  os.extsep + 'dat')

        utils.run_amber(parmchk, params)


    @report
    def create_top(self, boxtype = '', boxlength = '10.0', boxfile = None,
                   align = False, neutralize = False, addcmd = '',
                   addcmd2 = ''):
        """
        Generate an AMBER topology file via leap. Leap requires atom names in
        GAFF format to match against GAFF force field database.  Finally
        generated MOL2 file is in GAFF format.

        :param boxtype: rectangular, octahedron or set (set dimensions explicitly)
        :param boxlength: side length of the box
        :param boxfile: name of file containing box dimensions
        :param align: align solute along the principal axes
        :param neutralize: neutralise the system
        :param addcmd: inject additional leap commands
        :type boxtype: string
        :type boxlength: float
        :type boxfile: string
        :type align: bool
        :type neutralize: bool
        :type addcmd: string
        """

        # we allow the user to have their own leap input file which is used
        # instead of the autogenerated one
        if os.access(const.LEAP_IN, os.F_OK):
            self.amber_top = const.LEAP_IN + self.TOP_EXT
            self.amber_crd = const.LEAP_IN + self.RST_EXT
            self.amber_pdb = const.LEAP_IN + const.PDB_EXT

            utils.run_leap(self.amber_top, self.amber_crd, program = 'tleap',
                           script = const.LEAP_IN)

            return

        if self.mol_fmt == 'pdb':
            load_cmd = 'loadpdb "%s"' % self.mol_file
        elif self.mol_fmt == 'mol2':
            if self.mol_atomtype != 'gaff':
                mol_file = const.GAFF_MOL2_FILE
                antechamber = utils.check_amber('antechamber')

                utils.run_amber(antechamber,
                                '-i %s -fi ac '
                                '-o %s -fo mol2 '
                                '-at gaff -s 2 -pf y' %
                                (const.LIGAND_AC_FILE, mol_file) )
            else:
                mol_file = self.mol_file

            load_cmd = 'loadmol2 "%s"' % mol_file
        else:
            raise errors.SetupError('Leap unsupported input format: %s (only '
                                    'mol2 and pdb)' % self.mol_fmt)

        leapin = '''%s
source leaprc.gaff
%s
%s
mods = loadAmberParams "%s"
s = %s
savemol2 s leap.mol2 1
%s\n''' % (self.ff_cmd, self.solvent_load, addcmd, self.frcmod, load_cmd,
           addcmd2)

        leapin += self._amber_top_common(boxtype, boxlength, boxfile, align,
                                         neutralize)


        # Strangely, sleap does not create sander compatible top files with
        # TIP4P but tleap does.  Sleap also crashes when @<TRIPOS>SUBSTRUCTURE
        # is missing.  Sleap has apparently been abandonded.
        utils.run_leap(self.amber_top, self.amber_crd, 'tleap', leapin)

        # create DL_FIELD UDFF/PDB for vacuum case
        if not boxtype:
            amber = Sire.IO.Amber()

            try:
                mols = amber.readCrdTop(self.amber_crd, self.amber_top)[0]
            except UserWarning as error:
                raise errors.SetupError('error opening %s/%s: %s' %
                                        (self.amber_crd, self.amber_top, error) )

            # there should be only one molecule
            nmols = mols.nMolecules()
            if nmols > 1:
                return                  # FIXME: don't write this when pert top
                raise errors.SetupError('BUG: only one molecule expected, '
                                        'found %i' % nmols)

            lig = mols.molNums()[0]
            dlfield.dlf_write(mols.at(lig).molecule(), '_AG')

        self.model.add_file(self.mol_file)
        self.model.add_file(self.amber_top)
        self.model.add_file(self.amber_crd)

        self.model['crd.filename'] = self.amber_crd
        self.model['crd.filetype'] = 'amber'
        self.model['top.filename'] = self.amber_top
        self.model['top.filetype'] = 'amber'


    @report
    def mk_esp(self, program = 'gauss', gkeys = '', gmem = '', gnproc = '',
               gus_header = GUS_HEADER):
        """
        Create an input file for ESP calculation using the MK scheme.
        The input is written for either Gaussian or Gamess-US.

        :param program: ab initio QM program, either gauss or gus
        :type program: string
        :param gkeys: addition keys for Gaussian (antechamber)
        :type gkeys: string
        :param gmem: memory information for Gaussian (antechamber)
        :type gmem: string
        :param gnproc: number of processors for Gaussian (antechamber)
        :type gnproc: string
        :param gus_header: Gamess-US ESP control parameters
        :type gus_header: string
        """

        if program == 'gauss':
            antechamber = utils.check_amber('antechamber')

            ac_cmd = [
                '-i %s' % self.mol_file,
                '-fi %s' % self.mol_fmt,
                '-o %s' % GAUSS_INP,
                '-fo gcrt -pf y'
                ]

            if gkeys:
                ac_cmd.append('-gk "%s"' % gkeys)
            if gmem:
                ac_cmd.append('-gm "%s"' % gmem)
            if gnproc:
                ac_cmd.append('-gn "%s"' % gnproc)

            utils.run_amber(antechamber, ' '.join(ac_cmd) )
        elif program == 'gus':
            conv = ob.OBConversion()
            conv.SetInAndOutFormats(self.mol_fmt, 'gamin')

            obm = ob.OBMol()
            conv.ReadFile(obm, self.mol_file)

            inp = conv.WriteString(obm)
            nl = inp.find('\n') + 1     # skip first line
            
            with open(GUS_INP, 'w') as gus:
                gus.writelines (gus_header + inp[nl:])
                
        else:
            raise errors.SetupError('Unknow QM program %s' % program)


    def get_charge(self):
        """
        Read the ligand charge from the charge file.
        """

        charge_file = os.path.join(self.dst, const.CHARGE_FILE)

        with open(charge_file, 'r') as chf:
            self.charge = float(chf.read() )


    def get_topcrd(self):
        """
        :returns: file names of current topology and rst files
        """
        return self.amber_top, self.amber_crd


    def set_atomtype(self, atomtype):
        """
        :param atomtype: set the current force field atom type to either 'gaff'
           or 'sybyl'
        :type atomtype: string
        """

        if atomtype not in const.KNOWN_MOL2_ATOMTYPES:
            raise errors.SetupError('Only %s atom types are supported' %
                                    str(const.KNOWN_MOL2_ATOMTYPES) )

        self.mol_atomtype = atomtype


    def new_model(self, name = 'unnamed'):
        """
        Create a new model.

        :param name: model name
        :type name: string
        """

        
        self.model['name'] = name
        self.model.filename = name + const.MODEL_EXT

        self.model.write()

        shutil.move(self.model.filename, os.path.join(self.topdir,
                                                      self.model.filename) )
