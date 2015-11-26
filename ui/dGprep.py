#  Copyright (C) 2013-2015  Hannes H Loeffler
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
Simple user interface for free energy simulation preparation.
"""

from __future__ import print_function
import os, sys

if ('%x' % sys.hexversion)[:3] != '207':
    print ('Python version 2.7 required, you have %i.%i' %
           sys.version_info[:2], file=sys.stderr)
    sys.exit(1)


__revision__ = "$Id$"
__version__ = '0.7.1'


vstring = 'FESetup SUI version: %s' % __version__
print ('\n=== %s ===\n' % vstring)


import os, shutil, argparse, signal, copy, glob, atexit
from collections import OrderedDict, namedtuple

import warnings                         

import FESetup.prepare as prep
from FESetup import const, errors, create_logger, logger
from FESetup.ui.iniparser import IniParser

# FIXME: That's here solely to suppress a warning over a fmcs/Sire double
# data type registration collision.  Impact limited as much as possible but
# still potentially dangerous.  Fix actual problem instead!
with warnings.catch_warnings():
    warnings.filterwarnings('ignore',
                            'to-Python converter for.*already registered')

    from FESetup import mutate



LIST_SEP = ','
MORPH_PAIR_SEP = '>'
COM_PAIR_SEP = const.PROT_LIG_SEP

SECT_DEF = 'globals'
SECT_LIG = 'ligand'
SECT_PROT = 'protein'
SECT_COM = 'complex'

ALL_SECTIONS = (SECT_DEF, SECT_LIG, SECT_PROT, SECT_COM)

DONE_LOG = '_mols.done'

PARAM_CMDS = {
    'frcmod': 'loadAmberParams',
    'prep': 'loadAmberPrep',
    'lib': 'loadOff'
    }


class dGprepError(Exception):
    pass


def check_dict(section, opts):
    """Check the sanity of options.  Signal unset values."""

    for key, val in opts[section].iteritems():
        if val == None:
            raise dGprepError('key [%s] %s is mandatory but is not set '
                              'in input file' % (section, key) )


def prelude(opts):
    lff = len(opts[SECT_DEF]['forcefield'])
    deff = defaults[SECT_DEF]['forcefield'][0]
    
    if lff < 1 or lff > 5:
        raise dGprepError('malformed "forcefield" key, must have 1-5 strings '
                          'required')

    if len(opts[SECT_DEF]['mdengine']) != 2:
        raise dGprepError('malformed "mdengine" key, 2 strings required')


    create_logger(opts[SECT_DEF]['logfile'])
    logger.write('\n%s\n\n' % vstring)
    atexit.register(lambda : logger.finalize() )

    ff_opts = opts[SECT_DEF]['forcefield']
    lu = len(ff_opts)
    ff_opts[lu:] = deff[lu:]

    ff_opts.append(opts[SECT_DEF]['ff_addons'])
    ff_opts.append(opts[SECT_DEF]['mdengine'][0])
    ff_opts.append(opts[SECT_DEF]['parmchk_version'])

    return prep.ForceField(*ff_opts)


def _param_glob(cmds):
    load_cmds = []
    
    for ext in cmds:
        files = glob.glob('*' + os.extsep + ext)

        for f in files:
            if f == const.LIGAND_FRCMOD_FILE:
                continue

            load_cmds.append(cmds[ext] + ' ' + os.path.abspath(f) + '\n')

    return ''.join(load_cmds)


def _minmd_done(dico):
    for key in dico:
        if '.nsteps' in key and dico[key]:
            return True

    return False


def make_ligand(name, ff, opts, short = False):
    """
    Prepare ligands for simulation: charge parameters, vacuum top/crd,
    confomer search + alignment (both optional), optionally hydrated
    top/crd and minimisation/MD simulation.
    """
    
    lig = opts[SECT_LIG]

    if not lig['basedir']:
        raise dGprepError('[%s] "basedir" must be set' % SECT_LIG)

    if not lig['file.format']:
        fmt = os.path.splitext(lig['file.name'])[1][1:]
    else:
        fmt = lig['file.format']


    with ff.Ligand(name, lig['basedir'], lig['file.name'], fmt,
                   overwrite = opts[SECT_DEF]['overwrite']) as ligand:

        if not os.access(lig['file.name'], os.F_OK):
            raise errors.SetupError('start file %s does not exist in %s' %
                                    (lig['file.name'], os.getcwd() ) )

        load_cmds = ''

        if opts[SECT_DEF]['user_params']:
            load_cmds = _param_glob(PARAM_CMDS)

        if short:
            return ligand, load_cmds

        if lig['skip_param']:
            if fmt != 'pdb' and fmt != 'mol2':
                raise dGprepError('When parameterisation is skipped, the input '
                                  'format must be PDB or MOL2')

            ligand.prepare('', lig['add_hydrogens'], lig['calc_charge'],
                           lig['correct_for_pH'], lig['pH'])
        elif not os.access(os.path.join(ligand.dst, const.GAFF_MOL2_FILE),
                           os.F_OK):
            # IMPORTANT: do not allow OpenBabel to add Hs, it may mess up
            # everything
            ligand.prepare('mol2', lig['add_hydrogens'], lig['calc_charge'],
                           lig['correct_for_pH'], lig['pH'])
            ligand.param(lig['gb_charges'])
        else: # FIXME: ugly
            ligand.prepare('', lig['add_hydrogens'], lig['calc_charge'],
                           lig['correct_for_pH'], lig['pH'])
            ligand.mol_file = const.GAFF_MOL2_FILE
            ligand.mol_fmt = 'mol2'

        ligand.create_top(boxtype = '', addcmd = load_cmds)

        if opts[SECT_DEF]['MC_prep']:
            ligand.flex()

        nconf = lig['conf_search.numconf']

        if nconf > 0:
            ligand.conf_search(numconf = nconf,
                               geomsteps = lig['conf_search.geomsteps'],
                               steep_steps = lig['conf_search.steep_steps'],
                               steep_econv = lig['conf_search.steep_econv'],
                               conj_steps = lig['conf_search.conj_steps'],
                               conj_econv = lig['conf_search.conj_econv'],
                               ffield = lig['conf_search.ffield'])
            ligand.align()

        # FIXME: also check for boxlength and neutralize
        if lig['box.type']:
            ligand.create_top(boxtype = lig['box.type'],
                              boxlength = lig['box.length'],
                              neutralize = lig['neutralize'],
                              addcmd = load_cmds, remove_first = False)

            restr_force = lig['min.restr_force']
            nsteps = lig['min.nsteps']

            ligand.setup_MDEngine(opts[SECT_DEF]['mdengine'][1],
                                  opts[SECT_DEF]['mdengine.prefix'],
                                  opts[SECT_DEF]['mdengine.postfix'])

            if nsteps > 0:
                do_min(ligand, lig)

            #ligand.md('%SHRINK', 200, 5.0, 1.0, ':LIG', 5.0, wrap = True)

            press_done = False
            nsteps = lig['md.heat.nsteps']

            if nsteps > 0:
                restr_force = lig['md.heat.restr_force']
                do_md(ligand, lig, 'heat')

            nsteps = lig['md.constT.nsteps']

            if nsteps > 0:
                restr_force = lig['md.constT.restr_force']
                do_md(ligand, lig, 'constT')

            nsteps = lig['md.press.nsteps']

            if nsteps > 0:
                restr_force = lig['md.press.restr_force']
                do_md(ligand, lig, 'press')
                press_done = True

            nrestr = lig['md.relax.nrestr']

            if nrestr > 0:
                # FIXME: unify with AMBER mdengine
                if opts[SECT_DEF]['mdengine'][0] == 'namd':
                    ligand.md('%RELRES', lig['md.relax.nsteps'],
                              lig['md.relax.T'], lig['md.relax.p'],
                              lig['md.relax.restraint'], restr_force,
                              nrestr, wrap = True)
                else:
                    sp = restr_force / (nrestr - 1)

                    for k in range(nrestr - 2, -1, -1):
                        if press_done:
                            nmlist = '%PRESS'
                        else:
                            nmlist = '%CONSTT'

                        ligand.md(nmlist, lig['md.relax.nsteps'],
                                  lig['md.relax.T'], lig['md.relax.p'],
                                  lig['md.relax.restraint'], sp * k,
                                  wrap = True)

            if opts[SECT_DEF]['mdengine'][0] != 'amber':
                if _minmd_done(lig):
                   ligand.to_rst7()

    return ligand, load_cmds


def make_protein(name, ff, opts, short = False):
    """
    Prepare proteins for simulation:
    """

    prot = opts[SECT_PROT]

    if not prot['basedir']:
        raise dGprepError('[%s] "basedir" must be set' % SECT_PROT)

    with ff.Protein(name, prot['basedir'], prot['file.name'],
                    overwrite = opts[SECT_DEF]['overwrite']) as protein:
        load_cmds = ''

        if opts[SECT_DEF]['user_params']:
            load_cmds = _param_glob(PARAM_CMDS)

        if short:
            return protein, load_cmds

        if prot['propka']:
            protein.protonate_propka(pH = prot['propka.pH'])

        protein.get_charge()    # must be done explicitly
        protein.create_top(boxtype = '')

        # FIXME: also check for boxlength and neutralize
        if prot['box.type']:
            protein.create_top(boxtype = prot['box.type'],
                               boxlength = prot['box.length'],
                               neutralize = prot['neutralize'],
                               align = prot['align_axes'],
                               addcmd = load_cmds, remove_first = True)

            restr_force = prot['min.restr_force']
            nsteps = prot['min.nsteps']

            protein.setup_MDEngine(opts[SECT_DEF]['mdengine'][1],
                                   opts[SECT_DEF]['mdengine.prefix'],
                                   opts[SECT_DEF]['mdengine.postfix'])

            if nsteps > 0:
                do_min(protein, prot)

            press_done = False

            #protein.md('%SHRINK', 200, 5.0, 1.0, ':LIG', 5.0, wrap = True)

            nsteps = prot['md.heat.nsteps']

            if nsteps > 0:
                restr_force = prot['md.heat.restr_force']
                do_md(protein, prot, 'heat')

            nsteps = prot['md.constT.nsteps']

            if nsteps > 0:
                restr_force = prot['md.constT.restr_force']
                do_md(protein, prot, 'constT')

            nsteps = prot['md.press.nsteps']

            if nsteps > 0:
                restr_force = prot['md.press.restr_force']
                do_md(protein, prot, 'press')
                press_done = True

            nrestr = prot['md.relax.nrestr']

            if nrestr > 0:
                if opts[SECT_DEF]['mdengine'][0] == 'namd':
                    protein.md('%RELRES', prot['md.relax.nsteps'],
                               prot['md.relax.T'], prot['md.relax.p'],
                               prot['md.relax.restraint'], restr_force,
                               nrestr, wrap = True)
                else:
                    sp = restr_force / (nrestr - 1)

                    for k in range(nrestr - 2, -1, -1):
                        if press_done:
                            nmlist = '%PRESS'
                        else:
                            nmlist = '%CONSTT'

                        protein.md(nmlist, prot['md.relax.nsteps'],
                                   prot['md.relax.T'],
                                   prot['md.relax.p'],
                                   prot['md.relax.restraint'], sp * k,
                                   wrap = True)

            if opts[SECT_DEF]['mdengine'][0] != 'amber':
                if _minmd_done(prot):
                   protein.to_rst7()

    return protein, load_cmds


def make_complex(prot, lig, ff, opts, load_cmds, short = False):

    com = opts[SECT_COM]

    with ff.Complex(prot, lig,
                    overwrite = opts[SECT_DEF]['overwrite']) as complex:

        if short:
            return complex, load_cmds

        complex.ligand_fmt = lig.mol_fmt
        complex.create_top(boxtype = '', addcmd = load_cmds)

        # FIXME: also check for boxlength and neutralize
        if com['box.type']:
            complex.create_top(boxtype = com['box.type'],
                               boxlength = com['box.length'],
                               neutralize = com['neutralize'],
                               align = com['align_axes'],
                               addcmd = load_cmds, remove_first = True)

            restr_force = com['min.restr_force']
            nsteps = com['min.nsteps']

            complex.setup_MDEngine(opts[SECT_DEF]['mdengine'][1],
                                   opts[SECT_DEF]['mdengine.prefix'],
                                   opts[SECT_DEF]['mdengine.postfix'])

            if nsteps > 0:
                do_min(complex, com)

            if opts[SECT_DEF]['MC_prep']:
                complex.prot_flex()
                complex.flatten_rings()

            press_done = False

            #complex.md('%SHRINK', 200, 5.0, 1.0, 'bb_lig', 5.0, wrap = True)

            nsteps = com['md.heat.nsteps']

            if nsteps > 0:
                restr_force = com['md.heat.restr_force']
                do_md(complex, com, 'heat')

            nsteps = com['md.constT.nsteps']

            if nsteps > 0:
                restr_force = com['md.constT.restr_force']
                do_md(complex, com, 'constT')

            nsteps = com['md.press.nsteps']

            if nsteps > 0:
                restr_force = com['md.press.restr_force']
                do_md(complex, com, 'press')
                press_done = True

            nrestr = com['md.relax.nrestr']

            if nrestr > 0:
                if opts[SECT_DEF]['mdengine'][0] == 'namd':
                    complex.md('%RELRES', com['md.relax.nsteps'],
                               com['md.relax.T'], com['md.relax.p'],
                               com['md.relax.restraint'], restr_force,
                               nrestr, wrap = True)
                else:
                    sp = restr_force / (nrestr - 1)

                    for k in range(nrestr - 2, -1, -1):
                        if press_done:
                            nmlist = '%PRESS'
                        else:
                            nmlist = '%CONSTT'

                        complex.md(nmlist, com['md.relax.nsteps'],
                                   com['md.relax.T'], com['md.relax.p'],
                                   com['md.relax.restraint'], sp * k,
                                   wrap = True)

            if opts[SECT_DEF]['mdengine'][0] != 'amber':
                if _minmd_done(com):
                    complex.to_rst7()

    return complex, load_cmds


def do_min(what, opts):
    #FIXME: unify
    if options[SECT_DEF]['mdengine'][0] == 'amber':
        prot = '%ALL'
    else:
        prot = '%STD'

    what.minimize(prot, opts['min.nsteps'], opts['min.ncyc'],
                  opts['min.restraint'], opts['min.restr_force'])

def do_md(what, opts, how):
    what.md('%%%s' % how.upper(), opts['md.%s.nsteps' % how],
            opts['md.%s.T' % how], opts['md.%s.p' % how],
            opts['md.%s.restraint' % how],
            opts['md.%s.restr_force' % how], wrap = True)



defaults = {}

# All valid keys with defaults
# None is used to signal values required to be set by the user
defaults[SECT_DEF] = {
    'logfile': ('dGprep.log', None),
    'forcefield': (['amber', 'ff14SB', 'tip3p', 'hfe'], ('list', LIST_SEP) ),
    'ff_addons': ([], ('list', LIST_SEP) ),
    'mdengine': (['amber', 'sander'], ('list', LIST_SEP) ),
    'mdengine.prefix': ('', None),
    'mdengine.postfix': ('', None),
    'parmchk_version': (2, (int, ) ),
    'FE_type': ('Sire', None),
    'softcore_type': ('', None),
    'remake': (False, ('bool', ) ),
    'mcs.timeout': (60, (int, ) ),      # int because of FMCS/C++
    'mcs.match_by': ('', None),
    'overwrite': (False, ('bool', ) ),
    'user_params': (False, ('bool', ) ),
    'MC_prep': (False, ('bool', ) ),
    }

defaults[SECT_LIG] = {
    'basedir': ('', None),
    'file.name': ('ligand.pdb', None),
    'file.format': ('', None),
    'molecules': ('', ('list', LIST_SEP) ),
    'morph_pairs': ('', ('pairlist', LIST_SEP, MORPH_PAIR_SEP) ),
    'box.type': ('', None),
    'box.length': (10.0, (float, ) ),
    'neutralize': (False, ('bool', ) ),
    'conf_search.numconf': (0, (int, ) ),
    'conf_search.geomsteps': (5, (int, ) ),
    'conf_search.steep_steps': (100, (int, ) ),
    'conf_search.steep_econv': (1.0E-4, (float, ) ),
    'conf_search.conj_steps': (250, (int, ) ),
    'conf_search.conj_econv': (1.0E-6, (float, ) ),
    'conf_search.ffield': ('mmff94', None),
    'calc_charge': (False, ('bool', ) ),
    'gb_charges': (False, ('bool', ) ),
    'add_hydrogens': (False, ('bool', ) ),
    'correct_for_pH': (False, ('bool', ) ),
    'pH': (7.4, (float, ) ),
    'skip_param': (False, ('bool', ) )
    }

defaults[SECT_PROT] = {
    'basedir': ('', None),
    'file.name': ('protein.pdb', None),
    'molecules': ('', ('list', LIST_SEP) ),
    'box.type': ('', None),
    'box.length': (10.0, (float, ) ),
    'neutralize': (False, ('bool', ) ),
    'align_axes': (False, ('bool', ) ),
    'propka': (False, ('bool', ) ),
    'propka.pH': (7.0, (float, ) ),
    }

defaults[SECT_COM] = {
    'pairs': ('', ('pairlist', LIST_SEP, COM_PAIR_SEP) ),
    'box.type': ('', None),
    'box.length': (10.0, (float, ) ),
    'neutralize': (False, ('bool', ) ),
    'align_axes': (False, ('bool', ) ),
    'min.nsteps': (100, (int, ) ),
    'min.ncyc': (10, (int, ) ),
    'min.restraint': ('notsolvent', None),
    'min.restr_force': (10.0, (float, ) ),
    }

_minmd = {
    'min.nsteps': (0, (int, ) ),
    'min.ncyc': (10, (int, ) ),
    'min.restraint': ('notsolvent', None),
    'min.restr_force': (10.0, (float, ) ),
    'md.heat.nsteps': (0, (int, ) ),
    'md.heat.T': (300.0, (float, ) ),
    'md.heat.p': (1.0, (float, ) ),
    'md.heat.restraint': ('notsolvent', None),
    'md.heat.restr_force': (10.0, (float, ) ),
    'md.constT.nsteps': (0, (int, ) ),
    'md.constT.T': (300.0, (float, ) ),
    'md.constT.p': (1.0, (float, ) ),
    'md.constT.restraint': ('notsolvent', None),
    'md.constT.restr_force': (10.0, (float, ) ),
    'md.press.nsteps': (0, (int, ) ),
    'md.press.T': (300.0, (float, ) ),
    'md.press.p': (1.0, (float, ) ),
    'md.press.restraint': ('notsolvent', None),
    'md.press.restr_force': (10.0, (float, ) ),
    'md.relax.nrestr': (0, (int, ) ),
    'md.relax.nsteps': (0, (int, ) ),
    'md.relax.T': (300.0, (float, ) ),
    'md.relax.p': (1.0, (float, ) ),
    'md.relax.restraint': ('notsolvent', None)
}

defaults[SECT_LIG].update(_minmd)
defaults[SECT_PROT].update(_minmd)
defaults[SECT_COM].update(_minmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs = '?',
                        help = 'input file in INI format, if not given then '
                        'just output defaults')
    parser.add_argument('--tracebacklimit', type = int, default = 0,
               help = 'set the Python traceback limit (for debugging)')
    args = parser.parse_args()

    options = IniParser(copy.deepcopy(defaults))

    if not args.infile:
        options.output()
        sys.exit(0)

    sys.tracebacklimit = args.tracebacklimit

    options.parse(args.infile, 'globals')

    for section in ALL_SECTIONS:
        check_dict(section, options)

    ff = prelude(options)

    logger.write('Force field and MD engine:\n%s\n' % ff)


    done = set()

    if not options[SECT_DEF]['remake']:
        try:
            with open(DONE_LOG, 'r') as done_log:
                for line in done_log:
                    done.add(line.rstrip() )
        except IOError:
            pass

        mode = 'a'
    else:
        mode = 'w'

    done_log = open(DONE_LOG, mode, 1)

    def sigusr2_handler(*args):
        done_log.flush()
        os.fsync(done_log.fileno() )

    signal.signal(signal.SIGUSR2, sigusr2_handler)


    # FIXME: We keep all molecule objects in memory.  For 2000 morph pairs
    #        this may mean more than 1 GB on a 64 bit machine.

    ### proteins

    proteins = {}
    prot_failed = []

    for prot_name in options[SECT_PROT]['molecules']:
        if prot_name not in done:
            print ('Making protein %s...' % prot_name)

            try:
                protein, cmds = make_protein(prot_name, ff, options)
                proteins[protein] = cmds
                done_log.write(prot_name + '\n')
            except errors.SetupError as why:
                prot_failed.append(prot_name)
                print ('ERROR: %s failed: %s' % (prot_name, why))
        else:
            protein, cmds = make_protein(prot_name, ff, options, True)
            proteins[protein] = cmds


    ### ligands

    ligands = {}
    Ligdata = namedtuple('Ligdata', ['ref', 'leapcmd'])
    lig_failed = []

    morph_pairs = options[SECT_LIG]['morph_pairs']

    if morph_pairs:
        mols = [val for pairs in morph_pairs for val in pairs]  # flatten
        uniq = list(OrderedDict( (val, None) for val in mols) )
        options[SECT_LIG]['molecules'] = uniq

    for lig_name in options[SECT_LIG]['molecules']:
        if lig_name not in done:
            print ('Making ligand %s...' % lig_name)

            try:
                ligand, cmds = make_ligand(lig_name, ff, options)
                ligands[lig_name] = Ligdata(ligand, cmds)
                done_log.write(lig_name + '\n')
            except errors.SetupError as why:
                lig_failed.append(lig_name)
                print ('ERROR: %s failed: %s' % (lig_name, why))
        else:
            if not options[SECT_LIG]['file.format']:
                fmt = os.path.splitext(options[SECT_LIG]['file.name'])[1][1:]
            else:
                fmt = options[SECT_LIG]['file.format']

            ligand, cmds = make_ligand(lig_name, ff, options, True)
            ligands[lig_name] = Ligdata(ligand, cmds)


    ### ligand morphs

    if morph_pairs:
        print('Morphs will be generated for %s' % options[SECT_DEF]['FE_type'])
        logger.write('Morphs will be generated for %s\n' %
                     options[SECT_DEF]['FE_type'])

    morphs = []
    morph_failed = []
    for pair in morph_pairs:
        try:
            l1 = ligands[pair[0] ]
            l2 = ligands[pair[1] ]
        except KeyError as why:
            name = pair[0] + const.MORPH_SEP + pair[1]
            morph_failed.append(name)
            print ('ERROR: %s failed: %s' % (name, why))
            continue

        ligand1 = l1.ref
        ligand2 = l2.ref

        cmd1 = l1.leapcmd
        cmd2 = l2.leapcmd

        with mutate.Morph(ligand1, ligand2, ff,
                          options[SECT_DEF]['FE_type'],
                          options[SECT_DEF]['softcore_type'],
                          options[SECT_DEF]['mcs.timeout'],
                          options[SECT_DEF]['mcs.match_by']) as morph:

            print ('Morphing %s to %s...' % pair)

            if (pair[1], pair[0]) in morph_pairs:
                rev = ligand2
            else:
                rev = None

            try:
                morph.setup(cmd1, cmd2)

                if options[SECT_LIG]['box.type']:
                    morph.create_coords(ligand1, cmd1, cmd2, rev)
            except errors.SetupError as why:
                morph_failed.append(morph.name)
                print ('ERROR: %s failed: %s' % (morph.name, why))

            morphs.append(morph)


    ### complexes

    complexes = {}
    com_failed = []

    # NOTE: does a complex for individual ligands need to be built when
    #       complex morphs are requested?
    for protein in proteins.keys():
        if proteins[protein] in prot_failed:
            print ('WARNING: not making complex with protein %s because '
                   'build failed' %  protein.mol_name)
            continue

        for lig_name in ligands:
            if lig_name in lig_failed:
                print ('WARNING: not making complex with ligand %s '
                       'because build failed' % lig_name)
                continue
                
            if options[SECT_COM]['pairs']:
                bfound = False

                for pair in options[SECT_COM]['pairs']:
                    if (pair[0] == protein.mol_name and \
                        pair[1] == lig_name) or \
                        (pair[0] == lig_name and \
                         pair[1] == protein.mol_name):
                        bfound = True
                        break
            else:
                bfound = True

            cmds = ligands[lig_name].leapcmd + proteins[protein]

            if bfound:
                name = '%s:%s' % (protein.mol_name, lig_name)

                if name not in done:
                    print ('Making complex from protein %s and ligand %s... ' %
                           (protein.mol_name, lig_name) )

                    try:
                        complex, cmds = make_complex(protein,
                                                     ligands[lig_name].ref, ff,
                                                     options, cmds)

                        # NOTE: no real need here to limit the list of complexes
                        #       as it is checked again below
                        for pair in morph_pairs:
                            if complex.ligand.mol_name == pair[0]:
                                complexes[complex] = cmds

                        done_log.write(name + '\n')
                    except errors.SetupError as why:
                        com_failed.append(name)
                        print ('ERROR: %s failed: %s' % (name, why))
                else:
                    complex, cmds = make_complex(protein.mol_name, lig_name, ff,
                                                 options, cmds, True)
                    complexes[complex] = cmds


    ### complex morphs

    for morph in morphs:
        for complex in complexes.keys():
            name = complex.mol_name + '/' + morph.name

            # FIXME: Complex has no ligand component after restart
            if complex.ligand.mol_name == morph.initial_name and \
                name not in done:

                print('Creating complex %s with ligand morph %s...' %
                      (complex.mol_name, morph.name) )

                try:
                    morph.create_coords(complex, complexes[complex], '')
                    #done_log.write(name + '\n')
                except errors.SetupError as why:
                    morph_failed.append(name)
                    print ('ERROR: complex %s with ligand morph %s failed: %s'
                           % (complex.mol_name, morph.name, why) )


    ### final message

    done_log.close()
    success = True

    for failed in ( (prot_failed, 'proteins'), (lig_failed, 'ligands'),
                    (com_failed, 'complexes'),
                    (morph_failed, 'morphs') ):
        if failed[0]:
            print ('ERROR: The following %s have failed:' % failed[1])

            for name in failed[0]:
                print (' %s' % name)

            success = False

    if success:
        print ('\n=== All molecules built successfully ===\n')


