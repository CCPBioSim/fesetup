#  Copyright (C) 2013-2017  Hannes H Loeffler
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
import sys

if ('%x' % sys.hexversion)[:3] != '207':
    print ('Python version 2.7 required, you have %i.%i' %
           sys.version_info[:2], file=sys.stderr)
    sys.exit(1)


from FESetup import _release

__revision__ = "$Id$"
__version__ = '0.8.3'

vstring = 'FESetup release %s, SUI version: %s' % (_release.release, __version__)
istring = ('Please cite: HH Loeffler, J Michel, C Woods, '
           'J Chem Inf Mod, 55, 2485\n'
           '             DOI 10.1021/acs.jcim.5b00368\n'
           'For more information please visit '
           'http://www.ccpbiosim.ac.uk/fesetup/')

import os
import argparse
import shutil
import glob
import copy
import atexit
import warnings
from collections import OrderedDict, namedtuple

import FESetup.prepare as prep
from FESetup import const, errors, create_logger, logger, DirManager
from FESetup.ui.iniparser import IniParser
from FESetup.modelconf import ModelConfig

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
    logger.write('\n%s\n\n%s\n' % (vstring, istring))
    atexit.register(lambda : logger.finalize() )

    ff_opts = list(opts[SECT_DEF]['forcefield'])
    lu = len(ff_opts)
    ff_opts[lu:] = deff[lu:]

    ff_opts.append(opts[SECT_DEF]['ff_addons'])
    ff_opts.append(opts[SECT_DEF]['mdengine'][0])
    ff_opts.append(opts[SECT_DEF]['parmchk_version'])
    ff_opts.append(opts[SECT_DEF]['gaff'])

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


def _search_for_model(names, workdir):
    """
    Check if a model file is available.

    :param names: list of model names to be searcher for, highest ranking
                  first
    :type names: str
    :param workdir: the workdir where the model file is located
    :type name: str
    :type workdir: str

    :returns: the path to the model file or None if model file not found
    """

    for name in names:
        path = os.path.join(workdir, name)

        if os.access(path, os.F_OK):
            return path

    return None


def read_model(filename):
    """
    Read a ModelConfig model.

    :param filename: name of file to be saved to
    :type filename: str

    :returns: a ModelConfig model
    """

    model = ModelConfig()
    model.read(filename)

    try:
        model.check_data(model['data.hash'], model['data.hash_type'])

        if not model['is.valid']:
            # FIXME: change exception type
            raise errors.SetupError('invalid model %s' % model['name'])
    except KeyError:
        raise errors.SetupError('invalid model %s' % model['name'])

    model.check_keys()

    return model


def save_model(model, mol, filename, dest_dir):
    """
    Save the current ModelConfig model.

    :param model: the model to be saved
    :type model: ModelConfig
    :param mol: the 'molecule' class
    :param mol: Common
    :param filename: name of file to be saved to
    :type filename: string
    :param dest_dir: name of destination directory
    :type dest_dir: string
    """

    model['mdengine'] = options[SECT_DEF]['mdengine']

    # add only latest info
    model['crd.filename'] = mol.amber_crd
    model['crd.filetype'] = 'amber-rst7'
    model['top.filename'] = mol.amber_top
    model['top.filetype'] = 'amber-parm7'

    model.add_file(mol.amber_crd)
    model.add_file(mol.amber_top)

    # FIXME: check for minimisation/equilibration
    if mol.box_dims:
        box_dims = []

        for v in mol.box_dims:
            box_dims.append(float(v))

        model['box.dimensions'] = box_dims
        model['box.format'] = 'bla'  # FIXME: boxlengths-angle

    if mol.ssbond_file and os.path.isfile(mol.ssbond_file):
        model['top.ssbond_file'] = mol.ssbond_file
        model.add_file(mol.ssbond_file)

    # the final blessing
    model['is.valid'] = 1

    model.write(filename)

    mname = os.path.join(dest_dir, filename)

    if os.access(mname, os.F_OK):
        os.remove(mname)

    shutil.move(filename, dest_dir)


def make_ligand(name, ff, opts):
    """
    Prepare ligands for simulation: charge parameters, vacuum top/crd,
    confomer search + alignment (both optional), optionally hydrated
    top/crd and minimisation/MD simulation.

    :param name: the name of the ligandx
    :type name: str
    :param ff: the ForceField class holding all relevant data for setup and
               also simulation
    :type ff: ForceField
    :param opts: the name of the ligandx
    :type opts: IniParser
    """

    logger.write('*** Working on %s ***\n' % name)

    lig = opts[SECT_LIG]

    load_cmds = ''

    if opts[SECT_DEF]['user_params']:
        load_cmds = _param_glob(PARAM_CMDS)

    vac_model_filename = name + const.MODEL_EXT
    sol_model_filename = 'solv_' + name + const.MODEL_EXT
    from_scratch = True

    workdir = os.path.join(os.getcwd(), const.LIGAND_WORKDIR, name)

    if not opts[SECT_DEF]['remake']:
        model_path = \
                   _search_for_model([sol_model_filename, vac_model_filename],
                                     const.LIGAND_WORKDIR)

        # FIXME: check for KeyError
        if model_path:
            model = read_model(model_path)
            name = model['name']

            # FIXME: only extract when const.LIGAND_WORKDIR not present?
            model.extract(direc = workdir)

            ligand = ff.Ligand(name)

            # invoke context manager to preset directory for absolute file
            # creation below
            with DirManager(workdir):
                logger.write('Found model %s, extracting data' % name)

                ligand.charge = float(model['charge.total'])
                ligand.gaff = model['forcefield']
                ligand.amber_top = model['top.filename']
                ligand.amber_crd = model['crd.filename']
                ligand.orig_file = model['crd.original']
                ligand.mol_file = ligand.orig_file
                ligand.mol_fmt ='mol2'

                # this file will not be created when skip_param = True
                try:
                    ligand.frcmod = model['frcmod']
                except KeyError:
                    pass

                if 'box.dimensions' in model:
                    ligand.box_dims = [float(b) for b in
                                       model['box.dimensions'].strip('[]')\
                                       .split(',')]

                if lig['morph.absolute'] and \
                       opts[SECT_DEF]['AFE.type'] == 'Sire':
                    logger.write('Creating input files for absolute '
                                 'transformations with Sire')
                    ligand.create_absolute_Sire()

            if os.path.basename(model_path) == vac_model_filename:
                from_scratch = False
            else:
                return ligand, load_cmds

    print('Making ligand %s...' % name)

    if not lig['basedir']:
        raise dGprepError('[%s] "basedir" must be set' % SECT_LIG)

    if not lig['file.format']:
        fmt = os.path.splitext(lig['file.name'])[1][1:]
    else:
        fmt = lig['file.format']

    if from_scratch:
        model = ModelConfig(name)
        ligand = ff.Ligand(name, lig['file.name'], fmt)

    # this file will not be created when skip_param = True
    if os.path.isfile(ligand.frcmod):
        model['frcmod'] = ligand.frcmod
        model.add_file(ligand.frcmod)

    if os.path.isabs(lig['basedir']):
        src = os.path.join(lig['basedir'], name)
    else:
        src = os.path.join(os.getcwd(), lig['basedir'], name)

    with DirManager(workdir):
        if from_scratch:
            ligand.copy_files((src,), None, opts[SECT_DEF]['overwrite'])

            if not os.access(lig['file.name'], os.F_OK):
                raise errors.SetupError('start file %s does not exist in %s' %
                                        (lig['file.name'], os.getcwd() ) )

            if lig['skip_param']:
                if fmt != 'pdb' and fmt != 'mol2':
                    raise dGprepError('When parameterisation is skipped, the input '
                                      'format must be PDB or MOL2')

                ligand.prepare('', lig['add_hydrogens'], lig['calc_charge'],
                               lig['correct_for_pH'], lig['pH'])
            elif not os.access(os.path.join(workdir, const.GAFF_MOL2_FILE),
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

            ligand.prepare_top()
            ligand.create_top(boxtype='', addcmd=load_cmds,
                              write_dlf=lig['write_dlf'])

            # this file will not be created when skip_param = True
            if os.path.isfile(const.LIGAND_AC_FILE):
                model['charge.filename'] = const.LIGAND_AC_FILE
                model.add_file(const.LIGAND_AC_FILE)

            model['charge.total'] = ligand.charge
            model['charge.filetype'] = 'ac'
            model['charge.method'] = 'AM1-BCC'
            model['forcefield'] = ligand.gaff
            model['molecule.type'] = 'ligand'

            model.add_file(ligand.mol_file)
            model['crd.original'] = ligand.mol_file

            save_model(model, ligand, vac_model_filename, '..')

            if opts[SECT_DEF]['MC_prep']:
                ligand.flex()

            nconf = lig['conf_search.numconf']

            if lig['morph.absolute'] and opts[SECT_DEF]['AFE.type'] == 'Sire':
                ligand.create_absolute_Sire()

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
            ligand.prepare_top()
            ligand.create_top(boxtype = lig['box.type'],
                              boxlength = lig['box.length'],
                              neutralize = lig['neutralize'],
                              addcmd = load_cmds, remove_first = False)

            if lig['ions.conc'] > 0.0:
                ligand.create_top(boxtype = lig['box.type'],
                                  boxlength = lig['box.length'],
                                  neutralize = 2,
                                  addcmd = load_cmds, remove_first = False,
                                  conc = lig['ions.conc'],
                                  dens = lig['ions.dens'])

            restr_force = lig['min.restr_force']
            nsteps = lig['min.nsteps']

            ligand.setup_MDEngine(opts[SECT_DEF]['mdengine'][1],
                                  opts[SECT_DEF]['mdengine.prefix'],
                                  opts[SECT_DEF]['mdengine.postfix'])

            if nsteps > 0:
                do_min(ligand, lig)

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

            save_model(model, ligand, sol_model_filename, '..')

    return ligand, load_cmds


def make_protein(name, ff, opts):
    """
    Prepare proteins for simulation.
    """

    logger.write('*** Working on %s ***\n' % name)

    prot = opts[SECT_PROT]

    load_cmds = ''

    if opts[SECT_DEF]['user_params']:
        load_cmds = _param_glob(PARAM_CMDS)

    vac_model_filename = name + const.MODEL_EXT
    sol_model_filename = 'solv_' + name + const.MODEL_EXT
    from_scratch = True

    workdir = os.path.join(os.getcwd(), const.PROTEIN_WORKDIR, name)

    if not opts[SECT_DEF]['remake']:
        model_path = _search_for_model([sol_model_filename, vac_model_filename],
                                       const.PROTEIN_WORKDIR)
        # FIXME: check for KeyError
        if model_path:
            model = read_model(model_path)
            name = model['name']

            logger.write('Found model %s, extracting data' % name)

            # FIXME: only extract when const.PROTEIN_WORKDIR not present?
            model.extract(direc = workdir)

            protein = ff.Protein(name, prot['basedir'])

            protein.charge = float(model['charge.total'])
            protein.amber_top = model['top.filename']
            protein.amber_crd = model['crd.filename']
            protein.orig_file = model['crd.original']

            if 'top.ssbond_file' in model:
                protein.ssbond_file = model['top.ssbond_file']

            if 'box.dimensions' in model:
                protein.box_dims = model['box.dimensions']

            if os.path.basename(model_path) == vac_model_filename:
                from_scratch = False
            else:
                return protein, load_cmds

    print('Making biomolecule %s...' % name)

    if not prot['basedir']:
        raise dGprepError('[%s] "basedir" must be set' % SECT_PROT)

    if os.path.isabs(prot['basedir']):
        src = os.path.join(prot['basedir'], name)
    else:
        src = os.path.join(os.getcwd(), prot['basedir'], name)

    if from_scratch:
        model = ModelConfig(name)
        protein = ff.Protein(name, prot['file.name'])
        model['crd.original'] = protein.mol_file
        model.add_file(protein.mol_file)

    with DirManager(workdir):
        if from_scratch:
            protein.copy_files((src,), None, opts[SECT_DEF]['overwrite'])

            if prot['propka']:
                protein.protonate_propka(pH = prot['propka.pH'])

            protein.get_charge()    # must be done explicitly
            protein.prepare_top()
            protein.create_top(boxtype = '')

            model['charge.total'] = protein.charge
            model['forcefield'] = 'AMBER'    # FIXME
            model['molecule.type'] = 'biomolecule'

            save_model(model, protein, vac_model_filename, '..')

        # FIXME: also check for boxlength and neutralize

        if prot['box.type']:
            protein.prepare_top()
            protein.create_top(boxtype = prot['box.type'],
                               boxlength = prot['box.length'],
                               neutralize = prot['neutralize'],
                               align = prot['align_axes'],
                               addcmd = load_cmds, remove_first = True)

            if prot['ions.conc'] > 0.0:
                protein.create_top(boxtype = prot['box.type'],
                                   boxlength = prot['box.length'],
                                   neutralize = 2,
                                   align = prot['align_axes'],
                                   addcmd = load_cmds, remove_first = False,
                                   conc = prot['ions.conc'],
                                   dens = prot['ions.dens'])

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

            save_model(model, protein, sol_model_filename, '..')

    return protein, load_cmds


def make_complex(prot, lig, ff, opts, load_cmds):

    com = opts[SECT_COM]

    name = prot.mol_name + const.PROT_LIG_SEP + lig.mol_name
    vac_model_filename = name + const.MODEL_EXT
    sol_model_filename = 'solv_' + name + const.MODEL_EXT
    from_scratch = True

    workdir = os.path.join(os.getcwd(), const.COMPLEX_WORKDIR, name)

    if not opts[SECT_DEF]['remake']:
        model_path = _search_for_model([sol_model_filename, vac_model_filename],
                                       const.COMPLEX_WORKDIR)

        if model_path:
            model = read_model(model_path)
            name = model['name']

            logger.write('Found model %s, extracting data' % name)

            # FIXME: only extract when const.COMPLEX_WORKDIR not present?
            model.extract(direc = workdir)

            complex = ff.Complex(prot, lig)

            complex.charge = float(model['charge.total'])
            complex.amber_top = model['top.filename']
            complex.amber_crd = model['crd.filename']

            if 'top.ssbond_file' in model:
                complex.ssbond_file = model['top.ssbond_file']

            if 'box.dimensions' in model:
                complex.box_dims = model['box.dimensions']

            if os.path.basename(model_path) == vac_model_filename:
                from_scratch = False
            else:
                return complex, load_cmds

    print('Making complex from %s and %s...' % (prot.mol_name, lig.mol_name))

    if from_scratch:
        model = ModelConfig(name)

    if not model_path:
        complex = ff.Complex(prot, lig)

    lig_src = os.path.join(os.getcwd(), const.LIGAND_WORKDIR, lig.mol_name)
    prot_src = os.path.join(os.getcwd(), const.PROTEIN_WORKDIR, prot.mol_name)

    with DirManager(workdir):
        if from_scratch:
            complex.copy_files((lig_src, prot_src),
                               (ligand.orig_file, ligand.frcmod, protein.orig_file,
                                const.LIGAND_AC_FILE, const.SSBOND_FILE),
                               opts[SECT_DEF]['overwrite'])

            complex.ligand_fmt = lig.mol_fmt
            complex.prepare_top(gaff=options[SECT_DEF]['gaff'])
            complex.create_top(boxtype='', addcmd=load_cmds)

            model['name'] = complex.complex_name
            model['charge.total'] = complex.charge
            model['forcefield'] = 'AMBER'   # FIXME
            model['molecule.type'] = 'complex'  # FIXME

            save_model(model, complex, vac_model_filename, '..')

        # FIXME: also check for boxlength and neutralize
        if com['box.type']:
            complex.prepare_top(gaff=options[SECT_DEF]['gaff'])
            complex.create_top(boxtype=com['box.type'],
                               boxlength=com['box.length'],
                               neutralize=com['neutralize'],
                               align=com['align_axes'],
                               addcmd=load_cmds, remove_first = True)

            if com['ions.conc'] > 0.0:
                complex.create_top(boxtype=com['box.type'],
                                   boxlength=com['box.length'],
                                   neutralize=2,
                                   align=com['align_axes'],
                                   addcmd=load_cmds, remove_first=False,
                                   conc=com['ions.conc'],
                                   dens=com['ions.dens'])

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

            save_model(model, complex, sol_model_filename, '..')

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
    'gaff': ('gaff1', None),
    'mdengine': (['amber', 'sander'], ('list', LIST_SEP) ),
    'mdengine.prefix': ('', None),
    'mdengine.postfix': ('', None),
    'parmchk_version': (2, (int, ) ),
    'FE_type': ('', None),
    'AFE.type': ('Sire', None),
    'AFE.separate_vdw_elec': (True, ('bool', ) ),
    'softcore_type': ('ignored', None),
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
    'morph.absolute': (False, ('bool', ) ),
    'box.type': ('', None),
    'box.length': (10.0, (float, ) ),
    'neutralize': (False, ('bool', ) ),
    'ions.conc': (0.0, (float, ) ),
    'ions.dens': (1.0, (float, ) ),
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
    'skip_param': (False, ('bool', ) ),
    'write_dlf': (False, ('bool', ) )
    }

defaults[SECT_PROT] = {
    'basedir': ('', None),
    'file.name': ('protein.pdb', None),
    'molecules': ('', ('list', LIST_SEP) ),
    'box.type': ('', None),
    'box.length': (10.0, (float, ) ),
    'neutralize': (False, ('bool', ) ),
    'ions.conc': (0.0, (float, ) ),
    'ions.dens': (1.0, (float, ) ),
    'align_axes': (False, ('bool', ) ),
    'propka': (False, ('bool', ) ),
    'propka.pH': (7.0, (float, ) ),
    }

defaults[SECT_COM] = {
    'pairs': ('', ('pairlist', LIST_SEP, COM_PAIR_SEP) ),
    'box.type': ('', None),
    'box.length': (10.0, (float, ) ),
    'neutralize': (False, ('bool', ) ),
    'ions.conc': (0.0, (float, ) ),
    'ions.dens': (1.0, (float, ) ),
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
    parser.add_argument('infile', nargs='?',
                        help='input file in INI format, if not given then '
                        'just output defaults')
    parser.add_argument('-v', '--version', action='version',
                        version=vstring,
                        help='full version information')
    parser.add_argument('--tracebacklimit', metavar='N', type=int, default=0,
                        help='set the Python traceback limit (for debugging)')
    args = parser.parse_args()

    print('\n=== %s ===\n\n%s\n' % (vstring, istring))

    options = IniParser(copy.deepcopy(defaults))

    if not args.infile:
        print('\n'.join(options.format()))
        sys.exit(0)

    sys.tracebacklimit = args.tracebacklimit

    options.parse(args.infile, 'globals')

    # backward compatibility
    if options[SECT_DEF]['FE_type']:
        options[SECT_DEF]['AFE.type'] = options[SECT_DEF]['FE_type']

    if options[SECT_DEF]['gaff'] == 'gaff1':
        options[SECT_DEF]['gaff'] = 'gaff'


    for section in ALL_SECTIONS:
        check_dict(section, options)

    ff = prelude(options)

    logger.write('Command line: %s\n\nOptions:\n--------' % ' '.join(sys.argv))
    logger.write('\n'.join(options.format()))
    logger.write('--------\n\nForce field and MD engine:\n%s\n' % ff)


    # FIXME: We keep all molecule objects in memory.  For 2000 morph pairs
    #        this may mean more than 1 GB on a 64 bit machine.

    ### proteins

    proteins = {}
    prot_failed = []

    for prot_name in options[SECT_PROT]['molecules']:
        try:
            protein, cmds = make_protein(prot_name, ff, options)
            proteins[protein] = cmds
        except errors.SetupError as why:
            prot_failed.append(prot_name)
            print ('ERROR: %s failed: %s' % (prot_name, why))


    ### ligands

    ligands = {}
    Ligdata = namedtuple('Ligdata', ['ref', 'leapcmd'])
    lig_failed = []

    morph_pairs = copy.deepcopy(options[SECT_LIG]['morph_pairs'])
    molecules = copy.deepcopy(options[SECT_LIG]['molecules'])
    morph_maps = {}

    if morph_pairs:
        temp_pairs = []

        # FIXME: better error handling, parsing through IniParser?
        for pair in morph_pairs:
            p1 = pair[1].split('/')
            temp_map = {}

            if len(p1) > 1:
                p1[0] = p1[0].strip()

                for idx in p1[1:]:
                    if idx:
                        if idx.startswith('!'):
                            a = idx[1:]
                            b = -1
                        else:
                            a, b = idx.split('=')

                        try:
                            temp_map[int(a)] = int(b)
                        except ValueError:
                            print('Error: map contains non-integers')
                            sys.exit(1)

                morph_maps[pair[0],p1[0]] = temp_map

            temp_pairs.append( (pair[0], p1[0]) )

        morph_pairs = temp_pairs

        mols = [val for pairs in morph_pairs for val in pairs]  # flatten
        uniq = list(OrderedDict( (val, None) for val in mols) )
        molecules = uniq

    for lig_name in molecules:
        try:
            ligand, cmds = make_ligand(lig_name, ff, options)
            ligands[lig_name] = Ligdata(ligand, cmds)
        except errors.SetupError as why:
            lig_failed.append(lig_name)
            print('ERROR: %s failed: %s' % (lig_name, why))


    ### ligand morphs

    if morph_pairs:
        print('Morphs will be generated for %s' % options[SECT_DEF]['AFE.type'])
        logger.write('Morphs will be generated for %s\n' %
                     options[SECT_DEF]['AFE.type'])

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

        isotope_map = {}

        if (pair[0], pair[1]) in morph_maps:
            isotope_map = morph_maps[pair[0], pair[1]]

        basedir = os.path.join(os.getcwd(), options[SECT_LIG]['basedir'])
        wd1 = os.path.join(os.getcwd(), const.LIGAND_WORKDIR, pair[0])
        wd2 = os.path.join(os.getcwd(), const.LIGAND_WORKDIR, pair[1])

        with mutate.Morph(ligand1, ligand2, wd1, wd2, ff,
                          options[SECT_DEF]['AFE.type'],
                          options[SECT_DEF]['AFE.separate_vdw_elec'],
                          options[SECT_DEF]['mcs.timeout'],
                          options[SECT_DEF]['mcs.match_by'],
                          options[SECT_DEF]['gaff']) as morph:

            print ('Morphing %s to %s...' % pair)

            if (pair[1], pair[0]) in morph_pairs:
                rev = ligand2
            else:
                rev = None

            try:
                morph.setup(cmd1, cmd2, basedir, isotope_map)

                if options[SECT_LIG]['box.type']:
                    morph.create_coords(ligand1, 'solvated', wd1,
                                        cmd1, cmd2, rev, wd2)
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

                try:
                    complex, cmds = make_complex(protein,
                                                 ligands[lig_name].ref,
                                                 ff, options, cmds)

                    # NOTE: no real need here to limit the list of complexes
                    #       as it is checked again below
                    for pair in morph_pairs:
                        if complex.ligand.mol_name == pair[0]:
                            complexes[complex] = cmds
                except errors.SetupError as why:
                    com_failed.append(name)
                    print ('ERROR: %s failed: %s' % (name, why))


    ### complex morphs

    # FIXME: For some reason the cwd is _perturbations/type/morph for the
    #        second morph in the loop below as if a chdir didn't take place(?).
    #        So remember the cwd before entering create_coords()
    cwd = os.getcwd()

    for morph in morphs:
        for complex in complexes.keys():
            name = complex.mol_name + '/' + morph.name

            # FIXME: Complex has no ligand component after restart
            if complex.ligand.mol_name == morph.initial_name:
                print('Creating complex %s with ligand morph %s...' %
                      (complex.mol_name, morph.name) )

                #print(os.getcwd())
                wd = os.path.join(cwd, const.COMPLEX_WORKDIR, complex.mol_name)

                try:
                    morph.create_coords(complex, 'complex', wd,
                                        complexes[complex], '')
                except errors.SetupError as why:
                    morph_failed.append(name)
                    print ('ERROR: complex %s with ligand morph %s failed: %s'
                           % (complex.mol_name, morph.name, why) )


    ### final message

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
