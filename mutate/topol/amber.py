#  Copyright (C) 2016  Hannes H Loeffler
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
AMBER support functions for alchemical free energy simulations with either
pmemd or sander.
"""


__revision__ = "$Id$"


import Sire.Mol
import Sire.Units

from FESetup import const, errors
from FESetup.mutate import util



def dummy(lig_morph, con_morph, lig_final, atom_map):
    """
    Support for AMBER/dummy which creates information for leap to deal with
    its valency check.  Also, create the final state from the morph.

    :param lig_morph: the morph molecule
    :type lig_morph: Sire.Mol.CutGroup
    :param con_morph: the connectivity of the morph
    :type con_morph: Sire.Mol.Connectivity
    :param lig_final: the final state molecule
    :type lig_final: Sire.Mol.Molecule
    :param atom_map: the forward atom map
    :type atom_map: dict of _AtomInfo to _AtomInfo
    :returns: initial state molecule, final state molecule
    :rtype: Sire.Mol.Molecule, Sire.Mol.Molecule
    """

    state1_m = Sire.Mol.Molecule(lig_morph)
    state1 = state1_m.edit()    # MolEditor

    pert0_info = []
    pert1_info = []
    H_bonded = set()

    for iinfo, finfo in atom_map.items():
        istr = iinfo.name.value()
        fstr = finfo.name.value()

        new = state1.atom(iinfo.index)  # AtomEditor

        # Dummies bonded to H exceed leap's valency limit of 2 and leap will
        # only keep the first bond encountered.  Solution: "pert=true" and
        # recreate _ALL_ bonds to H.  Do not rely on which one leap had kept
        # so delete all bonds first because leap terminates on prmtop writing
        # when the one existing bond had been defined twice.

        if istr.startsWith('H') and con_morph.nConnections(iinfo.index) > 1:
            for bond_index in con_morph.connectionsTo(iinfo.index):
                atom = lig_morph.select(bond_index)
                name = '%s' % atom.name().value()
                pert0_info.append((str(istr), str(name)))

        if fstr.startsWith('H') and con_morph.nConnections(iinfo.index) > 1:
            for bond_index in con_morph.connectionsTo(iinfo.index):
                atom1 = lig_morph.select(bond_index)
                rname = util.search_atominfo(atom1.index(), atom_map)
                name = '%s' % rname.name.value()
                pert1_info.append((str(fstr), str(name)))

        if not finfo.atom:
            charge = 0.0 * Sire.Units.mod_electron
            ambertype = const.DUMMY_TYPE
        else:
            base = lig_final.atoms().select(finfo.index)

            charge = base.property('charge')
            ambertype = '%s' % base.property('ambertype')

            new.setProperty('element', Sire.Mol.Element(istr) )

        # tag atom name because Sire throws exception when duplicate atom names
        new.rename(Sire.Mol.AtomName('\x00%s\x00' % fstr ) )

        new.setProperty('ambertype', ambertype)
        new.setProperty('charge', charge)

        # NOTE: state1 is already MolEditor but if we do not call molecule() on
        #       it, Sire segfaults...
        state1 = new.molecule()

    # second pass to remove atom name tags
    for atom in state1.molecule().atoms():
        new = state1.atom(atom.index() )

        name = atom.name().value().replace('\x00', '')

        new.rename(Sire.Mol.AtomName(name) )
        state1 = new.molecule()

    state1 = state1.commit()

    return state1, pert0_info, pert1_info


def softcore(lig_morph, lig_final, atom_map):
    """
    Support for AMBER/softcore which deals with rearranging atoms in the
    final state to correspond to the initial state and deletes dummies.  The
    final state is created here.

    :param lig_morph: the morph molecule
    :type lig_morph: Sire.Mol.CutGroup
    :param lig_final: the final state molecule
    :type lig_final: Sire.Mol.Molecule
    :param atom_map: the forward atom map
    :type atom_map: dict of _AtomInfo to _AtomInfo
    :returns: initial state molecule, final state molecule
    :rtype: Sire.Mol.Molecule, Sire.Mol.Molecule
    """

    state1_m = Sire.Mol.Molecule(lig_morph)
    state1 = state1_m.edit()            # MolEditor

    idels = []
    fdels = []

    for iinfo, finfo in atom_map.iteritems():
        istr = iinfo.name.value()
        fstr = finfo.name.value()
        iidx = iinfo.index

        if not iinfo.atom:
            # when indices are deleted they will be recomputed
            idels.insert(0, iidx)

        if not finfo.atom:
            fdels.insert(0, iidx)
        else:
            new = state1.atom(iidx)     # AtomEditor

            base = lig_final.atoms().select(finfo.index)
            charge = base.property('charge')
            ambertype = '%s' % base.property('ambertype')

            # tag atom names because Sire throws exception for duplicates
            new.rename(Sire.Mol.AtomName('\x00%s\x00' % fstr) )

            new.setProperty('ambertype', ambertype)
            new.setProperty('charge', charge)
            new.setProperty('element', Sire.Mol.Element(istr) )

            state1 = new.molecule()     # MolEditor

    # second pass to remove atom name tags
    for atom in state1.molecule().atoms():
        new = state1.atom(atom.index() )

        name = atom.name().value().replace('\x00', '')

        new.rename(Sire.Mol.AtomName(name) )
        state1 = new.molecule()


    state0_m = Sire.Mol.Molecule(lig_morph)
    state0 = state0_m.edit()

    if idels:
        for idx in idels:
            state0 = state0.remove(idx).commit().edit()

    if fdels:
        state0 = state0.edit()
        max_idx = len(atom_map) - len(idels)
        mxi = Sire.Mol.AtomIdx(max_idx)

        for idx in fdels:
            state0 = state0.atom(idx).reindex(mxi).molecule()
            state1 = state1.remove(idx).commit().edit()

    return state0.commit(), state1.commit()


ONESTEP_MDIN = 'onestep%s.in'
DECHARGE_MDIN = 'decharge%s.in'
VDW_MDIN = 'vdw%s.in'
RECHARGE_MDIN = 'recharge%s.in'
CHARGE_MDIN = 'charge%s.in'
GROUP_FILE = '%s.group'

PMEMD_TEMPLATE = """\
timask1='{timask1}', timask2='{timask2}',
scmask1='{scmask1}',
scmask2='{scmask2}',"""

SANDER_TEMPLATE = "scmask='{scmask}',"

COMMON_TEMPLATE = '''TI/FEP, NpT, {title}
 &cntrl
 ! please adjust namelist parameters to your needs!

 ! parameters for general MD
 imin = 0, nstlim = 500000, irest = %%R1%%, ntx = %%R2%%, dt = {dt},
 ntt = 3, temp0 = 298.0, gamma_ln = 2.0, ig = -1,
 ntb = {ntb},
 {press},
 ntwe = 10000, ntwx = 10000, ntpr = 500, ntwr = 50000, ntave = 50000,
 ioutfm = 1, iwrap = {wrap}, ntxo = 2,

 ! parameters for alchemical free energy simulation
 ntc = 2, ntf = 1,
 noshakemask = '{noshakemask}',

 icfe = 1, ifsc = {ifsc}, clambda = %%L%%, scalpha = 0.5, scbeta = 12.0,
 ifmbar = 1, bar_intervall = 500, bar_l_min = 0.0, bar_l_max = 1.0,
 bar_l_incr = 0.1,   ! ntpr = bar_intervall for alchemical analysis tool

 %s
 crgmask = '{crgmask}',
 /
 &ewald
 /
'''

GROUP_FILE_TEMPLATE = '''\
# this is a template group file for multi-sander, edit output files as needed
-i {mdin_a} -p {base_a}.parm7 -c {base_a}.rst7 -O -inf ti001_0.info -o ti001_0.out -r ti001_0.rst7 -e ti001_0.en -x ti001_0.nc
-i {mdin_b} -p {base_b}.parm7 -c {base_b}.rst7 -O -inf ti001_1.info -o ti001_1.out -r ti001_1.rst7 -e ti001_1.en -x ti001_1.nc
'''

#FIXME: one vs two topology files
def write_mdin(atoms_initial, atoms_final, atom_map, prog, style='', vac=True):
    """
    Create mdin input file(s) with proper masks.

    :param atoms_initial: set of initial atoms
    :type atoms_initial: Sire.Mol.Selector_Atom
    :param atoms_final: set of final atoms
    :type atoms_final: Sire.Mol.Selector_Atom
    :param atom_map: the forward atom map
    :type atom_map: dict of _AtomInfo to _AtomInfo
    :param style: softcoreN or dummyN
    :type style: str
    :param prog: pmemd or sander
    :type prog: str
    :param vac: create vacuum input file
    :type vac: bool
    :raises: SetupError
    """

    mask0 = []
    mask1 = []
    dummies0 = False
    dummies1 = False
    ifsc = 1
    dt = 0.002

    for iinfo, finfo in atom_map.items():
        # FIXME: also set noshakemask if dt is set to smaller value
        if iinfo.atom and finfo.atom:
            el0 = atoms_initial.select(iinfo.index).property('element').symbol()
            el1 = atoms_final.select(finfo.index).property('element').symbol()

            if (el0 == 'H' and el1 != 'H') or (el1 == 'H' and el0 != 'H'):
                dt = 0.001

        # we should never have a dummy to dummy mapping
        if not iinfo.atom:
            dummies0 = True
            mask1.append(str(finfo.name.value()))
        elif not finfo.atom:
            dummies1 = True
            mask0.append(str(iinfo.name.value()))

    # softcores are not necessary if there are no appearing/disappearing atoms
    if not mask0 and not mask1:
        ifsc = 0

    if vac:
        ntb = 0
        press = 'ntp = 0, cut = 9999.0, igb = 6, ! igb for pmemd'
        wrap = 0
    else:
        ntb = 2
        press = 'ntp = 1, barostat = 1, pres0 = 1.01325, taup = 2.0'
        wrap = 1

    add_str0 = ''
    add_str1 = ''

    if 'dummy' in style and dummies0:
        add_str0 = 'DU='

    if 'dummy' in style and dummies1:
        add_str1 = 'DU='

    mode = 'w'

    if style == 'softcore' or style == 'dummy':
        mask_str0 = ','.join(mask0)
        mask_str1 = ','.join(mask1)

        title='1-step transformation' if ifsc == 1 else \
               'linear transformation only',

        if mask_str0:
            mask_str0 = ':1@%s%s' % (mask_str0, ',' + add_str0 if add_str0
                                     else '')
        elif add_str0:
            mask_str0 = ':1@%s' % add_str0

        if mask_str1:
            mask_str1 = ':2@%s%s' % (mask_str1, ',' + add_str1 if add_str1
                                     else '')
        elif add_str1:
            mask_str1 = ':2@%s' % add_str1

        if prog == 'pmemd':
            tmpl = COMMON_TEMPLATE % PMEMD_TEMPLATE

            with open(ONESTEP_MDIN % '', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title,
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        timask1=':1', timask2=':2',
                        noshakemask=':1,2', crgmask='', ifsc=ifsc,
                        scmask1=mask_str0, scmask2=mask_str1))
        else:
            tmpl = COMMON_TEMPLATE % SANDER_TEMPLATE

            for filen, mask_str in ( (ONESTEP_MDIN % '_a', mask_str0),
                                     (ONESTEP_MDIN % '_b',
                                      mask_str1.replace(':2', ':1', 1))):
                with open(filen, mode) as stfile:
                    stfile.write(
                        tmpl.format(
                            title=title,
                            dt=dt, ntb=ntb, press=press, wrap=wrap,
                            noshakemask=':1',
                            crgmask='', ifsc=ifsc,
                            scmask=mask_str))

            with open(GROUP_FILE % 'onestep', mode) as gfile:
                gfile.write(
                    GROUP_FILE_TEMPLATE.format(
                        mdin_a=ONESTEP_MDIN % '_a', mdin_b=ONESTEP_MDIN % '_b',
                        base_a='state0', base_b='state1'))

    elif style == 'softcore2' or style == 'dummy2':
        mask_str0 = ','.join(mask0)
        mask_str1 = ','.join(mask1)

        if mask_str0:
            mask_str0 = ':1@%s%s' % (mask_str0, ',' + add_str0 if add_str0
                                     else '')
        elif add_str0:
            mask_str0 = ':1@%s' % add_str0

        if mask_str1:
            mask_str1 = ':2@%s%s' % (mask_str1, ',' + add_str1 if add_str1
                                     else '')
        elif add_str1:
            mask_str1 = ':2@%s' % add_str1

        if dummies1:
            ifsc1, ifsc2 = 0, 1
            m0, m1, m2, m3 = '', '', mask_str0, mask_str1
            step1_filename = CHARGE_MDIN
            step1_name = 'charge%s'
            step2_filename = VDW_MDIN
            step2_name = 'vdw%s'
            title1 = 'charge transformation'
            title2 = 'vdW+bonded transformation'
        else:
            ifsc1, ifsc2 = 1, 0
            m0, m1, m2, m3 = mask_str0, mask_str1, '', ''
            step1_filename = VDW_MDIN
            step1_name = 'vdw%s'
            step2_filename = CHARGE_MDIN
            step2_name = 'charge%s'
            title1 = 'vdW+bonded transformation'
            title2 = 'charge transformation'

        if prog == 'pmemd':
            tmpl = COMMON_TEMPLATE % PMEMD_TEMPLATE

            with open(step1_filename % '' , mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title1,
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        timask1=':1', timask2=':2',
                        noshakemask=':1,2', crgmask='', ifsc=ifsc1,
                        scmask1=m0, scmask2=m1))

            with open(step2_filename % '', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title2,
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        timask1=':1', timask2=':2',
                        noshakemask=':1,2', crgmask='', ifsc=ifsc2,
                        scmask1=m2, scmask2=m3))
        else:
            tmpl = COMMON_TEMPLATE % SANDER_TEMPLATE

            with open(step1_filename % '_a' , mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title1 + ', step a',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        noshakemask=':1', crgmask='',
                        ifsc=ifsc1, scmask=m0))

            with open(step1_filename % '_b' , mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title1 + ', step b',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        noshakemask=':1', crgmask='',
                        ifsc=ifsc1, scmask=m1.replace(':2', ':1', 1)))

            fname = step1_name.replace('%s', '')
            with open(GROUP_FILE % fname, mode) as gfile:
                gfile.write(
                    GROUP_FILE_TEMPLATE.format(
                        mdin_a=step1_name % '_a', mdin_b=step1_name % '_b',
                        base_a='state0', base_b='state_int'))

            with open(step2_filename % '_a', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title2 + ', step a',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        noshakemask=':1', crgmask='',
                        ifsc=ifsc2, scmask=m2))

            with open(step2_filename % '_b', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title=title2 + ', step b',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        noshakemask=':1', crgmask='',
                        ifsc=ifsc2, scmask=m3.replace(':2', ':1', 1)))

            fname = step2_name.replace('%s', '')
            with open(GROUP_FILE % fname, mode) as gfile:
                gfile.write(
                    GROUP_FILE_TEMPLATE.format(
                        mdin_a=step2_name % '_a', mdin_b=step2_name % '_b',
                        base_a='state_int', base_b='state1'))

    elif style == 'softcore3' or style == 'dummy3':
        mask_str0 = ','.join(mask0)
        mask_str1 = ','.join(mask1)

        if add_str0:
            add_str0 = ',%s' % add_str0

        if add_str1:
            add_str1 = ',%s' % add_str1

        if prog == 'pmemd':
            tmpl = COMMON_TEMPLATE % PMEMD_TEMPLATE

            # FIXME: partial de/recharging with scmask for crgmask?
            with open(DECHARGE_MDIN % '', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='decharge transformation',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        timask1=':1', timask2=':2',
                        crgmask=':2', noshakemask=':1,2',
                        ifsc=0, scmask1='', scmask2=''))

            with open(VDW_MDIN % '', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='vdW+bonded transformation',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        timask1=':1', timask2=':2',
                        crgmask=':1,2', noshakemask=':1,2', ifsc=ifsc,
                        scmask1=':1@%s' % mask_str0 + add_str0,
                        scmask2=':2@%s' % mask_str1 + add_str1))

            with open(RECHARGE_MDIN % '', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='recharge transformation',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        timask1=':1', timask2=':2',
                        crgmask=':1', noshakemask=':1,2',
                        ifsc=0, scmask1='', scmask2=''))
        else:
            tmpl = COMMON_TEMPLATE % SANDER_TEMPLATE

            with open(DECHARGE_MDIN % '_a', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='decharge transformation, step a',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        crgmask='', noshakemask=':1',
                        ifsc=0, scmask=''))

            with open(DECHARGE_MDIN % '_b', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='decharge transformation, step b',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        crgmask=':1', noshakemask=':1',
                        ifsc=0, scmask=''))

            with open(GROUP_FILE % 'decharge', mode) as gfile:
                gfile.write(
                    GROUP_FILE_TEMPLATE.format(
                        mdin_a=DECHARGE_MDIN % '_a',
                        mdin_b=DECHARGE_MDIN % '_b',
                        base_a='state0', base_b='state0'))

            with open(VDW_MDIN % '_a', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='vdW+bonded transformation, step a',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        crgmask=':1', noshakemask=':1',
                        ifsc=ifsc, scmask=':1@%s' % mask_str0 + add_str0))

            with open(VDW_MDIN % '_b', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='vdW+bonded transformation, step b',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        crgmask=':1', noshakemask=':1',
                        ifsc=ifsc, scmask=':1@%s' % mask_str1 + add_str1))

            with open(GROUP_FILE % 'vdw', mode) as gfile:
                gfile.write(
                    GROUP_FILE_TEMPLATE.format(
                        mdin_a=VDW_MDIN % '_a', mdin_b=VDW_MDIN % '_b',
                        base_a='state0', base_b='state1'))

            with open(RECHARGE_MDIN % '_a', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='recharge transformation, step a',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        crgmask=':1', noshakemask=':1',
                        ifsc=0, scmask=''))

            with open(RECHARGE_MDIN % '_b', mode) as stfile:
                stfile.write(
                    tmpl.format(
                        title='recharge transformation, step b',
                        dt=dt, ntb=ntb, press=press, wrap=wrap,
                        crgmask='', noshakemask=':1',
                        ifsc=0, scmask=''))

            with open(GROUP_FILE % 'recharge', mode) as gfile:
                gfile.write(
                    GROUP_FILE_TEMPLATE.format(
                        mdin_a=RECHARGE_MDIN % '_a',
                        mdin_b=RECHARGE_MDIN % '_b',
                        base_a='state1', base_b='state1'))
    else:
        raise errors.SetupError('Unknown FE type: %s' % style)
