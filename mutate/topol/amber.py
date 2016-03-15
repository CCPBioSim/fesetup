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
Common AMBER support.
"""


__revision__ = "$Id$"


import Sire.Mol
import Sire.Units

from FESetup import const
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


    for iinfo, finfo in atom_map.items():
        istr = iinfo.name.value()
        fstr = finfo.name.value()

        new = state1.atom(iinfo.index)  # AtomEditor

        # Dummy may be bonded to H which exceeds leap's valency limit.  Leap
        # will only keep the first bond encountered.  Solution:
        # "pert=true" and explicit bonding to DU, recreate all bonds to not
        # rely on which one leap had kept
        if not iinfo.atom:
            for idx in con_morph.connectionsTo(iinfo.index):
                atom = lig_morph.select(idx)
                name = '%s' % atom.name().value()
                ambertype = '%s' % atom.property('ambertype')

                if ambertype.upper().startswith('H') and \
                       name.startswith('H'):
                    pert0_info.append( (str(istr), str(name)) )

        if fstr.startsWith('H') and con_morph.nConnections(iinfo.index) > 1:
            for bond_index in con_morph.connectionsTo(iinfo.index):
                atom1 = lig_morph.select(bond_index)
                rname = util.search_atominfo(atom1.index(), atom_map)
                name = '%s' % rname.name.value()
                pert1_info.append( (str(fstr), str(name)) )


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
