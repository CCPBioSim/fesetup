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

r'''
Generates a UDFF file and associated PDB file for one molecule.

This writes a UDFF file for DL_FIELD. An UDFF file is a template for a molecule(s)
to be read by DL_FIELD which will generate from it a FIELD file for use in DL_POLY.
'''

__revision__ = "$Id$"



from collections import OrderedDict

import Sire.IO
import Sire.MM

from FESetup import const, errors, logger



def dlf_write(mol, postfix = '', pdb_name = const.LIGAND_NAME):
    """
    Extract topology and coordinate information from mol and convert to
    UDFF and PDB format.

    :param mol: molecule to be translated into UDFF/PDB
    :type mol: Sire.Mol.Molecule
    :param postfix: postfix for atom types (DL_FIELD is not case sensitive)
    :type postfix: string
    :param pdb_name: molecule PDB name
    :type pdb_name: string
    """

    try:
        params = mol.property('amberparameters')   # Sire.Mol.AmberParameters
    except UserWarning:
        raise errors.SetupError("molecule doesn't have any parameters")

    amber_types = OrderedDict()
    tot_natoms = mol.nAtoms()
    mw = 0

    Sire.IO.PDB().write(mol, const.DLFIELD_PDB_NAME)

    connects = OrderedDict()
    atom_info = []
    charges = []

    connectivity = mol.property('connectivity')

    for atom in mol.atoms():
        # NOTE: atom types are not case sensitive in DL_FIELD
        ambertype = str(atom.property('ambertype') )
        charge = atom.property('charge').value()
        mass = atom.property('mass').value()
        lj = atom.property('LJ')

        atom_name = str(atom.name().value() )
        mw += mass

        atom_info.append( (atom_name, ambertype) )
        charges.append(charge)

        try:
            element = atom.property('element').symbol()
        except UserWarning:
            element = 'XX'

        ambertype += postfix
        # DL_FIELD wants r* as in the original AMBER force field database
        amber_types[ambertype] = (element, mass,
                                  lj.sigma().value() * const.RSTAR_CONV,
                                  lj.epsilon().value())

        cl = []

        for con in connectivity.connectionsTo(atom.number() ):
            cl.append(str(mol.select(con).name().value() ) )

        connects[atom_name] = cl

    udff = open(const.DLFIELD_UDFF_NAME, 'w')

    udff.write('UNIT kcal/mol\nPOTENTIAL AMBER\n\n')

    udff.write('ATOM_TYPE\n')

    for typ in amber_types:
        udff.write('%-2s %s %s %.4f\n' %
                   (typ.replace(postfix, ''), typ,
                    amber_types[typ][0], amber_types[typ][1]) )

    udff.write('END ATOM_TYPE\n\n')

    udff.write('MOLECULE_TYPE\n')
    udff.write('ligand %s %.3f\n' % (pdb_name, mw) )
    udff.write('END MOLECULE_TYPE\n\n')


    # eliminate negative zero in total charge
    total_charge = float('%.3f' % sum(charges) ) + 0.0
    udff.write('MOLECULE ligand %i %.1f\n' % (mol.nAtoms(), total_charge) )

    for names, charge in zip(atom_info, charges):
        udff.write('%-2s %-2s %9.5f\n' % (names[0], names[1], charge) )

    for con in connects.iteritems():
        if len(con) > 1:
            udff.write('CONNECT %s > %s\n' %
                       (con[0], ' '.join('%s' % s for s in con[1:][0]) ) )


    improper_params = {}

    try:
        mol.property('improper')
    except UserWarning:
        pass
    else:
        improper_names = []
        improper_params = {}

        for dihedral in params.getAllImpropers():
            at0 = dihedral.atom0()
            at1 = dihedral.atom1()
            at2 = dihedral.atom2()
            at3 = dihedral.atom3()

            at0s = mol.select(at0)
            at1s = mol.select(at1)
            at2s = mol.select(at2)
            at3s = mol.select(at3)

            at0_name = at0s.name().value()
            at1_name = at1s.name().value()
            at2_name = at2s.name().value()
            at3_name = at3s.name().value()

            at0_type = str(at0s.property('ambertype') ) + postfix
            at1_type = str(at1s.property('ambertype') ) + postfix
            at2_type = str(at2s.property('ambertype') ) + postfix
            at3_type = str(at3s.property('ambertype') ) + postfix

            pk, pn, phase = params.getParams(dihedral)
            improper_names.append( (at0_name, at1_name, at2_name, at3_name) )
            improper_params[at0_type, at1_type, at2_type, at3_type] = pk, phase, pn

        for name in improper_names:
            udff.write('IMPROPER %s %s %s %s\n' % name)

    udff.write('END MOLECULE\n\n')

    try:
        mol.property('bond')
    except UserWarning:
        pass
    else:
        bond_info = {}

        for bond in params.getAllBonds():   # Sire.Mol.BondID
            at0 = bond.atom0()              # Sire.Mol.AtomIdx!
            at1 = bond.atom1()

            at0_type = str(mol.select(at0).property('ambertype') ) + postfix
            at1_type = str(mol.select(at1).property('ambertype') ) + postfix

            k, r0 = params.getParams(bond)
            bond_info[at0_type, at1_type] = k, r0

        udff.write('BOND\n')
        visited = set()

        for b in bond_info.keys():
            visited.add(b)

            if b[1] == b[0] or (b[1], b[0]) not in visited:
                udff.write('%s %s %.2f %.3f\n' %
                           (b[0], b[1], bond_info[b][0], bond_info[b][1]) )

        udff.write('END BOND\n\n')

    try:
        mol.property('angle')
    except UserWarning:
        pass
    else:
        angle_info = {}

        for angle in params.getAllAngles():  # Sire.Mol.AngleID
            at0 = angle.atom0()  # Sire.Mol.AtomIdx!
            at1 = angle.atom1()
            at2 = angle.atom2()

            at0_type = str(mol.select(at0).property('ambertype') ) + postfix
            at1_type = str(mol.select(at1).property('ambertype') ) + postfix
            at2_type = str(mol.select(at2).property('ambertype') ) + postfix

            k, theta0 = params.getParams(angle)
            angle_info[at0_type, at1_type, at2_type] = k, theta0 * const.RAD2DEG

        udff.write('ANGLE\n')
        visited = set()

        for a in angle_info.keys():
            visited.add(a)

            if a[0] == a[2] or (a[2], a[1], a[0]) not in visited:
                udff.write('%s %s %s %.2f %.3f\n' %
                           (a[0], a[1], a[2], angle_info[a][0], angle_info[a][1]) )

        udff.write('END ANGLE\n\n')

    try:
        mol.property('dihedral')
    except UserWarning:
        pass
    else:
        propers = {}
        intrascale = mol.property('intrascale')

        for dihedral in params.getAllDihedrals():  # Sire.Mol.DihedralID
            at0 = dihedral.atom0()  # Sire.Mol.AtomIdx!
            at1 = dihedral.atom1()
            at2 = dihedral.atom2()
            at3 = dihedral.atom3()

            at0_type = str(mol.select(at0).property('ambertype') ) + postfix
            at1_type = str(mol.select(at1).property('ambertype') ) + postfix
            at2_type = str(mol.select(at2).property('ambertype') ) + postfix
            at3_type = str(mol.select(at3).property('ambertype') ) + postfix

            at0idx = at0.value() + 1
            at3idx = at3.value() + 1

            p = params.getParams(dihedral)
            x = []

            n = 3
            for i in range(0, len(p), n):       # pk, np, phase
                x.append(p[i:i+n])

            # sorting should not be necessary
            #x.sort(key = lambda c: c[1], reverse = True)

            sf = intrascale.get(at0, at3)
            scee = sf.coulomb()
            scnb = sf.lj()

            propers[at0_type, at1_type, at2_type, at3_type] = x, scee, scnb

        # work around for DL_FIELD which doesn't like empty DIHEDRAL
        if propers:
            udff.write('DIHEDRAL\n')

            for a, t in propers.iteritems():
                npar = len(t[0])
                cnt = 0

                for dih in t[0]:
                    cnt += 1

                    if cnt < npar:
                        np = -dih[1]
                    else:
                        np = dih[1]

                    # IDIFV has no relevance in the topology file
                    udff.write('%s %s %s %s 1 '
                               '%.4f %.3f %i %.6f %.2f\n' %
                               (a[0], a[1], a[2], a[3],
                                dih[0], dih[2] * const.RAD2DEG, np, t[1], t[2]) )

            udff.write('END DIHEDRAL\n\n')

    if improper_params:
        udff.write('IMPROPER\n')

        for t, p in improper_params.iteritems():
            udff.write('%s %s %s %s %.3f %.3f %i\n' %
                            (t[0], t[1], t[2], t[3],
                             p[0], p[1] * const.RAD2DEG, p[2]) )


        udff.write('END IMPROPER\n\n')

    udff.write('VDW\n')

    for t in amber_types:
        udff.write('%s %.4f %.4f\n' %
                       (t, amber_types[t][2], amber_types[t][3]) )

    udff.write('END VDW\n\n')

    udff.close()
