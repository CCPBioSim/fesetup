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
Generates a Gromacs topology file and optionally a gro coordinate file from a
AMBER parmtop and inpcrd.
'''

__revision__ = "$Id$"



import os, math
from collections import OrderedDict

import Sire.IO
import Sire.MM

from FESetup import const, errors, logger



ATOM_PREFIX = 'x'
BOX_BUFFER = 3.0

water_atom_names = {
    'O': 'OW',
    'H1': 'HW1',
    'H2': 'HW2',
    'EPW': 'EP'
    }

default_header = '''
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333
'''

WATER_POSRES = const.GROMACS_POSRES_PREFIX + 'WAT' + const.GROMACS_ITP_EXT

TIP3P_header = '''
[ moleculetype ]
WAT 2

[ atoms ]
1 OW 1 WAT OW  1 -0.834 16.000
2 HW 1 WAT HW1 1  0.417  1.008
3 HW 1 WAT HW2 1  0.417  1.008

#ifndef FLEXIBLE
[ settles ]
1 1 0.09572 0.15139

[ exclusions ]
1 2 3
2 1 3
3 1 2
#else
[ bonds ]
1 2 1 0.09572 502416.0
1 3 1 0.09572 502416.0

[ angles ]
2 1 3 1 104.52 628.02
#endif

#ifdef POSRES
#include "%s"
#endif
''' % WATER_POSRES

TIP4PEW_header = '''
[ moleculetype ]
WAT 2

[ atoms ]
1 OW 1 WAT OW  1 0.0     16.00000
2 HW 1 WAT HW1 1 0.52422  1.00800
3 HW 1 WAT HW2 1 0.52422  1.00800
4 EP 1 WAT EP  1 -1.04844 0.00000

#ifndef FLEXIBLE
[ settles ]
1 1 0.09572 0.15139
#else
[ bonds ]
1 2 1 0.09572 502416.0
1 3 1 0.09572 502416.0
 
[ angles ]
2 1 3 1 104.52 628.02
#endif

[ virtual_sites3 ]
4 1 2 3 1 0.106676721 0.106676721

[ exclusions ]
1 2 3 4
2 1 3 4
3 1 2 4
4 1 2 3

#ifdef POSRES
#include "%s"
#endif
''' % WATER_POSRES


SPCE_header = '''
[ moleculetype ]
WAT 2

[ atoms ]
1 OW 1 WAT OW  1 -0.8476 15.9994
2 HW 1 WAT HW1 1  0.4238  1.008
3 HW 1 WAT HW2 1  0.4238  1.008

#ifndef FLEXIBLE
[ settles ]
1 1 0.1 0.1633

[ exclusions ]
1 2 3
2 1 3
3 1 2
#else
[ bonds ]
1 2 1 0.1 345000.0
1 3 1 0.1 345000.0

[ angles ]
2 1 3 1 109.47 383.0
#endif

#ifdef POSRES
#include "%s"
#endif
''' % WATER_POSRES



# FIXME: review data structure
class TopContainer(object):
    """Topology container."""

    class MoleculeType(object):
        """[ moleculetype ]"""
        __slots__ = ['molname', 'atoms', 'bonds', 'pairs', 'angles',
                     'propers', 'impropers']

        def __init__(self, molname):
            self.molname = molname
            self.atoms = []
            self.bonds = []
            self.pairs = set()
            self.angles = []
            self.propers = []
            self.impropers = []

    def __init__(self):
        self.atomtypes = {}
        self.moleculetype = []
        self.specials = ""
        self.length = 0


    def pushMolType(self, molname):
        self.moleculetype.append(self.MoleculeType(molname) )

        return self.moleculetype[-1]


    def __len__(self):
        if not self.length:
            for comp in self.moleculetype:
                self.length += (len(comp.molname) + len(comp.atoms) +
                                len(comp.bonds) + len(comp.angles) +
                                len(comp.propers) + len(comp.impropers) )
        return self.length

    
class GromacsTop(object):
    """Basic Gromacs topology writer."""

    def __init__(self):
        self.coords = []
        self.box_dims = None
        self.top = TopContainer()
        self.parmtop = ""
        self.inpcrd = ""

        self.moltypes = []
        self.tot_natoms = 0
        self.nwat = 0

        self.molidx = {}

        self.mol_numbers = None
        self.mols = None


    def readParm(self, parmtop, inpcrd):
        """
        Extract topology and coordinate information from AMBER parmtop and
        inpcrd.  Store data internally.

        :param parmtop: parmtop file name
        :type parmtop: string
        :param inpcrd: inpcrd file name
        :type inprcd: string
        :raises: SetupError
        """

        amber = Sire.IO.Amber()

        try:
            # (Sire.Mol.Molecules,  Sire.Vol.PeriodicBox or Sire.Vol.Cartesian)
            mols, perbox = amber.readCrdTop(inpcrd, parmtop)
        except UserWarning as error:
            raise errors.SetupError('error opening %s/%s' % (parmtop, inpcrd) )

        self.parmtop = parmtop
        self.inpcrd = inpcrd

        resnames = OrderedDict()
        i = 0

        mol_numbers = mols.molNums()
        mol_numbers.sort()

        # first pass to get total number of atoms
        self.tot_natoms = sum(mols.at(num).molecule().nAtoms()
                              for num in mol_numbers)

        mcnt = 0

        # second pass to get atomtypes: grompp allows only one such section
        for num in mol_numbers:
            mol = mols.at(num).molecule()
            natoms = mol.nAtoms()
            key = tuple([str(res.name().value()) for res in mol.residues()])

            # store unique molecules because only topological data is needed
            try:
                resnames[key][1] += 1
            except KeyError:
                resnames[key] = [num, 1, [], natoms]

            idx = []

            if key[0] != 'WAT':
                if len(key) == 1:
                    mol_name = '%s' % key[0]
                else:
                    # FIXME: analyse names e.g. if protein, etc.
                    mcnt += 1
                    mol_name = 'MOL%i' % mcnt

                self.moltypes.append( (mol_name, 1) )

            for atom in mol.atoms():
                ambertype = str(atom.property('ambertype') )

                # silly Gromacs doesn't get along with type starting with digit
                if ambertype[0].isdigit():
                    ambertype = ATOM_PREFIX + ambertype

                mass = atom.property('mass').value()
                lj = atom.property('LJ')
                coords = atom.property('coordinates')

                atom_name = str( atom.name().value() )
                resname = str(atom.residue().name().value() )
                resnum = atom.residue().number().value() % 99999

                if resname == 'WAT':
                    atom_name = water_atom_names[atom_name]

                self.coords.append( (resnum, resname, atom_name,
                                     coords[0] * const.A2NM,
                                     coords[1] * const.A2NM,
                                     coords[2] * const.A2NM) )

                try:
                    nprot = atom.property('element').nProtons()
                except UserWarning:
                    nprot = 0

                # FIXME: check if duplicates are really the same?
                self.top.atomtypes[ambertype] = ( (mass,
                                          lj.sigma().value() * const.A2NM,
                                          lj.epsilon().value() * const.CAL2J) )

                idx.append(i)
                i += 1

            resnames[key][2].extend(idx)

        # FIXME: only orthorombic box
        try:
            dims = perbox.dimensions() # Sire.Maths.Vector

            x = dims.x()
            y = dims.y()
            z = dims.z()
        except:                         # FIXME: which exception?
            x, y, z = 0.0, 0.0, 0.0

        self.box_dims = x * const.A2NM, y * const.A2NM, z * const.A2NM


        mols_tmp = []
        mcnt = 0

        # FIXME: we need to preserve molecule order, e.g. ions may be
        #        alternating
        for seq, num in resnames.iteritems():
            slen = len(seq)

            if seq[0] == 'WAT' and slen == 1:
                mol_name = 'WAT'
                water_mol_number, self.nwat, d, d = num
            else:
                if slen == 1:
                    mol_name = '%s' % seq[0]
                else:
                    # FIXME: analyse names e.g. if protein, etc.
                    mcnt += 1
                    mol_name = 'MOL%i' % mcnt

                mols_tmp.append( (num[0], mol_name, num[1]) )

            self.molidx[mol_name] = (num[2], num[3])

        # third pass only over unique molecules (molecule types)
        for num, mol_name, mol_cnt in mols_tmp:
            mol = mols.at(num).molecule()
            is_atom = False

            try:
                params = mol.property('amberparameters') # Sire.Mol.AmberParameters
            except UserWarning:
                is_atom = True

            moltype = self.top.pushMolType(mol_name)

            for atom in mol.atoms():
                ambertype = str(atom.property('ambertype') )

                # FIXME: combine with loop above
                if ambertype[0].isdigit():
                    ambertype = ATOM_PREFIX + ambertype

                charge = atom.property('charge').value()
                mass = atom.property('mass').value()

                atom_name = str(atom.name().value() )
                resname = str(atom.residue().name().value() )

                moltype.atoms.append( (ambertype, resname, atom_name, charge,
                                       mass) )

            if is_atom:
                continue

            try:
                mol.property('bond')
            except UserWarning:
                continue

            for bond in params.getAllBonds():  # Sire.Mol.BondID
                at0 = bond.atom0()  # Sire.Mol.AtomIdx!
                at1 = bond.atom1()

                k, r = params.getParams(bond)

                # x100 because kJ/nm
                moltype.bonds.append( (at0.value() + 1, at1.value() + 1,
                                       r * const.A2NM, 200.0 * k * const.CAL2J) )

            moltype.bonds.sort()


            try:
                mol.property('angle')
            except UserWarning:
                continue

            for angle in params.getAllAngles():  # Sire.Mol.AngleID
                at0 = angle.atom0()  # Sire.Mol.AtomIdx!
                at1 = angle.atom1()
                at2 = angle.atom2()

                k, theta = params.getParams(angle)
                moltype.angles.append( (at0.value() + 1, at1.value() + 1,
                                        at2.value() + 1, theta * const.RAD2DEG,
                                        2.0 * k * const.CAL2J) )

            moltype.angles.sort()


            try:
                mol.property('dihedral')
            except UserWarning:
                continue

            intrascale = mol.property('intrascale')

            for dihedral in params.getAllDihedrals():  # Sire.Mol.DihedralID
                at0 = dihedral.atom0()  # Sire.Mol.AtomIdx!
                at1 = dihedral.atom1()
                at2 = dihedral.atom2()
                at3 = dihedral.atom3()

                at0idx = at0.value() + 1
                at3idx = at3.value() + 1

                p = params.getParams(dihedral)
                x = []

                n = 3
                for i in range(0, len(p), n):       # k, np, phase
                    x.append(p[i:i+n])

                moltype.propers.append( (at0idx, at1.value() + 1,
                                         at2.value() + 1, at3idx, x) )

                sf = intrascale.get(at0, at3)

                if not (sf.lj() == sf.coulomb() == 0.0):
                    if not (at3idx, at0idx) in moltype.pairs:
                        moltype.pairs.add( (at0idx, at3idx) )

            moltype.propers.sort()


            try:
                mol.property('improper')
            except UserWarning:
                continue

            for dihedral in params.getAllImpropers():
                at0 = dihedral.atom0()
                at1 = dihedral.atom1()
                at2 = dihedral.atom2()
                at3 = dihedral.atom3()

                pk, pn, phase = params.getParams(dihedral)

                moltype.impropers.append( (at0.value() + 1, at1.value() + 1,
                                           at2.value() + 1, at3.value() + 1,
                                           phase * const.RAD2DEG,
                                           pk * const.CAL2J, pn) )

            moltype.impropers.sort()


        if self.nwat:
            wmol = mols.at(water_mol_number).molecule()
            parms = wmol.property('amberparameters')

            # FIXME: rather weak tests
            if wmol.nAtoms() == 4:
                self.top.specials = TIP4PEW_header
            elif parms.getParams(parms.getAllBonds()[0])[0] == 553:
                self.top.specials = TIP3P_header
            else:
                self.top.specials = SPCE_header

        self.mol_numbers = mol_numbers
        self.mols = mols


    def addAtomTypes(self, atomtypes):
        """Add atom types.
        :param atomtypes: atom type list with atom type, mass, sigma, epsilon
        :type atomtypes: list
        """        

        for ambertype, mass, sigma, epsilon in atomtypes:
            self.top.atomtypes[ambertype] = ( (mass, sigma, epsilon) )

    def writeGro(self, filename):
        """Write gro coordinate file."""

        atnum = 0

        with open(filename, 'w') as gro:
            gro.write('FESetup version X\n%i\n' % self.tot_natoms)

            for line in self.coords:
                atnum += 1
                gro.write('%5d%-5.5s%5.5s%5d %13.8f %13.8f %13.8f\n' %
                          (line[:3] + (atnum,) + line[3:]) )
                atnum %= 99999

            gro.write('%.7f %.7f %.7f\n' % self.box_dims)


    # FIXME: #include 'atomtypes.itp'
    def writeTop(self, topname, typename='atomtypes.atp', pertname='',
                 itp=True, itp_inc_file=const.GROMACS_PERT_ITP):
        """Write top or itp file.
        :param topname: AMBER parmtop file name
        :type topname: string
        :param typename: filename for atom types, if empty include in TOP/ITP
        :type typename: string
        :param pertname: perturbed ligand name
        :type pertname: string
        :param itp: write ITP or TOP
        :type itp: bool
        :param itp_inc_file: ITP file to include
        :type itp_inc_file: string
        """

        with open(topname, 'w') as top:
            if not itp:
                top.write(default_header)

            top.write('\n[ atomtypes ]\n;name btype   mass   charge '
                      'ptype      sigma      epsilon\n')

            if not typename:
                handle = top
            else:
                handle = open(typename, 'w')

            for typ in sorted(self.top.atomtypes):
                handle.write('%-4s %-4s %8.3f   0.0000  A  %12.6e %12.6e\n' %
                          ( (typ, typ) + self.top.atomtypes[typ][:]) )

            if typename:
                handle.close()

            first = False

            for mt in self.top.moleculetype:
                if pertname and not first:  # FIXME: ligand first mol
                    if typename:
                        top.write('#include "%s"\n' % typename)

                    top.write('\n#include "%s"\n' % itp_inc_file)
                    mt.molname = pertname
                    first = True
                    continue

                top.write('\n[ moleculetype ]\n%s 3\n' % mt.molname)

                natom = 0
                cgnr = 0

                top.write('\n[ atoms ]\n;   nr  type resno resnm '
                          'atom    cgnr      charge        mass      '
                          'typeB    chargeB      massB\n')

                for atom in mt.atoms:
                    natom +=1
                    cgnr += 1

                    top.write('%7i %3s %6i %5s %4s %7i %11.6f %11.4f\n' %
                              ( (natom, atom[0], 1) + atom[1:3] + (cgnr,) +
                                atom[3:]) )

                if mt.bonds:
                    top.write('\n[ bonds ]\n;  ai        aj f           c0'
                              '           c1\n')

                for bond in mt.bonds:
                    top.write('%7i %7i 1 %12.6e %12.6e\n' % bond)


                if mt.angles:
                    top.write('\n[ angles ]\n;  ai        aj      ak f'
                              '           c0           c1\n')

                for angle in mt.angles:
                    top.write('%7i %7i %7i 1 %12.6e %12.6e\n' % angle)


                if mt.pairs:
                    top.write('\n[ pairs ]\n;    ai      aj f\n')

                for pair in sorted(mt.pairs):
                    top.write('%7i %7i 1\n' % pair)


                if mt.propers:
                    top.write('\n[ dihedrals ] ;propers\n;    i        j       '
                              'k       l f     phase     pk       pn\n')

                for proper in mt.propers:
                    for dih in proper[4]:
                        top.write('%7i %7i %7i %7i 9 %11.2f %11.5f %i\n' %
                                  (proper[:4] + (dih[2] * const.RAD2DEG,
                                                 dih[0] * const.CAL2J) +
                                   (dih[1],) ) )

                if mt.impropers:
                    top.write('\n[ dihedrals ] ;impropers\n;    i        j    '
                              '   k       l f      phase     pk       pn\n')

                for improper in mt.impropers:
                    top.write('%7i %7i %7i %7i 4 %11.2f %11.5f %i\n' % improper)

                posres_file = (const.GROMACS_POSRES_PREFIX + mt.molname +
                               const.GROMACS_ITP_EXT)

                top.write('\n#ifdef POSRES\n#include "%s"\n#endif\n\n' %
                          posres_file)

            if self.nwat:
                top.write(self.top.specials)

            # FIXME: system name
            if not itp:
                top.write('\n[ system ]\n%s\n\n[ molecules ]\n; Compound        '
                   'nmols\n' % 'SysX')

                for name, cnt in self.moltypes:
                    top.write('%s %i\n' % (name, cnt) )

                if self.nwat:
                    top.write('WAT %i\n' % self.nwat)

            top.write('\n\n')


    def unwrap(self, coords, boxx, boxy, boxz):
        """
        Unwrap coordinates because Gromacs uses atom-based wrapping.

        The current code assumes largest molecule is not larger than the
        box minus an arbitrary buffer.  Possible alternative: check all
        bonds an move atom if bond length is larger than half the box
        length.

        :param coords: flat list of coordinates
        :type coords: list
        :param boxx: box dimension in x direction
        :type boxx: float
        :param boxy: box dimension in y direction
        :type boxy: float
        :param boxz: box dimension in z direction
        :type boxz: float
        """

        offset = 0

        # FIXME: arbitrary buffer, unwrapping may fail?
        maxx = (boxx - BOX_BUFFER)**2
        maxy = (boxy - BOX_BUFFER)**2
        maxz = (boxz - BOX_BUFFER)**2

        boxx2 = boxx / 2.0
        boxy2 = boxy / 2.0
        boxz2 = boxz / 2.0

        for num in self.mol_numbers:
            mol = self.mols.at(num).molecule()
            natoms = mol.nAtoms()

            first = True

            for residue in mol.residues():
                for atom in residue.atoms():
                    idx = 3 * (atom.index().value() + offset)

                    if first:
                        x0 = coords[idx]
                        y0 = coords[idx+1]
                        z0 = coords[idx+2]

                        if x0 > boxx:
                            coords[idx] -= math.copysign(boxx, x1)

                        if y0 > boxy:
                            coords[idx+1] -= math.copysign(boxy, y1)

                        if z0 > boxz:
                            coords[idx+2] -= math.copysign(boxz, z1)

                        first = False
                        continue

                    x1 = coords[idx]
                    y1 = coords[idx+1]
                    z1 = coords[idx+2]

                    xd = (x0 - x1)**2
                    yd = (y0 - y1)**2
                    zd = (z0 - z1)**2

                    # FIXME: orthorombic only
                    if xd > maxx:
                        if x1 > boxx2:
                            coords[idx] -= boxx
                        else:
                            coords[idx] += boxx

                    if yd > maxy:
                        if y1 > boxy2:
                            coords[idx+1] -= boxy
                        else:
                            coords[idx+1] += boxy

                    if zd > maxz:
                        if z1 > boxz2:
                            coords[idx+2] -= boxz
                        else:
                            coords[idx+2] += boxz

            offset += natoms


    def __len__(self):
        return len(self.top)


def mixer(top0, top1, filename = const.GROMACS_PERT_ITP,
          typename = const.GROMACS_PERT_ATP):
    """
    Creates a perturbed Gromacs topology file from two top file mix-ins.

    :param parmtop: Gromacs topology object for state0
    :type parmtop: GromacsTop
    :param inpcrd: Gromacs topology object for state1
    :type inprcd: GromacsTop
    :param filename: itp file name for writing
    :type filename: string
    :raises: SetupError
    """

    import copy

    # FIXME: is this check sufficient?
    if len(top0) != len(top1):
        raise errors.SetupError('topologies of different lengths: %s %s' %
                                (top0.parmtop, top1.parmtop) )

    with open(typename, 'w') as atype:
        #atype.write('\n[ atomtypes ]\n;name btype   mass   charge '
        #            'ptype      sigma      epsilon\n')

        combined = copy.copy(top0.top.atomtypes)
        combined.update(top1.top.atomtypes)

        for typ in sorted(combined):
            atype.write('%-4s %-4s %8.3f   0.0000  A  %12.6e %12.6e\n' %
                        ( (typ, typ) + combined[typ][:]) )


    with open(filename, 'w') as itp:
        for mt0, mt1 in zip(top0.top.moleculetype, top1.top.moleculetype):
            # FIXME: sensible molecule name
            itp.write('\n[ moleculetype ]\n%s 3\n' % mt0.molname)

            natom = 0
            cgnr = 0

            itp.write('\n[ atoms ]\n;   nr  type resno resnm '
                      'atom    cgnr      charge        mass'
                      ' typeB   chargeB       massB\n')

            for atom0, atom1 in zip(mt0.atoms, mt1.atoms):
                natom +=1
                cgnr += 1

                itp.write('%7i %3s %6i %5s %4s %7i %11.6f %11.4f %3s '
                          '%11.6f %11.4f\n' %
                          ( (natom, atom0[0], 1) + atom0[1:3] + (cgnr,) +
                            atom0[3:] + (atom1[0],) + atom1[3:]) )


            if mt0.bonds:
                itp.write('\n[ bonds ]\n;  ai        aj f           r0'
                          '            k            r0B           kB\n')

            for bond0, bond1 in zip(mt0.bonds, mt1.bonds):
                if bond0[2] == 0.0 or bond1[2] == 0.0:
                    logger.write('Warning: zero bonds in %s, %s' %
                                 (' '.join(str(i) for i in bond0),
                                  ' '.join(str(i) for i in bond1)))

                itp.write('%7i %7i 1 %12.6e %12.6e   %12.6e %12.6e\n' %
                          (bond0 + bond1[2:]) )


            if mt0.angles:
                itp.write('\n[ angles ]\n;  ai        aj      ak f'
                          '           a0            k'
                          '            a0B           kB\n')

            for angle0, angle1 in zip(mt0.angles, mt1.angles):
                if angle0[3] == 0.0 or angle1[3] == 0.0:
                    logger.write('Warning: zero angles in %s, %s' %
                                 (' '.join(str(i) for i in angle0),
                                  ' '.join(str(i) for i in angle1)))

                itp.write('%7i %7i %7i 1 %12.6e %12.6e   %12.6e %12.6e\n' %
                          (angle0 + angle1[3:]) )


            if mt0.pairs:
                itp.write('\n[ pairs ]\n;    ai      aj f\n')

            if len(mt0.pairs) != len(mt1.pairs):
                raise errors.SetupError('pairs of different size')

            for pair0, pair1 in zip(sorted(mt0.pairs), sorted(mt1.pairs) ):
                if pair0 != pair1:
                    raise errors.SetupError('different pairs')
                    
                itp.write('%7i %7i 1\n' % pair0)


            if mt0.propers:
                itp.write('\n[ dihedrals ] ;propers\n;    i        j       '
                          'k       l f       phase          pk pn'
                          '        phase          pk pn\n')

            for proper0, proper1 in zip(mt0.propers, mt1.propers):
                for dih0, dih1 in zip(proper0[4], proper1[4]):
                    # for some reason, leap creates null entries
                    if dih0[0] == 0.0 and dih1[0] == 0.0:
                        continue

                    itp.write('%7i %7i %7i %7i 9 %11.2f %11.5f %i   '
                              '%11.2f %11.5f %i\n' %
                              (proper0[:4] + (dih0[2] * const.RAD2DEG,
                                              dih0[0] * const.CAL2J) +
                               (dih0[1],) + (dih1[2] * const.RAD2DEG,
                                             dih1[0] * const.CAL2J) +
                               (dih1[1],)) )


            if mt0.impropers:
                itp.write('\n[ dihedrals ] ;impropers\n;    i        j    '
                          '   k       l f       phase          pk pn'
                          '        phase          pk pn\n')

            for improper0, improper1 in zip(mt0.impropers, mt1.impropers):
                        itp.write('%7i %7i %7i %7i 4 %11.2f %11.5f %i   '
                                  '%11.2f %11.5f %i\n' %
                                  (improper0 + improper1[4:]) )

        itp.write('\n\n')



if __name__ == '__main__':
    import sys

    nargs = len(sys.argv)

    if nargs == 3:
        top = GromacsTop()
        top.readParm(sys.argv[1], sys.argv[2])
        top.writeTop('test.top', '', '', False)
        top.writeGro('test.gro')
    elif nargs == 5:
        top0 = GromacsTop()
        top0.readParm(sys.argv[1], sys.argv[2])
        top0.writeTop('0.itp')

        top1 = GromacsTop()
        top1.readParm(sys.argv[3], sys.argv[4])
        top1.writeTop('1.itp')

        mixer(top0, top1)
    else:
        sys.exit('Usage: %s parmtop1 inpcrd1 [parmtop2 inpcrd2]' % sys.argv[0])
