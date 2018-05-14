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
Generates a DL_POLY FIELD and CONFIG file from an AMBER parmtop and inpcrd.
'''

__revision__ = "$Id$"



import os, sys, math
from operator import itemgetter

import Sire.IO
import Sire.MM
import Sire.Maths

from FESetup import const, errors, logger


BOX_BUFFER = 3.0



class DLPolyField(object):
    """Basic DL_POLY topology writer."""

    def __init__(self):
        self.atoms = []
        self.posres = []
        self.bonds = []
        self.constraints = []
        self.angles = []
        self.rigids = []
        self.propers = []
        self.impropers = []
        self.vdw = {}
        
        self.coords = []

        self.mol_numbers = None
        self.mols = None

        self.box_dims = None

        self.parmtop = ""
        self.inpcrd = ""

        self.perbox = None

        
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
            mols, self.perbox = amber.readCrdTop(inpcrd, parmtop)
        except UserWarning as error:
            raise errors.SetupError('error opening %s/%s' % (parmtop, inpcrd) )

        self.parmtop = parmtop
        self.inpcrd = inpcrd

        mol_numbers = mols.molNums()
        mol_numbers.sort()

        offset = 1
        atomtypes = {}
        rigids = []
        
        for num in mol_numbers:
            mol = mols.at(num).molecule()
            natoms = mol.nAtoms()

            res = mol.residues()[0]

            # FIXME: always named WAT?
            if str(res.name().value() ) == 'WAT':
                rigids.append([atom.index().value() + offset
                               for atom in res.atoms() ] )

            for atom in mol.atoms():
                ambertype = str(atom.property('ambertype') )

                sfx = ''

                for ch in ambertype:
                    if ch.istitle():
                        sfx += 'U'
                    else:
                        sfx += 'l'

                charge = atom.property('charge').value()
                mass = atom.property('mass').value()
                coords = atom.property('coordinates')

                lj = atom.property('LJ')
                sigma = lj.sigma().value() * const.RSTAR_CONV
                epsilon = lj.epsilon().value()

                resname = str(atom.residue().name().value() )
                resnum = atom.residue().number().value()

                element = atom.property('element').symbol()

                # FIXME: really TIP4?, residue always named WAT?
                if ambertype == 'EP' and resname == 'WAT':
                    element = 'EP'       # for CONFIG comment

                atype = ambertype + '_' + sfx
                self.coords.append( (atype, resnum, resname, element,
                                     coords[0], coords[1], coords[2]) )

                self.atoms.append( (atype, mass, charge, resnum, resname) )

                atomtypes[atype] = (sigma, epsilon)


            try:
                params = mol.property('amberparameters') # Sire.Mol.AmberParameters
            except UserWarning:
                offset += natoms
                continue

            try:
                mol.property('bond')
            except UserWarning:
                offset += natoms
                continue

            for bond in params.getAllBonds():  # Sire.Mol.BondID
                at0 = bond.atom0()  # Sire.Mol.AtomIdx!
                at1 = bond.atom1()
                k, r = params.getParams(bond)

                at0sel = mol.select(at0)
                resn = str(at0sel.residue().name().value() )
                elem0 = at0sel.property('element').symbol()
                elem1 = mol.select(at1).property('element').symbol()

                idx0 = at0.value() + offset
                idx1 = at1.value() + offset

                self.bonds.append( (idx0, idx1, 2.0 * k, r) )

                # FIXME: make all-H vs water-only-H an option?; EP in TIP4?;
                #        consider rigid-body for TIP3P et al.
                if elem0 == 'H' or elem1 == 'H':
                    if resn == 'WAT':
                        rflag = 1
                    else:
                        rflag = 0

                    dist = Sire.Maths.Vector.distance(
                        mol.atom(at0).property("coordinates"),
                        mol.atom(at1).property("coordinates") )

                    self.constraints.append( (idx0, idx1, dist, rflag) )

            self.bonds.sort()


            try:
                mol.property('angle')
            except UserWarning:
                offset += natoms
                continue

            for angle in params.getAllAngles():  # Sire.Mol.AngleID
                at0 = angle.atom0()  # Sire.Mol.AtomIdx!
                at1 = angle.atom1()
                at2 = angle.atom2()
                k, theta = params.getParams(angle)

                self.angles.append( (at0.value() + offset, at1.value() + offset,
                                     at2.value() + offset, 2.0 * k,
                                     theta * const.RAD2DEG) )

            self.angles.sort()


            try:
                mol.property('dihedral')
            except UserWarning:
                offset += natoms
                continue

            intrascale = mol.property('intrascale')

            pairs = set()

            for dihedral in params.getAllDihedrals():  # Sire.Mol.DihedralID
                at0 = dihedral.atom0()  # Sire.Mol.AtomIdx!
                at1 = dihedral.atom1()
                at2 = dihedral.atom2()
                at3 = dihedral.atom3()

                idx0 = at0.value() + offset
                idx3 = at3.value() + offset

                sf = intrascale.get(at0, at3)

                # work-around for pairs double-counting bug in
                # Sire.IO.Amber().readCrdTop()
                if (idx0, idx3) in pairs or (idx3, idx0) in pairs:
                    scnb, scee = 0.0, 0.0
                else:
                    scee = sf.lj()
                    scnb = sf.coulomb()
                    pairs.add( (idx0, idx3) )
                    pairs.add( (idx3, idx0) )

                p = params.getParams(dihedral)

                for i in range(0, len(p), 3):
                    pk = p[i]
                    pn = p[i+1]
                    phase = p[i+2]

                    self.propers.append( (idx0, at1.value() + offset,
                                          at2.value() + offset, idx3,
                                          pk, phase * const.RAD2DEG, pn,
                                          scnb, scee) )

                    scnb, scee = 0.0, 0.0  # count multi-terms only once

            self.propers.sort(key = itemgetter(0, 1, 2, 3) )

            try:
                mol.property('improper')
            except UserWarning:
                offset += natoms
                continue

            for dihedral in params.getAllImpropers():
                at0 = dihedral.atom0()
                at1 = dihedral.atom1()
                at2 = dihedral.atom2()
                at3 = dihedral.atom3()

                pk, pn, phase = params.getParams(dihedral)

                self.impropers.append( (at0.value() + offset,
                                        at1.value() + offset,
                                        at2.value() + offset,
                                        at3.value() + offset,
                                        pk, phase * const.RAD2DEG, pn,
                                        0.0, 0.0) )

            self.impropers.sort()


            offset += natoms


        self.rigids = rigids
        self.mol_numbers = mol_numbers
        self.mols = mols

        div = math.pow(2.0, -1.0 / 6.0)

        for key0, a0 in atomtypes.iteritems():
            for key1, a1 in atomtypes.iteritems():
                vkey = (key0, key1)

                if vkey in self.vdw or (key1, key0) in self.vdw:
                    continue

                prod = a0[1] * a1[1]
                
                if prod == 0.0:         # DL_POLY does not check for zero!
                    continue

                # AMBER vdW mixing rules
                sigma = (a0[0] + a1[0]) * div
                epsilon = math.sqrt(prod)

                fac = math.pow(sigma, 6.0)
                B = 4.0 * epsilon * fac
                A = B * fac

                self.vdw[vkey] = (A, B)
                                              

    def writeConfig(self, configname = 'CONFIG', center = True):
        """
        Write DL_POLY CONFIG coordinate file.

        Format:
        =======

        record 1:
          header(A72)
        record 2:
          levcfg (integer): 0 = x, 1 = xv, 2 = xvf;
          imcon (integer): 0 = no PBC, 1 = cubic, 2 = orthorhombic,
                           3 = triclinic, 4 = truncated octahedron,
                           5 = rhombic, 6 = xz parallelogram/no z,
                           7 = hexagonal prism;
          megatm (integer): optional, total number of particles;
        record 3:
          cell(1-3: 3 floats)
        record 4:
          cell(4-6: 3 floats)
        record 5:
          cell(7-9: 3 floats)
        coordinates: (origin in centre of box)
          record I:
            atom name(A8); index
          record II (levcfg = 0):
            x, y, z
          record III (levcfg = 1):
            vx, vy, vz
          record IV (levcfg = 2):
            fx, fy, fz

        :param configname: name of the CONFIG file
        :type configname: string
        :param center: shift the box by half the box length = expect rst7
                       coordinate system with origin in box corner
        :type center: bool
         """

        atnum = 0

        try:
            dims = self.perbox.dimensions() # Sire.Maths.Vector

            x = dims.x()
            y = dims.y()
            z = dims.z()

            # FIXME: check for triclinic
            if x == y == z:
                imcon = 1
            else:
                imcon = 2
        except:
            imcon = 0


        delx = dely = delz = 0.0

        # FIXME: triclinic
        if center:
            minc = [ sys.float_info[0],  sys.float_info[0],  sys.float_info[0] ]
            maxc = [-sys.float_info[0], -sys.float_info[0], -sys.float_info[0] ]

            for coord in self.coords:
                if coord[4] < minc[0]:
                    minc[0] = coord[4]
                if coord[5] < minc[1]:
                    minc[1] = coord[5]
                if coord[6] < minc[2]:
                    minc[2] = coord[6]

                if coord[4]   > maxc[0]:
                    maxc[0] = coord[4]
                if coord[5] > maxc[1]:
                    maxc[1] = coord[5]
                if coord[6] > maxc[2]:
                    maxc[2] = coord[6]

            delx = (maxc[0] - minc[0]) / 2.0
            dely = (maxc[1] - minc[1]) / 2.0
            delz = (maxc[2] - minc[2]) / 2.0


        # FIXME: triclinic
        with open(configname, 'w') as cnf:
            cnf.write('Created by FESetup\n%10i%10i%10i\n' %
                      (0, imcon, len(self.coords) ) )

            if imcon:
                cnf.write('%20.12f%20.12f%20.12f\n' % (x, 0.0, 0.0) )
                cnf.write('%20.12f%20.12f%20.12f\n' % (0.0, y, 0.0) )
                cnf.write('%20.12f%20.12f%20.12f\n' % (0.0, 0.0, z) )

            for i, coord in enumerate(self.coords):
                cnf.write('%-10s %7i %i%s %s\n%20.8f%20.8f%20.8f\n' %
                          (coord[0], i+1, coord[1], coord[2], coord[3],
                           coord[4] - delx, coord[5] - dely, coord[6] - delz) )


    def writeField(self, topname, do_rigid = True):
        """
        Write DL_POLY FIELD file.

        :param topname: FIELD file name
        :type topname: string
        :param do_rigid: write out rigid sections
        :type do_rigid: bool
        """

        with open(topname, 'w') as top:
            # system is written out as one "molecule"
            top.write('Created by FESetup\nUnits kcal/mol\nMolecular types 1\n'
                      'Molecule name system\nnummols %i\n' % 1)

            top.write('atoms %i\n' % len(self.atoms) )

            for atom in self.atoms:
                # 1 is repeat counter (optional) for frozen atoms (=0 here)
                # residue number+residue name is ignored ("comment")
                top.write('%-10s %10.5f %12.5f %4d %4d %d%s\n' %
                          (atom[0], atom[1], atom[2], 1, 0, atom[3], atom[4]) )

            if self.posres:
                top.write('teth %i\n' % len(self.posres) )

                for restr in self.posres:
                    top.write('harm %6i %12.2f\n' % (restr[0], restr[1]) )

            top.write('bonds %i\n' % len(self.bonds) )

            for bond in self.bonds:
                top.write('harm %6i %6i %12.2f %8.5f\n' %
                          (bond[0], bond[1], bond[2], bond[3]) )

            if self.angles:
                top.write('angles %i\n' % len(self.angles) )

            for angle in self.angles:
                top.write('harm %6i %6i %6i %12.2f %8.2f\n' %
                          (angle[0], angle[1], angle[2], angle[3], angle[4]) )

            if self.constraints:
                if do_rigid:
                    nrigid = sum(c[3] for c in self.constraints)
                else:
                    nrigid = 0

                ncons = len(self.constraints) - nrigid

                if ncons > 0:
                    top.write('constraints %i\n' % ncons)

                    for cons in self.constraints:
                        if do_rigid and cons[3]:
                            continue

                        top.write('%i %i %.5f\n' % (cons[0], cons[1], cons[2]) )

            if do_rigid and self.rigids:
                top.write('rigid %i\n' % len(self.rigids) )

                for rigid in self.rigids:
                    top.write('%i' % len(rigid) )

                    for idx in rigid:
                        top.write(' %i' % idx)

                    top.write('\n')

            # non-bonded sections must always be written
            top.write('dihedrals %i\n' % (len(self.propers) +
                                          len(self.impropers) ) )

            for dihedral in self.propers + self.impropers:
                top.write('cos %6i %6i %6i %6i %8.3f %10.3f %5.2f '
                          '%7.5f %7.5f\n' %
                          (dihedral[0], dihedral[1], dihedral[2], dihedral[3],
                           dihedral[4], dihedral[5], dihedral[6],
                           dihedral[7], dihedral[8]) )

            top.write('finish\n')

            top.write('vdw %i\n' % len(self.vdw) )

            for v, p in self.vdw.iteritems():
                top.write('%-8s %-8s 12-6 %e %e\n' % (v[0], v[1], p[0], p[1]) )

            top.write('close\n')


    def unwrap(self, coords, boxx, boxy, boxz):
        """
        Unwrap coordinates because DL_POLY uses atom-based wrapping.

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
                        coords[idx] -= math.copysign(boxx, x1)

                    if yd > maxy:
                        coords[idx+1] -= math.copysign(boxy, y1)

                    if zd > maxz:
                        coords[idx+2] -= math.copysign(boxz, z1)

            offset += natoms


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 3:
        top = DLPolyField()
        top.readParm(sys.argv[1], sys.argv[2])
        top.writeField('FIELD')
        top.writeConfig('CONFIG')
    else:
        sys.exit('Usage: %s parmtop inpcrd' % sys.argv[0])
