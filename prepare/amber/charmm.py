#  Copyright (C) 2015  Hannes H Loeffler
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
Generates CHARMM prm, psf and cord files from AMBER parmtop and inpcrd.
Only the extended format is supported but CMAP and CHEQ are not.
'''

__revision__ = "$Id: gromacs.py 395 2015-01-26 10:45:12Z halx $"



import os

import Sire.IO

from FESetup import const, errors, logger


class CharmmTop(object):
    """Basic CHARMM prm and psf writer."""

    def __init__(self, parmtop, inpcrd):
        amber = Sire.IO.Amber()

        try:
            # (Sire.Mol.Molecules,  Sire.Vol.PeriodicBox or Sire.Vol.Cartesian)
            mols, perbox = amber.readCrdTop(inpcrd, parmtop)
        except UserWarning as error:
            raise errors.SetupError('error opening %s/%s' % (parmtop, inpcrd) )

        try:
            dims = perbox.dimensions() # Sire.Maths.Vector

            x = dims.x()
            y = dims.y()
            z = dims.z()
        except:                         # FIXME: which exception?
            x, y, z = 0.0, 0.0, 0.0

        self.box_dims = x, y, z

        self.mol_numbers = mols.molNums()
        self.mol_numbers.sort()

        self.tot_natoms = sum(mols.at(num).molecule().nAtoms()
                              for num in self.mol_numbers)

        self.mols = mols


    def writeCrd(self, filename):
        """Write crd coordinate file.
* Expanded format for more than 100000 atoms (upto 10**10) and with
upto 8 character PSF IDs. (versions c31a1 and later)
         title
         NATOM (I10)
         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
           I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10
        """

        fmt = '%10i%10i  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8s%20.10f\n'
        atomno = 0

        with open(filename, 'w') as crd:
            crd.write('* Created by FESetup\n*\n%10i  EXT\n' % self.tot_natoms)

            for num in self.mol_numbers:
                mol = self.mols.at(num).molecule()
                natoms = mol.nAtoms()

                for atom in mol.atoms():
                    atomno += 1
                    resno = atom.residue().number().value()
                    res = str(atom.residue().name().value() )
                    atom_type = str( atom.name().value() )
                    coords = atom.property('coordinates')
                    segid = 'X'  # FIXME
                    resid = str(resno)  # FIXME
                    weight = 0.0

                    crd.write(fmt % (atomno, resno, res, atom_type, coords[0],
                                     coords[1], coords[2], segid, resid,
                                     weight) )


    def writePrmPsf(self, prmname, psfname):
        """Write prm and psf files.

        ATOM                             (Flexible paramters only)
         MASS   code   type   mass       (Flexible paramters only)

        EQUIvalence                      (Flexible paramters only)
         group  atom [ repeat(atom) ]    (Flexible paramters only)

        BOND
         atom atom force_constant distance

        ANGLe or THETA
         atom atom atom force_constant theta_min UB_force_constant UB_rmin

        DIHE or PHI
         atom atom atom atom force_constant periodicity phase

        IMPRoper or IMPHI
         atom atom atom atom force_constant periodicity phase

        NBONd or NONB  [nonbond-defaults]
         atom* polarizability  e  vdW_radius -
              [1-4 polarizability  e  vdW_radius]

        NBFIX
         atom_i* atom_j*  emin rmin [ emin14 [ rmin14 ]]

        """

        with open(prmname, 'w') as prm:
            pass

        bonds = []
        angles = []
        dihedrals = []
        impropers = []
        groups = []

        offset = 0
        atomno = 0

        with open(psfname, 'w') as psf:
            psf.write('PSF EXT\n         1\n* Created by FESetup\n%10i\n' %
                      self.tot_natoms)

            afmt = '%10i %-8s %-8s %-8s %-8s %-6s %14.6g%14.6g%8i\n'

            for num in self.mol_numbers:
                mol = self.mols.at(num).molecule()
                natoms = mol.nAtoms()

                for atom in mol.atoms():
                    atomno += 1
                    resno = atom.residue().number().value()
                    res = str(atom.residue().name().value() )
                    atom_type = str( atom.name().value() )
                    amber_type = str(atom.property('ambertype') )
                    segid = 'X'  # FIXME
                    resid = str(resno)  # FIXME
                    charge = atom.property('charge').value()
                    mass = atom.property('mass').value()

                    psf.write(afmt % (atomno, segid, resid, res, atom_type,
                                     amber_type, charge, mass, 0.0) )

                params = mol.property('amberparameters') # Sire.Mol.AmberParameters

                for bond in params.getAllBonds():  # Sire.Mol.BondID
                    at0 = bond.atom0()  # Sire.Mol.AtomIdx!
                    at1 = bond.atom1()
                    bonds.append( (at0.value() + 1 + offset,
                                   at1.value() + 1 + offset) )

                for angle in params.getAllAngles():  # Sire.Mol.AngleID
                    at0 = angle.atom0()  # Sire.Mol.AtomIdx!
                    at1 = angle.atom1()
                    at2 = angle.atom2()
                    angles.append( (at0.value() + 1 + offset,
                                    at1.value() + 1 + offset,
                                    at2.value() + 1 + offset) )

                for dihedral in params.getAllDihedrals(): # Sire.Mol.DihedralID
                    at0 = dihedral.atom0()  # Sire.Mol.AtomIdx!
                    at1 = dihedral.atom1()
                    at2 = dihedral.atom2()
                    at3 = dihedral.atom2()
                    dihedrals.append( (at0.value() + 1 + offset,
                                       at1.value() + 1 + offset,
                                       at2.value() + 1 + offset,
                                       at3.value() + 1 + offset) )

                for dihedral in params.getAllImpropers():
                    at0 = dihedral.atom0()
                    at1 = dihedral.atom1()
                    at2 = dihedral.atom2()
                    at3 = dihedral.atom2()
                    impropers.append( (at0.value() + 1 + offset,
                                       at1.value() + 1 + offset,
                                       at2.value() + 1 + offset,
                                       at3.value() + 1 + offset) )

                # groups: base pointer charge type (1=neutral,2=charged),
                #         entire group fixed?
                for residue in mol.residues():
                    charge = 0.0
                    first = True

                    for atom in residue.atoms():
                        if first:
                            gp_base = atom.index().value() + offset
                            first = False

                        charge += atom.property('charge').value()

                    if charge > 0.0:
                        gp_type = 2
                    else:
                        gp_type = 1

                    groups.append( (gp_base, gp_type) )

                offset += natoms


            #####
            psf.write('%10i bonds\n' % len(bonds) )

            for i, bond_pair in enumerate(bonds):
                psf.write('%10i%10i' % (bond_pair[0], bond_pair[1]) )

                if not (i + 1) % 4:
                    psf.write('\n')

            if (i + 1) % 4:
                psf.write('\n')

            psf.write('%10i angles\n' % len(angles) )
    
            for i, angle_triple in enumerate(angles):
                psf.write('%10i%10i%10i' % (angle_triple[0],
                                            angle_triple[1],
                                            angle_triple[2]) )

                if not (i + 1) % 3:
                    psf.write('\n')

            if (i + 1) % 3:
                psf.write('\n')

            psf.write('%10i dihedrals\n' % len(dihedrals) )
    
            for i, dihedral_quadruple in enumerate(dihedrals):
                psf.write('%10i%10i%10i%10i' % (dihedral_quadruple[0],
                                                dihedral_quadruple[1],
                                                dihedral_quadruple[2],
                                                dihedral_quadruple[3]) )

                if not (i + 1) % 2:
                    psf.write('\n')

            if (i + 1) % 2:
                psf.write('\n')

            psf.write('%10i impropers\n' % len(impropers) )
    
            for i, dihedral_quadruple in enumerate(impropers):
                psf.write('%10i%10i%10i%10i' % (dihedral_quadruple[0],
                                                dihedral_quadruple[1],
                                                dihedral_quadruple[2],
                                                dihedral_quadruple[3]) )

                if not (i + 1) % 2:
                    psf.write('\n')

            if (i + 1) % 2:
                psf.write('\n')

            psf.write('%10i donors\n%10i acceptors\n' % (0, 0) )
            psf.write('%10i non-bonded exclusions\n' % 0)

            # iblo (exclusion pointer): one entry for each atom
            for i in range(0, self.tot_natoms):
               psf.write('%10i' % 0)

               if not (i + 1) % 8:
                    psf.write('\n')

            if (i + 1) % 8:
                psf.write('\n')

            psf.write('%10i%10i group and ST2\n' % (len(groups), 0) )

            for i, group in enumerate(groups):
                psf.write('%10i%10i%10i' % (group[0], group[1], 0) )

                if not (i + 1) % 3:
                    psf.write('\n')

            if (i + 1) % 3:
                psf.write('\n')

            # molnt?
            psf.write('%10i molnt\n' % 1)

            for i in range(0, self.tot_natoms):
               psf.write('%10i' % 1)

               if not (i + 1) % 8:
                    psf.write('\n')

            if (i + 1) % 8:
                psf.write('\n')

            psf.write('%10i%10i lone pairs\n' % (0, 0) )



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

        pass



if __name__ == '__main__':
    import sys

    nargs = len(sys.argv)

    if nargs != 3:
        sys.exit('Usage: %s prmtop inpcrd' % sys.argv[0])


    top = CharmmTop(sys.argv[1], sys.argv[2])
    top.writeCrd('test.crd')
    top.writePrmPsf('test.prm', 'test.psf')
 
