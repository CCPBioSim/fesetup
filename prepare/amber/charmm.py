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



import os, math

import Sire.IO
import Sire.MM

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


    def writePrmPsf(self, rtfname, prmname, psfname):
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

        atoms = []
        bonds = []
        angles = []
        dihedrals = []
        impropers = []
        groups = []

        atom_params = {}
        bond_params = {}
        angle_params = {}
        dihedral_params = {}
        improper_params = {}
        nonbonded_params = []

        offset = 0
        atomno = 0

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
                lj = atom.property('LJ')

                if amber_type[0].islower():
                  amber_type = 'G' + amber_type.upper()

                atoms.append( (atomno, segid, resid, res, atom_type,
                               amber_type, charge, mass) )
                atom_params[amber_type] = (mass, lj)


            params = mol.property('amberparameters') # Sire.Mol.AmberParameters

            for bond in params.getAllBonds():  # Sire.Mol.BondID
                at0 = bond.atom0()  # Sire.Mol.AtomIdx!
                at1 = bond.atom1()
                name0 = str(mol.select(at0).property('ambertype') )
                name1 = str(mol.select(at1).property('ambertype') )
                k, r = params.getParams(bond)

                if name0[0].islower():
                  name0 = 'G' + name0.upper()

                if name1[0].islower():
                  name1 = 'G' + name1.upper()

                bonds.append( (at0.value() + 1 + offset,
                               at1.value() + 1 + offset) )
                bond_params[name0, name1] = (k, r)

            bonds.sort()

            for angle in params.getAllAngles():  # Sire.Mol.AngleID
                at0 = angle.atom0()  # Sire.Mol.AtomIdx!
                at1 = angle.atom1()
                at2 = angle.atom2()
                name0 = str(mol.select(at0).property('ambertype') )
                name1 = str(mol.select(at1).property('ambertype') )
                name2 = str(mol.select(at2).property('ambertype') )
                k, theta = params.getParams(angle)

                if name0[0].islower():
                  name0 = 'G' + name0.upper()

                if name1[0].islower():
                  name1 = 'G' + name1.upper()

                if name2[0].islower():
                  name2 = 'G' + name2.upper()

                angles.append( (at0.value() + 1 + offset,
                                at1.value() + 1 + offset,
                                at2.value() + 1 + offset) )
                angle_params[name0, name1, name2] = (k, theta * const.RAD2DEG)

            angles.sort()

            for dihedral in params.getAllDihedrals(): # Sire.Mol.DihedralID
                at0 = dihedral.atom0()  # Sire.Mol.AtomIdx!
                at1 = dihedral.atom1()
                at2 = dihedral.atom2()
                at3 = dihedral.atom3()

                name0 = str(mol.select(at0).property('ambertype') )
                name1 = str(mol.select(at1).property('ambertype') )
                name2 = str(mol.select(at2).property('ambertype') )
                name3 = str(mol.select(at3).property('ambertype') )

                if name0[0].islower():
                  name0 = 'G' + name0.upper()

                if name1[0].islower():
                  name1 = 'G' + name1.upper()

                if name2[0].islower():
                  name2 = 'G' + name2.upper()

                if name3[0].islower():
                  name3 = 'G' + name3.upper()

                p = params.getParams(dihedral)
                terms = []

                n = 3
                for i in range(0, len(p), n):       # k, np, phase
                    terms.append(p[i:i+n])

                #sf = intrascale.get(at0, at3)

                dihedrals.append( (at0.value() + 1 + offset,
                                   at1.value() + 1 + offset,
                                   at2.value() + 1 + offset,
                                   at3.value() + 1 + offset) )

                dihedral_params[name0, name1, name2, name3] = terms

            dihedrals.sort()

            for dihedral in params.getAllImpropers():
                at0 = dihedral.atom0()
                at1 = dihedral.atom1()
                at2 = dihedral.atom2()
                at3 = dihedral.atom3()

                name0 = str(mol.select(at0).property('ambertype') )
                name1 = str(mol.select(at1).property('ambertype') )
                name2 = str(mol.select(at2).property('ambertype') )
                name3 = str(mol.select(at3).property('ambertype') )

                if name0[0].islower():
                  name0 = 'G' + name0.upper()

                if name1[0].islower():
                  name1 = 'G' + name1.upper()

                if name2[0].islower():
                  name2 = 'G' + name2.upper()

                if name3[0].islower():
                  name3 = 'G' + name3.upper()

                term = params.getParams(dihedral)

                impropers.append( (at0.value() + 1 + offset,
                                   at1.value() + 1 + offset,
                                   at2.value() + 1 + offset,
                                   at3.value() + 1 + offset) )

                improper_params[name0, name1, name2, name3] = term

            impropers.sort()

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

                if charge > 0.01: # FIXME
                    gp_type = 2
                else:
                    gp_type = 1

                groups.append( (gp_base, gp_type) )

            offset += natoms


        #####
        with open(psfname, 'w') as psf:
            psf.write('PSF EXT\n\n'
                      '        1 !NTITLE\n'
                      '* Created by FESetup\n'
                      '\n%10i !NATOM\n' %
                      self.tot_natoms)

            # I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6
            afmt = '%10i %-8s %-8s %-8s %-8s %-6s %14.6g%14.6g%8i\n'

            for atom in atoms:
                psf.write(afmt % (atom[0], atom[1], atom[2], atom[3], atom[4],
                                  atom[5], atom[6], atom[7], 0.0) )

            psf.write('\n%10i !NBOND\n' % len(bonds) )

            for i, bp in enumerate(bonds):
                psf.write('%10i%10i' % (bp[0], bp[1]) )

                if not (i + 1) % 4:
                    psf.write('\n')

            if (i + 1) % 4:
                psf.write('\n')

            psf.write('\n%10i !NTHETA\n' % len(angles) )
    
            for i, at in enumerate(angles):
                psf.write('%10i%10i%10i' % (at[0], at[1], at[2]) )

                if not (i + 1) % 3:
                    psf.write('\n')

            if (i + 1) % 3:
                psf.write('\n')

            psf.write('\n%10i !NPHI\n' % len(dihedrals) )
    
            for i, dq in enumerate(sorted(dihedrals) ):
                psf.write('%10i%10i%10i%10i' % (dq[0], dq[1], dq[2], dq[3]) )

                if not (i + 1) % 2:
                    psf.write('\n')

            if (i + 1) % 2:
                psf.write('\n')

            psf.write('\n%10i !NIMPHI\n' % len(impropers) )
    
            for i, dq in enumerate(impropers):
                psf.write('%10i%10i%10i%10i' % (dq[0], dq[1], dq[2], dq[3]) )

                if not (i + 1) % 2:
                    psf.write('\n')

            if (i + 1) % 2:
                psf.write('\n')

            psf.write('\n%10i !NDON\n%10i !NACC\n' % (0, 0) )
            psf.write('\n%10i !NNB\n\n' % 0)

            # iblo (exclusion pointer): one entry for each atom
            for i in range(0, self.tot_natoms):
               psf.write('%10i' % 0)

               if not (i + 1) % 8:
                    psf.write('\n\n')

            if (i + 1) % 8:
                psf.write('\n')

            psf.write('\n%10i%10i !NGRP NST2\n' % (len(groups), 0) )

            for i, group in enumerate(groups):
                psf.write('%10i%10i%10i' % (group[0], group[1], 0) )

                if not (i + 1) % 3:
                    psf.write('\n')

            if (i + 1) % 3:
                psf.write('\n')

            # molnt?
            psf.write('\n%10i !MOLNT\n' % 1)

            for i in range(0, self.tot_natoms):
               psf.write('%10i' % 1)

               if not (i + 1) % 8:
                    psf.write('\n')

            if (i + 1) % 8:
                psf.write('\n')

            psf.write('\n%10i%10i !NUMLP NUMLPH\n' % (0, 0) )

        atomno = 0

        with open(rtfname, 'w') as prm:
            prm.write('* created by FESetup\n'
                      '* minimal RTF\n'
                      '*\n'
                      '36 1\n\n')

            for atom, param in atom_params.iteritems():
                atomno += 1
                prm.write('MASS %5i %-6s %9.5f\n' % (atomno, atom, param[0]) )

            prm.write('\nEND\n')

        with open(prmname, 'w') as prm:
            prm.write('* created by FESetup\n')

            prm.write('\nBONDS\n')
            visited = set()

            for n, p in bond_params.iteritems():
                visited.add(n)

                if n[0] == n[1] or (n[1], n[0]) not in visited:
                    prm.write('%-6s %-6s %7.2f %10.4f\n' % (n[0], n[1],
                                                            p[0], p[1]) )

            prm.write('\nTHETAS\n')
            visited = set()

            for n, p in angle_params.iteritems():
                visited.add(n)

                if (n[2], n[1], n[0]) not in visited:
                    prm.write('%-6s %-6s %-6s %7.2f %10.4f\n' %
                              (n[0], n[1], n[2], p[0], p[1]) )

            prm.write('\nPHI\n')
            visited = set()

            for n, terms in dihedral_params.iteritems():
                visited.add(n)

                if (n[3], n[2], n[1], n[0]) not in visited:
                    for term in terms:
                        prm.write('%-6s %-6s %-6s %-6s %10.4f %4i %10.4f\n' %
                                  (n[0], n[1], n[2], n[3],
                                   term[0], term[1], term[2] * const.RAD2DEG) )

            prm.write('\nIMPHI\n')
            visited = set()

            for n, term in improper_params.iteritems():
                visited.add(n)

                if (n[3], n[2], n[1], n[0]) not in visited:
                    prm.write('%-6s %-6s %-6s %-6s %10.4f %4i %10.4f\n' %
                              (n[0], n[1], n[2], n[3],
                               term[0], term[1], term[2] * const.RAD2DEG) )

            prm.write('''
NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -
     CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4
''')
            for atom, param in atom_params.iteritems():
                epsilon = -param[1].epsilon().value()
                sigma = param[1].sigma().value() * const.RSTAR_CONV
                prm.write('%-6s %3.1f %9.6f %9.4f %3.1f %9.6f %9.4f\n' %
                          (atom, 0.0, epsilon, sigma,
                                 0.0, epsilon / 2.0, sigma) )

            prm.write('\nEND\n')


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
    top.writePrmPsf('test.rtf', 'test.prm', 'test.psf')
 
