#  Copyright (C) 2012-2016  Hannes H Loeffler, Julien Michel
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
A class to build a complex.  Derives from Common.

The complex Setup class composes a complex object from a protein and a ligand
object.  The class can create an AMBER topology file for the complex.
"""


__revision__ = "$Id$"



import FESetup
from FESetup import const, errors, logger
import utils

from ligand import Ligand
from protein import Protein
from common import *



class Complex(Common):
    """The complex setup class."""

    from FESetup.prepare.ligutil import flex as lig_flex

    SSBONDS_OFFSET = 1

    def __init__(self, protein, ligand):
        """
        :param protein: the protein for complex composition
        :type protein: Protein or string
        :param ligand: the ligand for complex composition
        :type ligand: Ligand or string
        :raises: SetupError
        """

        self.leap_added = False

        # FIXME: remove when ModelConfig is done
        #        this is still used for the Morph class
        if type(protein) == str and type(ligand) == str:
            super(Complex, self).__init__(protein + const.PROT_LIG_SEP +
                                          ligand)

            self.protein_file = protein
            self.ligand_file = ligand

            # FIXME: quick fix to allow dGprep to redo complex morph
            self.ligand = Ligand(ligand, '')

            return

        assert type(protein) == Protein
        assert type(ligand) == Ligand

        self.complex_name = protein.mol_name + const.PROT_LIG_SEP + \
                            ligand.mol_name

        super(Complex, self).__init__(self.complex_name)

        self.ligand_file = ligand.orig_file
        self.protein_file = protein.orig_file
        self.charge = protein.charge + ligand.charge

        if abs(self.charge) > const.TINY_CHARGE:
            logger.write('Warning: non-zero complex charge (%f)' % self.charge)

        self.protein = protein
        self.ligand = ligand
        self.frcmod = self.ligand.frcmod


    @report
    def prepare_top(self, gaff='gaff', pert=None, add_frcmods=[]):
        """
        Prepare for parmtop creation i.e. add molecules to Leap structure.
        This needs to be run before create_top() to ensure that the molecule
        has been added but also to not trigger parmtop generation to early.
        Pmemd needs to have a second molecule added in alchemical free
        energy setups.
        """

        # ensure ligand is in MOL2/GAFF format
        if os.access(const.LIGAND_AC_FILE, os.F_OK):
            mol_file = const.GAFF_MOL2_FILE
            antechamber = utils.check_amber('antechamber')
            utils.run_amber(antechamber,
                            '-i %s -fi ac -o %s -fo mol2 -j 1 -at %s -pf y' %
                            (const.LIGAND_AC_FILE, mol_file, gaff) )
            self.ligand_fmt = 'mol2'
        else:
            # antechamber has trouble with dummy atoms
            mol_file = self.ligand_file

            if self.ligand_fmt != 'mol2' and self.ligand_fmt != 'pdb':
                raise errors.SetupError('unsupported leap input format: %s ('
                                        'only mol2 and pdb)' % self.ligand_fmt)


        frcmods = [self.frcmod]

        if add_frcmods:
            frcmods.extend(add_frcmods)

        if not self.leap_added:
            self.leap.add_force_field(gaff)
            self.leap.add_mol(mol_file, self.ligand_fmt, frcmods, pert=pert)


    @report
    def create_top(self, boxtype='', boxlength=10.0, align=False,
                   neutralize=False, addcmd='', addcmd2='',
                   remove_first=False, conc = 0.0, dens = 1.0):
        """Generate an AMBER topology file via leap.

        :param boxtype: rectangular, octahedron or set (set dimensions explicitly)
        :param boxlength: side length of the box
        :param align: align solute along the principal axes
        :param neutralize: neutralise the system
        :param addcmd: inject additional leap commands
        :param remove_first: remove first unit/residue
        :param conc: ion concentration
        :type conc: float
        :param dens: expected target density
        :type dens: float
        :type boxtype: string
        :type boxlength: float
        :type align: bool
        :type neutralize: int
        :type addcmd: string
        :type remove_first: bool
        """

        if not self.leap_added:
            self.leap.add_mol(self.protein_file, 'pdb')
            self.leap_added = True

        # FIXME: there can be problems with the ordering of commands, e.g.
        #        when tip4pew is used the frcmod files are only loaded after
        #        reading PDB and MOL2
        leapin = self._amber_top_common(boxtype, boxlength,
                                        neutralize, align=align,
                                        remove_first=remove_first,
                                        conc=conc, dens=dens)

        utils.run_leap(self.amber_top, self.amber_crd, 'tleap', leapin)


    @report
    def prot_flex(self, cut_sidechain = 15.0, cut_backbone = 15.0):
        """
        Create a flexibility file for the protein describing how the input
        molecule can be moved by Sire.

        :param cut_sidechain: side chain cutoff
        :type cut_sidechain: float
        :param cut_backbone: backbone cutoff
        :type cut_backbone: float
        :raises: SetupError
        """


        if cut_sidechain < 0.0 or cut_backbone < 0.0:
            raise errors.SetupError('Cutoffs must be positive')


        import Sire.IO

        amber = Sire.IO.Amber()
        molecules, space = amber.readCrdTop(self.sander_crd, self.amber_top)

        moleculeNumbers = molecules.molNums()
        moleculeNumbers.sort()
        moleculeList = []

        for moleculeNumber in moleculeNumbers:
            molecule = molecules.molecule(moleculeNumber).molecule()
            moleculeList.append(molecule)

        ligand = moleculeList[0]
        not_ligand = moleculeList[1:]

        sc_bb_residues = []

        logger.write('Computing flexible protein residues from %s' %
                     self.sander_crd)

        for cut in cut_sidechain, cut_backbone:
            cut_residues = []
            cut2 = cut**2

            for molecule in not_ligand:

                #for residue in molecule.residues():
                nmolresidues = molecule.nResidues()
                for z in range(0,nmolresidues):
                    residue = molecule.residues()[z]
                    # FIXME: a better way to skip unwanted residues would be to
                    # examine amber.zmatrices directly
                    # .value() returns a QtString!
                    if (str(residue.name().value()) not in
                        const.AMBER_PROTEIN_RESIDUES):
                        continue

                    shortest_dist2 = float('inf')

                    #for resat in residue.atoms():
                    nresatoms = residue.nAtoms()
                    for x in range(0,nresatoms):
                        resat = residue.atoms()[x]

                        if (resat.property('mass').value() <
                            const.MAX_HYDROGEN_MASS):
                            continue

                        rescoords = resat.property('coordinates')

                        #for ligat in ligand.atoms():
                        nligatoms = ligand.nAtoms()
                        for y in range(0,nligatoms):
                            ligat = ligand.atoms()[y]

                            if (ligat.property('mass').value() <
                                const.MAX_HYDROGEN_MASS):
                                continue

                            ligcoords = ligat.property('coordinates')
                            dist2 = space.calcDist2(rescoords, ligcoords)

                            if dist2 < shortest_dist2:
                                shortest_dist2 = dist2

                    if shortest_dist2 < cut2:
                        cut_residues.append(residue)

            sc_bb_residues.append(cut_residues)

        lines = ['''# Flexible residues were only selected from the following list of residue names
# %s
# Cut-off for selection of flexible side chains %s Angstroms
# Number of residues with flexible side chains: %s
# Cut-off for selection of flexible backbone: %s Angstroms
# Number of residues with flexible backbone: %s
''' % (', '.join(const.AMBER_PROTEIN_RESIDUES), cut_sidechain,
       len(sc_bb_residues[0]), cut_backbone, len(sc_bb_residues[1])) ]

        htext = ['flexible sidechain', 'flexible backbone']

        for i in 0, 1:
            lines.append("%s\n" % htext[i])
            nums = []

            for residue in sc_bb_residues[i]:
                nums.append(residue.number().value() )

            nums.sort()

            line = ''
            for num in nums:
                if len(line) > 75:
                    lines.append('%s\n' % line)
                    line = ''
                line += ' %4d' % num

            lines.append('%s\n' % line)

        with open(const.PROTEIN_FLEX_FILE, 'w') as output:
            output.write(''.join(lines))



if __name__ == '__main__':
    pass
