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
Common methods for the force field ligand classes.  These methods are expected
to be directly imported into the Ligand class.
"""

__revision__ = "$Id$"



import os, sys, re

import pybel
import openbabel as ob

from FESetup import const, errors, report, logger



### not-imported functions ###
def ob_read_one(conv, filen, mol, ifmt, ofmt):
    """
    Use Openbabel to read one structure from a molecular structure file.

    :param conv: Openbanel conversion object
    :type conv: OBConversion
    :param filen: molecular structure file name
    :type filen: string
    :param mol: Openbabel molecule
    :type mol: OBMol
    :param fmt: molecular structure file format known to Openbabel
    :type fmt: string
    :raises: SetupError
    """

    if not conv.SetInAndOutFormats(ifmt, ofmt):
        raise errors.SetupError('conversion from %s to %s not supported '
                                'by Openbabel' % (ifmt, ofmt) )

    # set error level to avoid warnings about non-standard input
    errlev = ob.obErrorLog.GetOutputLevel()
    ob.obErrorLog.SetOutputLevel(0)

    try:
        success = conv.ReadFile(mol, filen)
    except IOError as why:
        raise errors.SetupError(why)

    ob.obErrorLog.SetOutputLevel(errlev)

    if not success:
        txt = 'cannot read %s (%s format), ' % (filen, ifmt)

        if not os.path.isfile(filen):
            raise errors.SetupError(txt + 'file does not exist')
        else:
            raise errors.SetupError(txt + 'check if file/format is valid')

    return mol, conv


### imported methods ###
@report
def prepare(self, to_format = 'mol2', addH = False, calc_charge = False,
            correct_for_pH = False, pH = 7.4):
    """
    Read a molecule structure with Openbabel and get basic information from
    it.  If requested modifiy the data and convert the file.

    :param to_format: format to convert to if not empty
    :type to_format: string
    :param addH: add hydrogens? (experimental)
    :type addH: bool
    :param calc_charge: force total formal charge calculation if True,
                        leave it to Openbabel if False
    :type calc_charge: bool
    :param correct_for_pH: correct for pH, i.e. determine protonation state?
                           (very experimental!)
    :type correct_for_pH: bool
    :param pH: pH to be considered when correct_for_pH is True
               (very experimental!)
    :type pH: float
    :raises: SetupError
    """

    conv = ob.OBConversion()
    mol = ob.OBMol()
    ob_read_one(conv, self.mol_file, mol, self.mol_fmt, to_format)

    # FIXME: does this test for the right thing?
    if mol.GetDimension() != 3:
        raise errors.SetupError('input cooridnates (%s) must have 3 dimensions'
                                % self.mol_file)

    orig_atoms = []

    etab = ob.OBElementTable()
    acnt = 0

    for atom in ob.OBMolAtomIter(mol):
        acnt += 1

        res = atom.GetResidue()     # when is this NULL?
        element = etab.GetSymbol(atom.GetAtomicNum() )
        orig_atoms.append( (res.GetAtomID(atom).strip(), element) )

        if to_format:
            res.SetAtomID(atom, '%s%d' % (element, acnt) )

    if addH:
        logger.write('Adding hydrogens')

        mol.DeleteHydrogens()

        if correct_for_pH:              # reimplementation from phmodel.cpp
            mol.SetAutomaticFormalCharge(correct_for_pH)

            transform = ob.OBChemTsfm()

            logger.write('Applying transforms')

            for reactant, product, pKa in const.TRANSFORM_SMARTS:
                success = transform.Init(reactant, product)

                if not success:
                    raise ValueError('BUG in SMARTS transform Init %s >> %s',
                                     (reactant, product) )

                if (transform.IsAcid() and pH > pKa) or \
                       (transform.IsBase() and pH < pKa):
                    success = transform.Apply(mol)

        mol.AddHydrogens(False, False, 0.0) # valence filling


    if mol.GetTotalSpinMultiplicity() != 1:
        raise errors.SetupError('only multiplicity=1 supported (%s)'
                                % self.mol_file)

    if to_format:
        logger.write('Converting/Writing %s (%s format) to %s format' %
                     (self.mol_file, self.mol_fmt, to_format) )

        # NOTE: this relies on a modified Openbabel MOL2 writer
        conv.AddOption('r', ob.OBConversion.OUTOPTIONS)  # do not append resnum

        self.mol_file = const.CONV_MOL2_FILE % to_format
        self.mol_fmt = to_format

        try:
            conv.WriteFile(mol, self.mol_file)
        except IOError as why:
            raise errors.SetupError(why)
    else:
        logger.write('Leaving %s unmodified.' % self.mol_file)


    # some formats allow formal charge definition:
    # PDB (col 79-80), SDF(M CHG, atom block), MOL2 (UNITY_ATOM_ATTR?)
    if calc_charge:
        logger.write('Computing total formal charge\n')

        # set formal charge for some functional groups by hand because
        # OpenBabel doesn't support it in C++ for formats like mol2
        formal_charge = 0
        sm = ob.OBSmartsPattern()

        for smarts, fchg in const.CHARGE_SMARTS:
            sm.Init(smarts)

            if sm.Match(mol):
                m = list(sm.GetUMapList() )
                formal_charge += fchg * len(m)

        mol.SetTotalCharge(formal_charge)
        self.charge = formal_charge
    else:
        logger.write('Total formal charge taken from coordinate file\n')
        self.charge = mol.GetTotalCharge()  # trust Openbabel...

    self.atomtype = 'sybyl'     # Openbabel will try to convert to Sybyl format

    conv.SetOutFormat('smi')
    conv.AddOption('n')
    conv.AddOption('c')
    smiles = conv.WriteString(mol).rstrip()

    conv.SetOutFormat('inchikey')
    errlev = ob.obErrorLog.GetOutputLevel()
    ob.obErrorLog.SetOutputLevel(0)

    inchi_key = conv.WriteString(mol).rstrip()

    ob.obErrorLog.SetOutputLevel(errlev)

    logger.write('''Ligand key data:
Formula: %s
SMILES: %s
InChIKey: %s
Net charge: %i
Molecular weight: %f\n''' % (mol.GetFormula(), smiles, inchi_key, self.charge,
                             mol.GetMolWt() ) )


@report
def preminimize(self, add_hyd = False, ffield = 'mmff94', nsteps = 10):
    """
    Preminimize mol2 file using the small molecule force fields of OpenBabel.
    Typically usage is to avoid convergence problems with sqm due to a 'bad'
    starting structure.

    :param add_hyd: add hydrogens
    :type add_hyd: bool
    :param ffield: small molecule force field supported by OpenBabel
    :type ffield: string
    :param nsteps: number of minimisation steps
    :type nsteps: integer
    :raises: SetupError
    """

    logger.write('Minimizing %s (%s format)' %
                 (self.mol_file, self.mol_fmt) )

    # NOTE: Openbabel may miscalculate charges from mol2 files
    try:
        mol = pybel.readfile(self.mol_fmt, self.mol_file).next()
    except IOError as why:
        raise errors.SetupError(why)

    if add_hyd:
        logger.write('Adding hydrogens')
        mol.addh()

    mol.localopt(forcefield = ffield, steps = nsteps)

    logger.write('Writing structure to %s (%s format)' %
                 (self.mol_file, self.mol_fmt) )

    try:
        mol.write(self.mol_fmt, self.mol_file, overwrite = True)
    except IOError as why:
        raise errors.SetupError(why)

    self.ref_file = self.mol_file
    self.ref_fmt = self.mol_fmt


@report
def flex(self, name = const.LIGAND_NAME, dobonds = True,
         doangles = True, dodihedrals = True):
    """
    Create a Sire flexibility file describing how the input ligand can be moved.

    :param name: ligand name
    :type name: string
    :param dobonds: include bonds
    :type dobonds: bool
    :param doangles: include angles
    :type doangles: bool
    :param dodihedrals: include dihedrals
    :type dodihedrals: bool
    :raises: SetupError, Sire error
    """

    import Sire.IO
    import Sire.Units
    import Sire.Mol


    re_sire_error = re.compile(const.RE_SIRE_ERROR_STR)

    amber = Sire.IO.Amber()
    molecules = amber.readCrdTop(self.amber_crd, self.amber_top)[0]

    nmol = molecules.molNums()
    nmol.sort()

    # we assume ligand is molecule 0
    solute = molecules.molecule(nmol[0]).molecule()
    solute = solute.edit().rename(name).commit()

    ### guess internals
    connectivity = solute.property('connectivity')
    all_bonds = connectivity.getBonds()
    all_angles = connectivity.getAngles()
    all_dihedrals = connectivity.getDihedrals()

    angle_deltas = {}
    bond_deltas = {}

    logger.write('Computing flexible ligand residues from %s' % self.amber_crd)

    # JM 02/13 With the new Sire code (r1832 onwards) there won't be a ring
    # error so some of the code below could be cleaned up, but I won't right
    # now to maintain compatibility with older versions

    # Redundant torsions are discarded according to the following algorithm
    # 1) Do not sample a torsion at0-at1-at2-at3 if a variable torsion has
    # already been defined around at1-at2 or at2-at1.
    # 2) Do not sample a torsion if it would break a ring
    #
    var_dihedrals = []

    if dodihedrals:
        for dihedral in all_dihedrals:
            at0 = dihedral.atom0()
            at1 = dihedral.atom1()
            at2 = dihedral.atom2()

            move = True

            # see if one of the variable dihedral already rotates around
            # the same torsion
            # --> Note to self This behavior will cause us to skip some
            # torsions (for instance won't sample separately
            # 1-2-3-4 , 1-2-3-5 and 1-2-3-6 with 4-5-6 attached to 3
            #
            # JM 03/13
            # Commented out code below. We will end up with multiple torsions
            # that are sampled around the same central bond Sire can choose to
            # either rotate all atoms around the bond, or sample one torsion
            # separately.
            #for vardih in var_dihedrals:
            #    if ( (at1 == vardih.atom1() and at2 == vardih.atom2() ) or
            #         (at2 == vardih.atom1() and at1 == vardih.atom2() ) ):
            #        # Yes so will not move this torsion
            #        move = False
            #        break
            for vardih in var_dihedrals:
                if dihedral == vardih or dihedral == vardih.mirror():
                    move = False
                    break

            # see if a rotation around this dihedral would break a ring
            if move:
                try:
                    dihbond = Sire.Mol.BondID(at1, at2)
                    solute.move().change(dihbond, 1.0 * Sire.Units.degrees)
                except UserWarning as error:
                    error_type = re_sire_error.search(str(error)).group(0)

                    if error_type == 'SireMol::ring_error':
                        #move = False
                        continue
                    else:
                        raise error

            if move:
                var_dihedrals.append(dihedral)

                # find out how many atoms would move
                gr0, gr1 = connectivity.split(at1, at2)
                ngr0 = gr0.nSelected()
                ngr1 = gr1.nSelected()

                if (ngr0 <= ngr1):
                    smallgroup = gr0
                else:
                    smallgroup = gr1

                smallgroup = smallgroup.subtract(at1)
                smallgroup = smallgroup.subtract(at2)

                nmoved = smallgroup.nSelected()
                if nmoved < 1:
                    nmoved = 1

                if ( connectivity.inRing(dihedral) ):
                    angle_deltas[dihedral] = const.BASE_DIHEDRALRING_FLEX / nmoved
                else:
                    angle_deltas[dihedral] = const.BASE_DIHEDRAL_FLEX / nmoved

    ## angles
    var_angles = []

    # JM 02/13. The strategy to not select angles that move at0 and at2 if they
    # can be both moved by other angles lead to problematic cases where
    # structural readjustments for some topologies become impossible or very
    # hard to achieve. Current approach is to avoid trying to be clever and
    # sample everything
    if doangles:
        for angle in all_angles:
            at0 = angle.atom0()
            at2 = angle.atom2()

            at1 = angle.atom1()

            # test if the angle breaks a ring
            try:
                solute.move().change(angle, 1.0 * Sire.Units.degrees)
            except UserWarning as error:
                error_type = re_sire_error.search(str(error)).group(0)

                if error_type == 'SireMol::ring_error':
                    continue
                else:
                    raise error

            var_angles.append(angle)

            gr0, gr1 = connectivity.split(at0, angle.atom1(), at2)
            ngr0 = gr0.nSelected()
            ngr1 = gr1.nSelected()

            if (ngr0 <= ngr1):
                smallgroup = gr0
            else:
                smallgroup = gr1

            if ( connectivity.inRing(angle) ):
                angle_deltas[angle] = const.BASE_ANGLERING_FLEX / smallgroup.nSelected()
            else:
                angle_deltas[angle] = const.BASE_ANGLE_FLEX / smallgroup.nSelected()

    ## bonds
    var_bonds = []

    if dobonds:
        for bond in all_bonds:
            try:
                solute.move().change(bond, 1.0 * Sire.Units.angstrom)
            except UserWarning as error:
                error_type = re_sire_error.search(str(error)).group(0)

                if error_type == 'SireMol::ring_error':
                    continue
                else:
                    raise error

            var_bonds.append(bond)

            gr0, gr1 = connectivity.split(bond.atom0(), bond.atom1() )
            ngr0 = gr0.nSelected()
            ngr1 = gr1.nSelected()

            if (ngr0 <= ngr1):
                smallgroup = gr0
            else:
                smallgroup = gr1

            if ( connectivity.inRing(bond) ):
                bond_deltas[bond] = const.BASE_BONDRING_FLEX / smallgroup.nSelected()
            else:
                bond_deltas[bond] = const.BASE_BOND_FLEX / smallgroup.nSelected()


    ### guess translation
    translation = const.BASE_TRANSLATION / (solute.nAtoms() / 5.0 + 1.0)


    ### guess rotation
    sphere_radius = solute.evaluate().boundingSphere().radius()
    rotation = const.BASE_ROTATION / sphere_radius**2


    ### write flexibility template
    solute_name = const.LIGAND_NAME

    lines = ['''version 1
molecule %s
rigidbody rotate %-5.3f translate %-5.3f
maximumbondvariables %-3d
maximumanglevariables %-3d
maximumdihedralvariables %-3d\n''' % (solute_name, rotation.to(Sire.Units.degrees),
translation.to(Sire.Units.angstrom), const.BASE_MAXBONDVAR, const.BASE_MAXANGLEVAR, const.BASE_MAXDIHEDRALVAR) ]

    for bond in var_bonds:
        at0name = solute.select(bond.atom0()).name().value()
        at1name = solute.select(bond.atom1()).name().value()
        delta = bond_deltas[bond].to(Sire.Units.angstrom)

        lines.append('bond %-4s %-4s flex %-5.3f\n' \
                     % (at0name, at1name, delta))

    for angle in var_angles:
        at0name = solute.select(angle.atom0()).name().value()
        at1name = solute.select(angle.atom1()).name().value()
        at2name = solute.select(angle.atom2()).name().value()
        delta = angle_deltas[angle].to(Sire.Units.degrees)

        lines.append('angle %-4s %-4s %-4s flex %-5.3f\n' \
                 % (at0name, at1name, at2name, delta))

    for dihedral in var_dihedrals:
        at0name = solute.select(dihedral.atom0()).name().value()
        at1name = solute.select(dihedral.atom1()).name().value()
        at2name = solute.select(dihedral.atom2()).name().value()
        at3name = solute.select(dihedral.atom3()).name().value()
        delta = angle_deltas[dihedral].to(Sire.Units.degrees)

        lines.append('dihedral %-4s %-4s %-4s %-4s flex %-5.3f\n' \
                 % (at0name, at1name, at2name, at3name, delta))

    with open(const.LIGAND_FLEX_FILE, 'w') as output:
            output.write(''.join(lines))


@report
def conf_search(self, do_ga = False, numconf = 30, geomsteps = 5,
                numchildren = 5, mutability = 5, convergence = 25,
                ffield = 'mmff94', steep_steps = 100, steep_econv = 1.0E-4,
                conj_steps = 250, conj_econv = 1.0E-6):
    """Conformer search via OpenBabel.

    The default method is a weighted rotor search with pre (steepest descent)
    and post (conjugate gradient) minimization steps.  The new genetic
    algorithm introduced in Openbabel 2.3 is accessible here too but the
    Python code crashes at present with segfaults when writing the output
    file.

    The output file is written to const.FREE_MOL2_FILE in MOL2/Sybyl format.

    :raises: SetupError
    """

    try:
        mol = pybel.readfile(self.mol_fmt, self.mol_file).next()
    except IOError as why:
        raise errors.SetupError(why)

    obmol = mol.OBMol

    self.mol_fmt = 'mol2'
    outmol = const.FREE_MOL2_FILE % self.mol_fmt

    if not do_ga:
        obff = ob.OBForceField.FindForceField(ffield)

        logger.write('Performing conformer search with pre and post '
                     'minimisation on %s' % self.mol_file)

        if not obff:
            raise errors.SetupError('Error: cannot find %s force field' %
                                    ffield)

        if not obff.Setup(obmol):
            raise errors.SetupError('Error: cannot setup mol with %s force '
                                    'field' % ffield)

        logger.write('Parameters are:\nffield = %s, numconf = %i, '
                     'geomsteps = %i, steep_steps = %i,\nsteep_econv = %g, '
                     'conj_steps = %i, conj_econv = %g\n' %
                     (ffield, numconf, geomsteps, steep_steps, steep_econv,
                      conj_steps, conj_econv) )

        obff.SteepestDescent(steep_steps, steep_econv)
        obff.WeightedRotorSearch(numconf, geomsteps)
        obff.ConjugateGradients(conj_steps, conj_econv)
        obff.GetCoordinates(obmol)
    else:
        search = ob.OBConformerSearch()
        cando = search.Setup(obmol, numconf, numchildren, mutability,
                         convergence)

        if not cando:
            logger.write('Not performing conformer search: only one '
                              'possible rotamer')
            return
        else:
            logger.write('Performing conformer search on %s' %
                              self.mol_file)

            score = ob.OBEnergyConformerScore()
            score.Preferred = 1     # enum{HighScore, LowScore}
            score.Convergence = 1   # enum{Highest, Lowest, Sum, Average}

            search.SetScore(score)
            search.Search()

            search.GetConformers(obmol)

            logger.write('%i conformers generated.  Writing highest '
                              'scored to %s' % (obmol.NumConformers(), outmol) )

            obmol.SetConformer(0)

    logger.write('Writing output file in %s format to %s' %
                 (self.mol_fmt, outmol) )

    try:
        mol.write(self.mol_fmt, outmol, overwrite = True)
    except IOError as why:
        raise errors.SetupError(why)

    self.mol_file = outmol
    self.mol_atomtype = 'sybyl'


@report
def align(self, inc_hyd = False, symmetry = False, filt = False):
    """
    Align two structures via OpenBabel.

    :param inc_hyd: include hydrogens
    :type inc_hyd: bool
    :param symmetry: consider symmetry of the molecule
    :type symmetry: bool
    :param filt: invoke primitive filter
    :type filt: bool
    :raises: SetupError
    """

    if not self.ref_file:
        raise errors.SetupError('reference file not set yet')

    if self.ref_file ==  self.mol_file:
        logger.write('Identical files: will not perform alignment')
        return


    conv = ob.OBConversion()
    conv.SetInAndOutFormats(self.ref_fmt, self.mol_fmt)

    ref = ob.OBMol()
    tgt = ob.OBMol()

    # ignore warning messages about non-standard PDB
    errlev = ob.obErrorLog.GetOutputLevel()
    ob.obErrorLog.SetOutputLevel(0)

    conv.ReadFile(ref, self.ref_file)
    conv.ReadFile(tgt, self.mol_file)

    ob.obErrorLog.SetOutputLevel(errlev)


    if filt:
        # delete unwanted atoms in reference
        delat = []

        for atom in ob.OBMolAtomIter(ref):
            resname = atom.GetResidue().GetName()

            # FIXME: replace with proper function
            if resname in const.IGNORE_RESIDUES:
                delat.append(atom)

        ref.BeginModify()

        for atom in delat:
            ref.DeleteAtom(atom, True)

        ref.EndModify()

        # copy wanted atoms from target to new OBMol
        cpy = ob.OBMol()

        for atom in ob.OBMolAtomIter(tgt):
            resname = atom.GetResidue().GetName()

            # FIXME: replace with proper function
            if resname not in const.IGNORE_RESIDUES:
                # NOTE: this copies only some info but incl. coordinates
                cpy.AddAtom(atom)

        molecs = ob.OBAlign(ref, cpy, inc_hyd, symmetry)
    else:
        molecs = ob.OBAlign(ref, tgt, inc_hyd, symmetry)

    logger.write('Aligning %s with %s as reference' %
                      (self.mol_file, self.ref_file) )

    if not molecs.Align():
        logger.write('Alignment failed')
        return

    logger.write('RMSD is %.2f' % molecs.GetRMSD() )

    if filt:
        rotate = molecs.GetRotMatrix()
        molecs.UpdateCoords(ref)
        first = True

        for atom in ob.OBMolAtomIter(tgt):
            tmpvec = ob.vector3(atom.GetVector())
            tmpvec *= rotate

            if first:
                # we obviously can't do this directly: OB's vector3 does not
                # support the '-' but only the '-=' operator.  But the
                # latter leads to a memory corruption bug...
                # shift = ref.GetAtom(1).GetVector()
                # shift -= cpy.GetAtom(1).GetVector()

                # we assume atoms 1 are equivalent in both molecules...
                v1 = ref.GetAtom(1).GetVector()
                v2 = cpy.GetAtom(1).GetVector()

                x = v1.GetX() - v2.GetX()
                y = v1.GetY() - v2.GetY()
                z = v1.GetZ() - v2.GetZ()

                shift = ob.vector3(x, y, z)
                first = False

            tmpvec += shift
            atom.SetVector(tmpvec)

    else:
        if not molecs.UpdateCoords(tgt):
            logger.write('Coordinate update failed')
            return

    try:
        conv.WriteFile(tgt, self.mol_file)
    except IOError as why:
        raise errors.SetupError(why)

    self.mol_atomtype = 'sybyl'
