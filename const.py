#  Copyright (C) 2012-2016  Hannes H Loeffler
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

"""Module with global constants for FESetup.

.. moduleauthor: Hannes Loeffler <Hannes.Loeffler@stfc.ac.uk>

"""

__revision__ = "$Id$"



import os

import Sire.Units



AMBER_BIN_PATH = ''

AMBER_PROTEIN_RESIDUES = frozenset(
    ('ALA', 'ARG', 'ASN', 'ASP', 'ASH', 'CYS', 'CYM', 'CYX', 'GLU', 'GLH',
     'GLN', 'GLY', 'HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'ILE',
     'LEU', 'LYS', 'LYN', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR',
     'VAL')
)

AMBER_NUCLEIC_RESIDUES = frozenset(
    ('A', 'A3', 'A5', 'AN', 'C', 'C3', 'C5', 'CN', 'DA', 'DA3', 'DA5', 'DAN',
     'DC', 'DC3', 'DC5', 'DCN', 'DG', 'DG3', 'DG5', 'DGN', 'DT', 'DT3', 'DT5',
     'DTN', 'G', 'G3', 'G5', 'GN', 'OHE', 'U', 'U3', 'U5', 'UN')
)

AROMATICS = frozenset(('TRP', 'TYR', 'PHE', 'HIS', 'HID', 'HIE', 'HIP'))

TORSIONS_TO_ZERO = frozenset(
    ('CZ-CE1-CD1-CG', 'CE2-CZ-CE1-CD1', 'CD2-CE2-CZ-CE1', 'NE2-CE1-ND1-CG',
     'CD2-NE2-CE1-ND1','CE2-NE1-CD1-CG', 'CZ3-CH2-CZ2-CE2', 'CE3-CZ3-CH2-CZ2',
     'CD2-CE3-CZ3-CH2')
)

TORSIONS_TO_ONEHUNDREDEIGHTY = frozenset(
    ('CE1-CD1-CG-CB', 'HD1-CD1-CG-CE1', 'HE1-CE1-CD1-CZ' 'HZ-CZ-CE1-CE2',
     'HE2-CE2-CZ-CD2', 'HD2-CD2-CE2-CG', 'OH-CZ-CE1-CE2', 'CE1-ND1-CG-CB',
     'HD1-ND1-CG-CE1', 'HE1-CE1-ND1-NE2', 'HD2-CD2-NE2-CG', 'HE2-NE2-CE1-CD2',
     'NE1-CD1-CG-CB', 'CZ2-CE2-NE1-CD1', 'CH2-CZ2-CE2-NE1', 'HD1-CD1-CG-NE1',
     'HE1-NE1-CD1-CE2', 'HZ2-CZ2-CE2-CH2', 'HH2-CH2-CZ2-CZ3', 'HZ3-CZ3-CH2-CE3',
     'HE3-CE3-CZ3-CD2')
)

LIGAND_NAME = 'LIG'
LIGAND0_NAME = 'L0'
LIGAND1_NAME = 'L1'
INT_NAME = 'INT'

LIGAND_INFO_FILE = 'ligand.info'
LIGAND_AC_FILE = 'ligand.ac'
CORR_AC_FILE = 'corr.ac'
CORR_CH_FILE = 'corr.ch'
SQM_PDB_FILE = 'sqm.pdb'
CONV_MOL2_FILE = 'ligand_conv' + os.extsep + '%s'
LIGAND_TMP = 'ligand_tmp'
FREE_MOL2_FILE = 'ligand_free' + os.extsep + '%s'
GB_PREFIX = 'gb_charge'
GB_FRCMOD_FILE = GB_PREFIX + os.extsep + 'frcmod'
LIGAND_FRCMOD_FILE = 'ligand.frcmod'
OFF_FILE = 'ligand.off'
GAFF_MOL2_FILE = 'gaff.mol2'

DUMMY_TYPE = 'du'

IGNORE_RESIDUES = frozenset((LIGAND_NAME, 'TIP', 'WAT', 'Cl-', 'Na+'))

KNOWN_MOL2_ATOMTYPES = ('gaff', 'gaff2', 'sybyl', 'amber')

LIGAND_WORKDIR = '_ligands'
PROTEIN_WORKDIR = '_proteins'
COMPLEX_WORKDIR = '_complexes'
MORPH_WORKDIR = '_perturbations'

PROTEIN_FLEX_FILE = 'protein.flex'
LIGAND_FLEX_FILE = 'ligand.flex'

BASE_DIHEDRALRING_FLEX = 2.5 * Sire.Units.degrees
BASE_DIHEDRAL_FLEX = 30.0 * Sire.Units.degrees
BASE_ANGLERING_FLEX = 0.10 * Sire.Units.degrees
BASE_ANGLE_FLEX = 0.25 * Sire.Units.degrees
BASE_BONDRING_FLEX = 0.020 * Sire.Units.angstrom
BASE_BOND_FLEX = 0.025 * Sire.Units.angstrom
BASE_TRANSLATION = 0.5 * Sire.Units.angstrom
BASE_ROTATION = 30.0 * Sire.Units.degrees
BASE_MAXVAR = 10
BASE_MAXBONDVAR = 5
BASE_MAXANGLEVAR = 5
BASE_MAXDIHEDRALVAR = 5

MAX_CHARGE = 1.0
MAX_CHARGE_DIFF = 0.01
TINY_CHARGE = 0.01

RE_SIRE_ERROR_STR = r'(Sire\w+?::\w+)'

SSBOND_FILE = 'ssbonds'
PROPKA_OUT_FILE = 'propka.dat'

LEAP_IN = 'leap.in'
LEAP_SOLVATED = 'solvated'
LEAP_VACUUM = 'vacuum'
LEAP_IONIZED = 'ionized'

MAX_HYDROGEN_MASS = 2.0

MCS_MAP_FILE = 'mcs_map.pkl'
MCS_MOL_FILE = 'mcs.mol2'
MORPH_NAME = 'MORPH'

PDB_EXT = os.extsep + 'pdb'
PRMTOP_EXT = os.extsep + 'parm7'
INPCRD_EXT = os.extsep + 'rst7'
MOL2_EXT = os.extsep + 'mol2'

MODEL_EXT = os.extsep + 'model'
MODEL_SOLV_PREFIX = 'solv_'

FLAT_RINGS_FILE = 'flatrings' + PDB_EXT
PROTONATED_PDB_FILE = 'protonated' + PDB_EXT

NOT_FIRST_TOP = 'not_first' + PRMTOP_EXT
NOT_FIRST_CRD = 'not_first' + INPCRD_EXT
NOT_FIRST_PDB = 'not_first' + PDB_EXT

DLFIELD_UDFF_NAME = 'dl_field.udff'
DLFIELD_PDB_NAME = 'dl_field' + PDB_EXT
RSTAR_CONV = 0.5612310241546865         # math.pow(2.0, 1.0 / 6.0) / 2.0

SIRE_ABS_PERT_FILE = 'MORPH.onestep.pert'
SIRE_ABS_PERT_EL_FILE = 'MORPH.elec.pert'
SIRE_ABS_PERT_VDW_FILE = 'MORPH.vdw.pert'

GROMACS_GRO_EXT = os.extsep + 'gro'
GROMACS_ITP_EXT = os.extsep + 'itp'
GROMACS_TOP_EXT = os.extsep + 'top'
GROMACS_POSRES_PREFIX = 'posres_'
GROMACS_PERT_ITP = 'pert' + GROMACS_ITP_EXT
GROMACS_PERT_ATP = 'pert' + os.extsep + 'atp'

RAD2DEG = 57.29577951308232             # 180.0 / math.pi
AMBER_VELCONV = 20.455                  # AMBER time in 1/20.455 ps
A2NM = 0.1
CAL2J = 4.184
CALA2JNM2 = CAL2J                       # kcal/mol/A2 to kJ/mol/nm2
ATM2BAR = 1.01325
AMU2GRAMS = 1.6605402                   # * 10^-24 grams

# ':' is a noop but only as a command so should be safe in filenames
PROT_LIG_SEP = ':'
# in the middle of a filename a tilde should be safe
MORPH_SEP = '~'

_phospho = '[#15D4](~[OD1])(~[OD1])(~[OD2])~[OD2]'
_sulfonyl = '[#16](~[OD1])(~[OD1])'

# FIXME: this compensates for missing charge determination by Openbabel, so
#        will fail to detect charges in cases not covered here'; relies on
#        explicit hydrogens!
CHARGE_SMARTS = (
    ('[#16;$([#16D4](~[OD1])(~[OD1])~[OD1])]', -1),  # sulfates, Sac
    ('[#16;$([#16D3](~[OD1])~[OD1])]', -1),  # sulfonates, Sac
    ('[#15;$([#15D4](~[OD1])(~[OD1])~[OD1])]', -2),  # phosphates, Pac
    ('[#15;$([#15D3](~[OD1])~[OD1])]', -2),  # phosphonates, Pac
    ('[#15;$(%s)]' % _phospho, -1),     # polyphosphates or phosphate diesters
    ('[CD3](~[NH2])(~[NH2])~[NH]', +1), # guanidinium
    ('O=C~[CD3,ND2]~C=O', -1),          # beta-dicarbonyls, imides
    ('%s~[ND2]~C(=O)~[NH]' % _sulfonyl, -1), # sulfonylureas
    ('%s~[ND2]~C(=O)~[#6]' % _sulfonyl, -1), # (sulfon)imides
    ('%s~[ND2]~c' % _sulfonyl, -1),     # N-aryl sulfonamides
    ('{0}~[ND2]~{0}'.format(_sulfonyl), -1), # sulfonimides 2
    ('[#6][ND3]([#6])~[#6H]~[ND3]([#6])[#6]', +1), # formamidiniums
    ('C(~[OD1])~[OD1]', -1),            # carbonates, Cac in atomtyp.txt
    ('[#6][#16D1]', -1),                # thiolates
    ('[#6][#8D1]', -1),                 # alkoxides (alcoholates)
    ('[#7D4]', +1),                     # atomtyp.txt checks for [#7X4]
    ('{0}~[ND3](~{0})~[CD3](~{0})~{0}'.format('[#6,#1]'), +1),  # imines
    )

# NOTE: we rely on explicit hydrogens here!
TRANSFORM_SMARTS = (
    ('O=C[OD1-0:1]', 'O=C[O-1:1]', 4.0),  # see phmodel.txt
    ('[N^3;!$(N~[!#6;!#1]):1]', '[N+:1]', 9.8),  # see phmodel.txt
    )
