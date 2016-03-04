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

r"""
Create perturbed topologies.
"""


__revision__ = "$Id$"


import os, sys
from collections import defaultdict

import sander
from FESetup import const, errors, logger
from FESetup.prepare.amber import gromacs



MORPH_TOP = 'morph.top'
MORPH_GRO = 'morph.gro'
VAC_MDP_FILE = '_vac.mdp'
SOL_MDP_FILE = '_sol.mdp'


class PertTopology(object):

    def __init__(self, FE_sub_type, sc_type, ff, con_morph, atoms_initial,
                 atoms_final, lig_initial, lig_final, atom_map,
                 reverse_atom_map, zz_atoms, gaff):

        self.sc_type = sc_type
        self.ff = ff
        self.gaff = gaff
        self.con_morph = con_morph
        self.atoms_initial = atoms_initial
        self.atoms_final = atoms_final
        self.lig_initial = lig_initial
        self.lig_final = lig_final
        self.atom_map = atom_map
        self.reverse_atom_map = reverse_atom_map
        self.zz_atoms = zz_atoms

        self.topol = None
        

    def setup(self, curr_dir, lig_morph, cmd1, cmd2):

        topol = sander.PertTopology('dummy', self.sc_type, self.ff,
                                    self.con_morph, self.atoms_initial,
                                    self.atoms_final, self.lig_initial,
                                    self.lig_final, self.atom_map,
                                    self.reverse_atom_map, self.zz_atoms,
                                    self.gaff)

        topol.setup(curr_dir, lig_morph, cmd1, cmd2)

        lig0 = topol.lig0._parm_overwrite
        lig1 = topol.lig1._parm_overwrite
        top = topol.lig0.TOP_EXT
        rst = topol.lig0.RST_EXT

        # top0 and top1 written mainly for debugging purposes
        top0 = gromacs.GromacsTop()
        top0.readParm(lig0 + top, lig0 + rst)
        top0.writeTop(lig0 + const.GROMACS_ITP_EXT, lig0 + '.atp')
        top0.writeGro(lig0 + const.GROMACS_GRO_EXT)

        top1 = gromacs.GromacsTop()
        top1.readParm(lig1 + top, lig1 + rst)
        top1.writeTop(lig1 + const.GROMACS_ITP_EXT, lig1 + '.atp')
        top1.writeGro(lig1 + const.GROMACS_GRO_EXT)

        gromacs.mixer(top0, top1)

        self.topol = topol

        # FIXME: ugly kludges!
        if not os.access(MORPH_GRO, os.F_OK):
            os.symlink(lig0 + const.GROMACS_GRO_EXT, MORPH_GRO)

        with open(MORPH_TOP, 'w') as mtop:
            mtop.write('''
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

[ atomtypes ]
#include "%s"

#include "%s"

[ system ]
FESetup

[ molecules ]
; Compound        nmols
%s 1
''' % (const.GROMACS_PERT_ATP, const.GROMACS_PERT_ITP, const.LIGAND_NAME) )

        with open(VAC_MDP_FILE, 'w') as mdp:
            mdp.write(VAC_MDP)


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2):
        """
        Create only topology file, not coordinates.
        """

        self.topol.create_coords(curr_dir, dir_name, lig_morph, pdb_file,
                                 system, cmd1, cmd2)

        # FIXME: ugly kludge, assuming the file is one level up
        if not os.access(const.GROMACS_PERT_ATP, os.F_OK):
            os.symlink('../%s' % const.GROMACS_PERT_ATP,
                       '%s' % const.GROMACS_PERT_ATP)

        # FIXME: ugly kludge, assuming the file is one level up
        if not os.access(const.GROMACS_PERT_ITP, os.F_OK):
            os.symlink('../%s' % const.GROMACS_PERT_ITP,
                       '%s' % const.GROMACS_PERT_ITP)

        top = gromacs.GromacsTop()
        top.readParm(self.topol.parmtop, self.topol.inpcrd)

        # FIXME: need to reparse ATP file
        with open(const.GROMACS_PERT_ATP, 'r') as atp:
            atomtypes = []

            for line in atp:
                tmp = line.split()
                atomtypes.append( (tmp[0], float(tmp[2]), float(tmp[5]),
                                   float(tmp[6]) ) )
 
        top.addAtomTypes(atomtypes)
        top.writeTop(MORPH_TOP, '', const.LIGAND_NAME, False)
        top.writeGro(MORPH_GRO)

        with open(SOL_MDP_FILE, 'w') as mdp:
            mdp.write(SOL_MDP)



VAC_MDP = '''; TI template for vacuum
integrator               = sd
ld-seed                  = -1
bd-fric                  = 0
dt                       = 0.001
nsteps                   = 2000000
nstcomm                  = 100
comm-mode                = Angular

nstxout                  = 10000
nstvout                  = 10000
nstfout                  = 0
nstlog                   = 10000
nstenergy                = 10000
nstxtcout                = 0

cutoff-scheme            = group
nstlist                  = 0
ns_type                  = simple
pbc                      = no
rlist                    = 0.0

coulombtype              = cut-off
coulomb-modifier         = none
rcoulomb                 = 0.0

vdwtype                  = cut-off
vdw-modifier             = none
rvdw                     = 0.0
DispCorr                 = no

tcoupl                   = v-rescale
nsttcouple               = 10
tc_grps                  = System
tau_t                    = 2.0
ref_t                    = 298.0

pcoupl                   = No

gen-vel                  = yes
gen-temp                 = 298.0
gen-seed                 = -1

constraints              = none
constraint_algorithm     = Lincs
continuation             = no
lincs_order              = 4
lincs_warnangle          = 30

free-energy              = yes
delta-lambda             = 0
init-lambda-state        = %L%
fep-lambdas              = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
nstdhdl                  = 100
calc-lambda-neighbors    = -1
sc-alpha                 = 0.5
sc-coul                  = yes
sc-power                 = 1.0
sc-r-power               = 6
sc-sigma                 = 0.3
dhdl-derivatives         = yes
dhdl-print-energy        = no
separate-dhdl-file       = yes
dh_hist_size             = 0
dh_hist_spacing          = 0.1
'''

SOL_MDP = '''; TI template for solution
integrator               = sd
ld-seed                  = -1
bd-fric                  = 0
dt                       = 0.001
nsteps                   = 500000
nstcomm                  = 1000
comm-mode                = Linear

nstxout                  = 10000
nstvout                  = 10000
nstfout                  = 0
nstlog                   = 10000
nstenergy                = 10000
nstxtcout                = 0

cutoff-scheme            = group
nstcalclr                = 1
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 0.8

coulombtype              = PME
coulomb-modifier         = none
rcoulomb                 = 0.8
fourierspacing           = 0.10
pme_order                = 4
ewald_rtol               = 1.0E-5
optimize_fft             = yes

vdwtype                  = cut-off
vdw-modifier             = none
rvdw                     = 0.8
DispCorr                 = AllEnerPres

tcoupl                   = v-rescale
tc_grps                  = System
tau_t                    = 2.0
ref_t                    = 298.0

pcoupl                   = Berendsen
pcoupltype               = isotropic
tau_p                    = 5.0
compressibility          = 4.5e-5
ref_p                    = 1.0
refcoord-scaling         = com

gen-vel                  = no

constraints              = h-bonds
constraint_algorithm     = Lincs
continuation             = yes
lincs_order              = 4
lincs_warnangle          = 30

free-energy              = yes
delta-lambda             = 0
init-lambda-state        = %L%
fep-lambdas              = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
mass-lambdas             = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
nstdhdl                  = 100
calc-lambda-neighbors    = -1
sc-alpha                 = 0.5
sc-coul                  = yes
sc-power                 = 1.0
sc-r-power               = 6
sc-sigma                 = 0.3
dhdl-derivatives         = yes
dhdl-print-energy        = no
separate-dhdl-file       = yes
dh_hist_size             = 0
dh_hist_spacing          = 0.1
'''
