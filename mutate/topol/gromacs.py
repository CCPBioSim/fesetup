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


import os
import sys
from collections import defaultdict

import sander
from FESetup import const, errors, logger
from FESetup.prepare.amber import gromacs



MORPH_TOP = 'morph.top'
MORPH1_TOP = 'morph1.top'
MORPH2_TOP = 'morph2.top'
MORPH_GRO = 'morph.gro'
VAC_MDP_FILE = '_vac.mdp'
SOL_MDP_FILE = '_sol.mdp'


def _lambda_paths(dummies0, dummies1):
    '''
    Compute lambda paths depending on what state the dummy atoms are in.
    Works for appearing or disappearing atoms only.

    :param initial: the initial state of the morph pair
    :type initial: either Ligand or Complex
    :param final: the final state of the morph pair
    :type final: either Ligand or Complex
    '''

    if not dummies0 and not dummies1:
        fepl = '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0'
        vdwl = '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0'
        masl = '0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
        seps = 'linear transformation only'
    elif dummies0:
        fepl = ('0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.4 '
                '0.6 0.8 1.0')
        vdwl = ('0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.0 1.0 '
                '1.0 1.0 1.0')
        masl = ('0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 '
                '0.0 0.0 0.0')
        seps = 'appearing atoms: vdW before q_on'
    elif dummies1:
        fepl = ('0.0 0.2 0.4 0.6 0.8 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 '
                '1.0 1.0 1.0')
        vdwl = ('0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 '
                '0.8 0.9 1.0')
        masl = ('0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 '
                '0.0 0.0 0.0')
        seps = 'disappearing atoms: q_off before vdW'
    else:
        raise errors.SetupError('BUG: gromacs morph tolopgy - unexpected'
                                'dummies in "dummy" protocol')

    return fepl, vdwl, masl, seps


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

        self.dummies0 = not all([a.atom for a in self.atom_map.keys()])
        self.dummies1 = not all([a.atom for a in self.atom_map.values()])

        self.topol = None
        

    def setup(self, curr_dir, lig_morph, cmd1, cmd2):
        """
        Create TOP and link to GRO.
        """

        if self.dummies0 and self.dummies1:
            self.FE_sub_type = 'dummy3'
        else:
            self.FE_sub_type = 'dummy'

        topol = sander.PertTopology(self.FE_sub_type, self.sc_type, self.ff,
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

        # FIXME: ugly kludges!
        if not os.access(MORPH_GRO, os.F_OK):
            os.symlink(lig0 + const.GROMACS_GRO_EXT, MORPH_GRO)

        if self.FE_sub_type == 'dummy':
            gromacs.mixer(top0, top1)

            with open(MORPH_TOP, 'w') as mtop:
                mtop.write(TOP_TMPL.format(title='one-step TI/FEP',
                                           atp=const.GROMACS_PERT_ATP,
                                           itp=const.GROMACS_PERT_ITP,
                                           ligname=const.LIGAND_NAME))

            fepl, vdwl, masl, seps = _lambda_paths(self.dummies0,
                                                   self.dummies1)

            with open(VAC_MDP_FILE, 'w') as mdp:
                mdp.write(
                    (VAC_MDP % (COMMON_MDP_TMPL,
                                FE_TMPL)).format(nsteps='2000000',
                                                 seps=seps,
                                                 fep_lambdas=fepl,
                                                 vdw_lambdas=vdwl,
                                                 mass_lambdas=masl))
                                                      
        elif self.FE_sub_type == 'dummy3':
            int_name = topol.int_state._parm_overwrite

            int_state = gromacs.GromacsTop()
            int_state.readParm(int_name + top, int_name + rst)
            int_state.writeTop(int_name + const.GROMACS_ITP_EXT,
                               int_name + '.atp')
            int_state.writeGro(int_name + const.GROMACS_GRO_EXT)

            n1 = 'pert1' + os.extsep
            n2 = 'pert2' + os.extsep
            itp1 = n1 + 'itp'
            atp1 = n1 + 'atp'
            itp2 = n2 + 'itp'
            atp2 = n2 + 'atp'

            gromacs.mixer(top0, int_state, itp1, atp1)
            gromacs.mixer(int_state, top1, itp2, atp2)

            with open(MORPH_TOP1, 'w') as mtop:
                mtop.write(TOP_TMPL.format(title='step 1/2: q_off and vdW '
                                           'on/off',
                                           atp=atp1, itp=itp1,
                                           ligname=const.LIGAND_NAME))

            with open(MORPH_TOP2, 'w') as mtop:
                mtop.write(TOP_TMPL.format(title='step 2/2: q_on',
                                           atp=atp2, itp=itp2,
                                           ligname=const.LIGAND_NAME))

            with open(VAC_MDP_FILE, 'w') as mdp:
                mdp.write(VAC_MDP)
        else:
            raise NotImplementedError

        self.topol = topol


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2):
        """
        Create only topology file but not GRO coordinates.
        """

        self.topol.create_coords(curr_dir, dir_name, lig_morph, pdb_file,
                                 system, cmd1, cmd2)

        if self.FE_sub_type == 'dummy':
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

            fepl, vdwl, masl, seps = _lambda_paths(self.dummies0,
                                                   self.dummies1)
            with open(SOL_MDP_FILE, 'w') as mdp:
                mdp.write(
                    (SOL_MDP % (COMMON_MDP_TMPL,
                                FE_TMPL)).format(nsteps='500000',
                                                 seps=seps,
                                                 fep_lambdas=fepl,
                                                 vdw_lambdas=vdwl,
                                                 mass_lambdas=masl))


TOP_TMPL = '''\
; {title}

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

[ atomtypes ]
#include "{atp}"
#include "{itp}"

[ system ]
FESetup

[ molecules ]
; Compound        nmols
{ligname} 1
'''

COMMON_MDP_TMPL = '''\
integrator               = sd
ld-seed                  = -1
bd-fric                  = 0
dt                       = 0.001
nsteps                   = {nsteps}
nstcomm                  = 100
comm-mode                = Angular

nstxout                  = 10000  ; != nstdhdl (in case of -rerun)
nstvout                  = 0
nstfout                  = 0
nstlog                   = 10000
nstenergy                = 10000
nstxtcout                = 0

tcoupl                   = v-rescale
nsttcouple               = 10
tc_grps                  = System
tau_t                    = 2.0
ref_t                    = 298.0
'''

FE_TMPL = '''\
free-energy              = yes
delta-lambda             = 0
init-lambda-state        = %L%
; lambda paths: {seps}
fep-lambdas              = {fep_lambdas}
vdw-lambdas              = {vdw_lambdas}
mass-lambdas             = {mass_lambdas}
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

VAC_MDP = '''; TI/FEP mdp template for vacuum
%s
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

pcoupl                   = No

gen-vel                  = yes
gen-temp                 = 298.0
gen-seed                 = -1

constraints              = none
continuation             = no

; TI/FEP parameters
%s
'''

SOL_MDP = '''; TI/FEP mdp template for solution
%s
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

; TI/FEP parameters
%s
'''
