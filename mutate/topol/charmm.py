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

r"""
Create perturbed topologies for CHARMM/PERT.
"""


__revision__ = "$Id$"


import os, sys
from collections import defaultdict

import Sire.MM
import Sire.Units

import sander
from FESetup import const, errors, logger
from FESetup.prepare.amber import charmm



VAC_INP_FILE = 'vac.inp'
SOL_INP_FILE = 'sol.inp'
STATE_INT = 'state_int'

VAC_TEMPLATE = '''\
* relative AFE in vacuo
* Please note that this input file is only intended as a rough template to
* illustrate how a TI simulation would need to be set up
*

bomlev -2
prnlev 5

read rtf card name combined.rtf
read param flex card name combined.prm

ioformat extended
read psf xplor card name "{state0}.psf"
read coor card name "{state0}.cor"

use amber

! {charge0}

format (i1)
calc old = @step - 1
set oldrst fe@old.rst

format (f15.7)

if @lambda .lt. 0 then set lambda 0
if @lambda .gt. 1 then set lambda 1

set lstart @lambda
set lstop @lambda

calc delta = 0.5 / (@nlambda - 1)
calc lstart = @lstart - @delta
calc lstop  = @lstop  + @delta

if @lambda .eq. 0 then set lstart 0.0
if @lambda .eq. 1 then set lstop 1.0

lower
set base fe@step
set rst @base.rst
set dcd @base.dcd
set en  @base.en
set whamf @base.wham

open unit 10 write form name @rst
open unit 11 write unform name @dcd
open unit 12 write form name @en
open unit 22 read form name @oldrst
open unit 54 write card name @whamf

nbonds atom cdie shift vatom vswith cutnb 999.0 ctofnb 997.0
update inbfrq -1 ihbfrq -1

pert select segid A end
delete atoms select all end

read psf xplor card name "{state1}.psf"
read coor card name "{state1}.cor"

! {charge1}

prnlev 5

scalar fbeta set 5.0 select all end

dynamics langevin leap start nstep 500000 timestep 0.001 -
  lstart @lstart lambda @lambda lstop @lstop wham 54 -
    pstart 0 pstop 500000 pwindow nopssp -
  tbath 300.0 rbuf 0.0 ilbfrq 10 -
  iunwri 10 iuncrd 11 kunit 12 iunrea 22 -
  nsavc 10000 nprint 5000 iprfreq 20000 isvfrq 10000 ntrfrq 1000

stop
'''

TR_TABLE = {'pert': 'dummy', 'pert2': 'dummy2', 'pert3': 'dummy3'}


class PertTopology(object):

    def __init__(self, FE_sub_type, separate, ff, con_morph, atoms_initial,
                 atoms_final, lig_initial, lig_final, atom_map,
                 reverse_atom_map, zz_atoms, gaff):

        if FE_sub_type[:4] != 'pert':
            raise errors.SetupError('Morph code only supports the '
                                    'CHARMM/PERT module')
        self.separate = separate
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
        self.itypes = []
        self.ftypes = []

        self.files_created = []
        
        self.dummies0 = not all([a.atom for a in self.atom_map.keys()])
        self.dummies1 = not all([a.atom for a in self.atom_map.values()])

        if self.separate and self.dummies0 and self.dummies1:
            self.FE_sub_type = 'pert3'
        elif self.separate and (self.dummies0 or self.dummies1):
            self.FE_sub_type = 'pert2'
        else:
            self.FE_sub_type = 'pert'

            if self.separate:
                logger.write('Warning: linear transformation, not separated '
                             'into vdw and electrostatic step\n')

        self.stype = self.FE_sub_type


    def setup(self, curr_dir, lig_morph, cmd1, cmd2):

        for atom in lig_morph.atoms():
            self.itypes.append(str(atom.property('initial_ambertype') ) )

        for iinfo, finfo in self.atom_map.items():
            if not finfo.atom:
                final_type = "du"
            else:
                base = self.atoms_final.select(finfo.index)
                final_type = "%s" % str(base.property("ambertype") )

            self.ftypes.append(final_type)

  
        topol = sander.PertTopology('_' + TR_TABLE[self.stype], self.separate,
                                    self.ff, self.con_morph, self.atoms_initial,
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
        top0 = charmm.CharmmTop(self.ftypes)
        top0.readParm(lig0 + top, lig0 + rst)
        top0.writePsf('state0.psf')
        top0.writeCrd('state0.cor')
        self.files_created.extend(('state0.psf', 'state0.cor'))

        # NOTE: pert2 - dummies zero q, vdW
        #       pert3 - disappearing zero q
        if self.stype == 'pert2':
            top_int = charmm.CharmmTop()
            top_int.readParm(STATE_INT + top, STATE_INT + rst)
            top_int.writePsf(STATE_INT + '.psf')
            top_int.writeCrd(STATE_INT + '.cor')
            self.files_created.extend((STATE_INT + '.psf', STATE_INT + '.cor'))

        # for pert3 we go down the sander route otherwise we would need
        # 4 topologies (0, 0', 1', 1 where ' means q_on/off only:
        # step1: state0 + state0(scalar charge set)
        # step2: state0(scalar charge set) + state1(scalar charge set)
        # step3: state1(scalar charge set) + state1

        top1 = charmm.CharmmTop(self.itypes)
        top1.readParm(lig1 + top, lig1 + rst)
        top1.writePsf('state1.psf')
        top1.writeCrd('state1.cor')

        top0.combine(top1)              # adds parms from top1 to top0!
        top0.writeRtfPrm('combined.rtf', 'combined.prm')

        self.files_created.extend(('state1.psf', 'state1.cor',
                                   'combined.rtf', 'combined.prm'))

        self.topol = topol

        with open(VAC_INP_FILE, 'w') as inp:
            inp.write(VAC_INP)


    def create_coords(self, curr_dir, dir_name, lig_morph, pdb_file, system,
                      cmd1, cmd2, boxdims):
        """
        Create only topology file, not coordinates.
        """

        # FIXME: support FE_sub_type
        self.topol.create_coords(curr_dir, dir_name, lig_morph, pdb_file,
                                 system, cmd1, cmd2, boxdims)

        lig0 = self.topol.lig0._parm_overwrite
        lig1 = self.topol.lig1._parm_overwrite
        top = self.topol.lig0.TOP_EXT
        rst = self.topol.lig0.RST_EXT

        # top0 and top1 written mainly for debugging purposes
        top0 = charmm.CharmmTop(self.ftypes)
        top0.readParm(lig0 + top, lig0 + rst)
        top0.writePsf('state0.psf')
        top0.writeCrd('state0.cor')

        if self.stype == 'pert2' or self.stype == 'pert3':
            top_int = charmm.CharmmTop()
            top_int.readParm(STATE_INT + top, STATE_INT + rst)
            top_int.writePsf(STATE_INT + '.psf')
            top_int.writeCrd(STATE_INT + '.cor')

        top1 = charmm.CharmmTop(self.itypes)
        top1.readParm(lig1 + top, lig1 + rst)
        top1.writePsf('state1.psf')
        top1.writeCrd('state1.cor')

        top0.combine(top1)              # adds parms from top1 to top0!
        top0.writeRtfPrm('combined.rtf', 'combined.prm')

        with open(SOL_INP_FILE, 'w') as inp:
            inp.write(SOL_INP)



VAC_INP = '''* vacuum
* FESetup
*

'''

SOL_INP = '''* solution
* FESEtup
*

'''
