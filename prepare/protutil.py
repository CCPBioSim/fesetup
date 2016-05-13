#  Copyright (C) 2012-2013  Hannes H Loeffler
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
Common methods for the force field protein classes.  These methods are expected
to be directly imported into the Protein class.
"""

__revision__ = "$Id$"



import os

from FESetup import const, report, CaptureOutput, logger


# NOTE: PROPKA 3.1 can actually also protonate ligands but for the time being
#       we only allow protein protonation
@report
def protonate_propka(self, pH = 7.0):
    """
    Protonate a protein with PROPKA 3.1.

    :param pH: desired pH for protein protonation
    :type pH: float
    """

    import os
    import sys
    import  StringIO

    import propka.lib as plib
    import FESetup.propka.newmc as pmc

    PROT_RES = ('HIS', 'ASP', 'GLU')
    DEPROT_RES = ('LYS', 'CYS', 'ARG', 'TYR')

    # add PDB file name to options to avoid warning of missing file, also
    # set config file path to module path to find propka.cfg
    options, dummy = plib.loadOptions( ['--pH', pH, '-q', self.mol_file] )
    options.parameters = os.path.join(os.path.dirname(pmc.__file__),
                                      options.parameters)

    with CaptureOutput() as output:
        mol = pmc.Molecular_container_new(self.mol_file, options)
        pKas = mol.calculate_pka()

    logger.write('%s%s' % (output[0], output[1]) )

    protres = []

    for resName, resSeq, chainID, pKa in pKas:
        # NOTE: currently we ignore termini 'N+' and 'C-'
        if pH < pKa:
            if resName in PROT_RES:
                protres.append( (resName, resSeq, chainID) )
        else:
            if resName in DEPROT_RES:
                protres.append( (resName, resSeq, chainID) )

    logger.write('pH = %.2f' % pH)

    msg_res = set()

    with open(const.PROTONATED_PDB_FILE, 'w') as newfile:
        with open(self.mol_file, 'r') as pdbfile:
            for line in pdbfile:
                if line[:6] in ('ATOM  ', 'HETATM'):
                    resName = line[17:21].strip()
                    resSeq = int(line[22:26])
                    chainID = line[20:22].strip()

                    if (resName, resSeq, chainID) in protres:
                        msg_res.add( (resName, resSeq, chainID) )
                        newfile.write(line[:17] +
                                      '{:3s} '.format(self.PROT_MAP[resName]) +
                                      line[21:])
                    else:
                        newfile.write(line)

    for resName, resSeq, chainID in msg_res:
        logger.write('Changing %s %i %s to %s' %
                          (resName, resSeq, chainID, self.PROT_MAP[resName]) )

    self.mol_file = const.PROTONATED_PDB_FILE
