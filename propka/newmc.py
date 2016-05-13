r"""
Derived Molecular_container class for better fine control.
"""

__revision__ = "$Id$"



import sys

from propka.molecular_container import *


class Molecular_container_new(Molecular_container):
    """Overwritten calculate_pka() for better fine control."""


    def calculate_pka(self):
        for name in self.conformation_names:
            self.conformations[name].calculate_pka(self.version, self.options)

        self.find_non_covalently_coupled_groups()
        self.average_of_conformations()

        pkas = []

        for residue_type in self.version.parameters.write_out_order:
            for group in self.conformations['AVR'].groups:
                if group.residue_type == residue_type:
                    # FIXME: we assume a protein label here but labels can be
                    #        more complicated, see propka.group
                    resName = group.label[:3].strip()
                    resSeq = int(group.label[3:7])
                    chainID = group.label[7:].strip()

                    pkas.append( (resName, resSeq, chainID, group.pka_value) )

        return pkas
