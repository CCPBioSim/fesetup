# FESetup

FESetup automates the setup of relative alchemical free energy (AFE) simulations such as thermodynamic integration (TI) and free energy perturbation (FEP).  Post–processing methods like MM–PBSA and LIE are supported as well.  FESetup can also be used for general simulation setup ("equilibration") through an abstract MD engine (currently supported MD engines are AMBER, GROMACS, NAMD and DL_POLY).  For relative AFE simulation the mapping of corresponding atoms between the two free energy states, that is their topological similarity, is computed via a maximum common substructure search (MCSS).  This enables a maximal single topology description of the perturbed molecule pair.  Ligand molecules can automatically be parameterised using the AMBER GAFF/AM1-BCC method.  Supported force fields for biomolecules are all the modern AMBER force fields.

The AFE simulation packages that are currently supported are Sire, AMBER, GROMACS and CHARMM/PERT.  All these codes implement AFE simulation by making use of a hybrid single/dual topology description of the perturbed region i.e. the mapped region (single topology) can be used simultaneously with an un–mapped, duplicated region (dual topology) from each state. There is also some support for NAMD's purely dual topology implementation but this requires an additional PDB file to mark appearing/vanishing atoms and possibly relative restraints to keep ligands spatially in place and/or together.

FESetup particularly aims at automation where it makes sense and is possible, ease of use and robustness of the code.  Users are very welcome to discuss on our forum, report issues and request new features.  The software is licensed under the GPL2 and such is a community effort: user contributions in any form are highly encouraged!

The basis of the current code was a collection of Python and shell scripts written previously by Julien Michel and Christopher Woods. The current code base is now mainly developed by Hannes Loeffler (STFC) with contributions from the original developers.

Please cite [DOI: 10.1021/acs.jcim.5b00368](https://dx.doi.org/10.1021/acs.jcim.5b00368) when you use FESetup.
