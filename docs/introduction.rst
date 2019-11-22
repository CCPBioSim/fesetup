============
Introduction
============ 

FESetup automates the setup of relative alchemical free energy (AFE) simulations such as thermodynamic integration (TI) and free energy perturbation (FEP).  Post–processing methods like MM–PBSA and LIE are supported as well.  FESetup can also be used for general simulation setup ("equilibration") through an abstract MD engine (currently supported MD engines are AMBER, GROMACS, NAMD and DL_POLY).  For relative AFE simulation the mapping of corresponding atoms between the two free energy states, that is their topological similarity, is computed via a maximum common substructure search (MCSS).  This enables a maximal single topology description of the perturbed molecule pair.  Ligand molecules can automatically be parameterised using the AMBER GAFF/AM1-BCC method.  Supported force fields for biomolecules are all the modern AMBER force fields.

The AFE simulation packages that are currently supported are Sire, AMBER, GROMACS and CHARMM/PERT.  All these codes implement AFE simulation by making use of a hybrid single/dual topology description of the perturbed region i.e. the mapped region (single topology) can be used simultaneously with an un–mapped, duplicated region (dual topology) from each state. There is also some support for NAMD's purely dual topology implementation but this requires an additional PDB file to mark appearing/vanishing atoms and possibly relative restraints to keep ligands spatially in place and/or together.

FESetup particularly aims at automation where it makes sense and is possible, ease of use and robustness of the code.  Users are very welcome to discuss on our forum, report issues and request new features.  The software is licensed under the GPL2 and such is a community effort: user contributions in any form are highly encouraged!

The basis of the current code was a collection of Python and shell scripts written previously by Julien Michel and Christopher Woods. The FESetup1.2 code base was mainly developed by Hannes Loeffler (STFC) with contributions from the original developers.

Please cite DOI: 10.1021/acs.jcim.5b00368 when you use FESetup.

----------------------------
Capabilities and Description
----------------------------

The sections here give a brief overview of FESetup from a general point of view i.e. topics like background information, concepts, implementation ideas and philosophies, etc. but also the limits of the current software.

------- 
General
-------

FESetup is currently (July 2017) mainly geared at automated setup for receptor-ligand simulations.  But in principle any biomolecular system can be set up for simulation provided the system is supported by the force field (AMBER and AMBER/GAFF) and the ligand is not covalently bound and can reasonably be parameterised with AM1/BCC.  The latter implies that the ligand is a relatively "small" organic molecule as typically considered in the drug development process.

Core functionality of FESetup are the parameterisation of (possibly large numbers of) ligands and the "morphing" of a pair of ligands.  Morphing describes how a molecule's force field parameters are transformed (or "perturbed") into a set of parameters descriptive for a second, topological related molecule.  As such FESetup is a perturbed topology creator.  Ligand parameterisation enables standard MD simulation including post-processing methods like MM-PBSA for simplified free energy estimates.  The topology file describing the morphing of the ligand pair enables alchemical free energy simulation like Thermodynamic Integration (TI) or Free Energy Perturbation (FEP, also called EXP for exponential formula).  Currently, FESetup supports such simulation setups for Sire, AMBER, Gromacs and the PERT module in CHARMM.  In the case of AMBER the dual topology softcore approach as well as the older explicit dummy atom approach, as implemented in both sander and pmemd, are supported.

FESetup reads a series of coordinates for ligands and proteins.  The ligands can then be parameterised with the AM1/BCC method and can be combined with a protein (or other molecule supported by the AMBER force fields or with user supplied parameters) into a complex.  This requires the complex to be set up properly beforehand, e.g. the complex may have been created via a previous docking run or by any other suitable method.  Protein and ligand coordinates must be supplied separately.  While the protein is expected in PDB format the ligand can be in any format supported by OpenBabel (but must have 3D coordinates).  A caveat though is that OpenBabel does not read total charge information from all file formats that support this and thus the most suitable formats are either PDB or SDF for the ligand if a charge other than zero is needed.  Internally the coordinates will be converted to the mol2 format.  If morph pairs are requested, FESetup will compute the "difference" between the two ligand molecules using a maximum common substructure search (MCSS) algorithm (from RDKit) as a distance metric and use this information to create appropriate topology files and control files for alchemical free energy simulation.

---------- 
MD engines
----------

FESetup makes use of an abstract MD engine, currently for the purpose of "equilibration" of simulation systems.  This means that various popular MD simulation software packages (at the moment: AMBER, DL_POLY, NAMD and GROMACS) can be used in a transparent fashion, that is without the need to know the specific control and command structures of a particular simulation software.  This is, of course, somewhat limited by MD software not always providing fully equivalent features as in other software packages.  But for the general purpose of carrying out setup processes like minimization, heating, pressurizing, restraint release, etc. this is sufficient and sensible starting structures for production MD simulation can be obtained without problems.

The MD engine mechanism allows the user to choose the specific binary (program) to be used for MD/minimisation, e.g. with AMBER that could be any variant of sander or pmemd, with GROMACS any variant of mdrun, etc.  Parallel versions are supported to a certain extent too as described in the following.  Multi-threaded binaries like mdrun with thread-MPI or the multi-core version of NAMD will run on any single machine or "node".  MPI versions like the AMBER binaries do also work but currently there is no support for job schedulers like PBS, LSF, etc.

Some of the limitations of abstraction are that software like GROMACS or NAMD appear to always need a minimization step of a solvated system prior to MD whereas AMBER is less sensitive to this.  FESetup currently includes standard input configuration settings for minimization, random velocity assignment, heating, constant T and pressurizing.  For all these, restraints can be defined by either keyword or an AMBER mask.  All corresponding run-time control parameters are created through templates and each of these steps is carried out as an individual simulation run.  This is also done for NAMD despite it having full scripting capabilities of its own.  However, the abstract interface hides the details of the actual MD engine away and thus is much easier to use.  In addition, the interface is unified i.e. it appears the same independent of the MD engine chosen.  AMBER and GROMACS do no have built-in scripting facilities and thus the separate steps are needed especially the step-wise release of restraints.

----------------------------------
Maximum Common Substructure Search
----------------------------------

Morph pairs are generated through a graph based maximum common substructure search (MCSS) algorithm, in particular the implementation in RDKit.  By comparing the graphs of two ligands the maximum match between their atoms, based on the topology i.e. connectivity, are found.  Currently we do not distinguish between atom types or bond types, i.e. the graph is unlabeled.  Rings, on the other hand, are required to match other rings and cannot be broken.  All atoms not part of the MCS mapping (equivalent atoms between the two ligands are recorded), are considered to be dummy atoms and coordinates are created accordingly from internal coordinates (this is not needed in the AMBER dual topology softcore apprach).

There are, however, a few caveats to keep in mind.  First, the comparison is carried out in graph space which is essentially only two dimensional.  Thus 3D features are not retained unless otherwise specified.  This includes the substitution patterns of stereogenic centres, e.g. the absolute configuration can be reverted with the scheme described above.  Symmetry related issues, e.g. molecules may have multiple MCSs but a "random" one is assigned to a given morph pair, may appear too.   Retaining specific binding modes is another problem which may be caused by symmetry but other origins, e.g. flipping by 180 degrees or parts of a molecules preferring another binding pocket when decorated differently, are possible too.  FESetup allows the user to provide explicit atom tagging of individual atoms to overcome these problems albeit at the expense of automation.  Tagging however allows the user to overwrite the default behaviour of FESetup and can thus guide the mapping to their own preferences.

Second, the time behaviour of MCSS algorithms may pose problems.  In the worst case a search would be exponential in time due to the need of an exhaustive search (in practice, the MCSS algorithms have clever shortcuts but there is no universal algorithm available because NP-complete) which implies that for every additional atom the search time would be doubled.  We have found several examples where the MCS search can range from many minutes to hours to days (in a very large system of a molecular weight of nearly 1500 with the old Python implementation of fmcs).  In practice, however, the MCS appears to be found within just a minute or so (more careful testing needed) and we have not found any problem cases yet where the MCS would not be as expected from visual inspection.  FESetup  provides a setting to limit the time spend on the MCS search.

----------------
Adding hydrogens
----------------

This is a particular serious matter and attitudes among researchers vary.  We strongly recommend that the user makes sure that hydrogens are added to the ligand and protonation state as well as tautomeric state are fully determined before FESetup is run.  It is the user's responsibility to get the chemistry right.

In principle, there is a setting in FESetup to allow addition of hydrogens in a simple valence filling fashion.  But we have found that this does not always work properly with OpenBabel, e.g. the N9 (binding to the ribose) in nucleosides appears to be always perceived as being located within a double-bonded or aromatic ring.  This results in addition of a hydrogen and a charge of +1 on N9.  However, charge parameterisation is achieved through sqm, a semi-empirical tool in the AmberTools, which crucially depends on having the charge properly assigned, or otherwise may terminate with an error or compute grossly wrong results.

Finding protonation and tautomeric states is even much more difficult because they depend strongly on the environment.  There is currently no support for assignment of such states in FESetup.  Any future work in this direction will have to await very thorough investigations into what is possible and what not.

----------------------- 
Charge parameterisation
-----------------------

The parameterisation steps in FESetup currently carries out charge calculation for the AMBER/GAFF force field at the AM1/BCC level with the help of the AmberTools toolchain antechamber.  This method derives Mulliken charges at the semiempirical AM1 level of theory (via sqm) and then applies bond charge corrections (BCC) to these Mulliken charges to finally obtain charges almost equivalent to the higher level HF 6-31G* RESP charge derivation scheme.  A subsequent run through parmchk(2) will assign missing bonded parameters through a similarity search in the first attempt or, if this cannot be achieved, an empirical approach is followed.  If both methods fail the correpsonding force field parameters will be set to zero.  This is actually exploited in the case of the dummy atoms as the corresponding atom type does not have an equivalent in the data base (the code relies on zero force field parameters for identifying dummy atoms at various places).  Lennard-Jones parameters and mixing rules will be applied as per the initial perception of GAFF atom types (atomtype and bondtype tools).

There are a few notable limitations.  Sqm carries out geometry minimizations in vacuum.  This can lead to distorted structures when highly charged groups are present, e.g. zwitterions. We have also found proton "shifts" between the two ends of a molecule and "decarboxylation", both in very rare cases.  FESetup provides a (undocumented) option to carry out geometry minimizations with the help of an implicit solvent model (through the QM/MM feature in sander without an actual MM part).  But we did not really find an improvement for the aforementioned problems except for more stable SCF convergence and a smoother geometry optimisation.

Larger molecules should probably be considered to be broken into sensible smaller fragments similar to how the biomolecules force field has been parameterised.  This is beyond the scope of FESetup currently, however.  More elaborate schemes, on the other hand. like the original RESP method or possibly specialised software like R.E.D. may be incorporated in the future.

