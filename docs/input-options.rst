=============
Input Options
=============

----------------
Full option list
----------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Options unique to each section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following tables list all options unique to each section.  Note that empty strings (denoted as 'none' in the table) means that the user has to use appropriate values. 'molecules' will be overwritten i.e. ignored when 'morph_pairs' are used in [ligand].  The default for complex building is to combine every protein with every ligand.  If you do no want that, you must explicitly list all pairs using the 'pairs' key.  Please note that your file and directory names must not contain the characters ':' (colon), '>' (right angle bracket), '"' (double quote) and '~' (tilde).  The comma ',' is permited as long as the filename is enclosed in double quotes, e.g.

    morph_pairs = "1,2-dichloroethane" > E-dichloroethene,
                  E-dichloroethene > "1,2-dichloroethane"

[globals]

+-----------------------+---------------------+---------+-------------------------------+
| Key                   | Values, default     | Type    |  Explanation                  |
|                       | listed first        |         |                               |  
+=======================+=====================+=========+===============================+
| AFE.type              | Sire, sander/dummy, | string  | free energy type, determines  |
|                       | sander/softcore,    |         | which MD package the input    |
|                       | gromacs,            |         | files are created for (for    |
|                       | pmemd/softcore,     |         | backwards compatability       |
|                       | pmemd/dummy         |         | AMBER = sander/dummy and      |
|                       | charmm/pert         |         | AMBER/softcore =              |
|                       |                     |         | sander/softcore)              |
+-----------------------+---------------------+---------+-------------------------------+
| AFE.separate_vdw_elec | True, False         | bool    | separate the Coulomb (charge) |
|                       |                     |         | transformation from the       |
|                       |                     |         | vdW+bonded transformation     |
+-----------------------+---------------------+---------+-------------------------------+
| forcefield            | amber, ff14SB,      | list of | ff family, subtype of ff,     |
|                       | tip3p, hfe          | strings | water ff, divalent ion set    |
+-----------------------+---------------------+---------+-------------------------------+
| ff_addons             | empty               | list of | additional force fields like  |
|                       |                     | strings | GLYCAM_06j-1 or lipid14       |
+-----------------------+---------------------+---------+-------------------------------+
| gaff                  | gaff1, gaff2        | string  | Choice for the small          |
|                       |                     |         | molecules forcefield, either  |
|                       |                     |         | GAFF 1.x or GAFF2.x           |
+-----------------------+---------------------+---------+-------------------------------+
| logfile               | dGprep.log          | string  | name of the debug log file    |
+-----------------------+---------------------+---------+-------------------------------+
| mdengine              | amber, sander;      | list of | program for minimisation and  |
|                       | amber,pmemd;        | 2       | MD, the first in the list is  |
|                       | gromacs, mdrun      | strings | the MD package, the second is |
|                       | namd, namd2         |         | the actual binary             |
+-----------------------+---------------------+---------+-------------------------------+
| mdengine.prefix       | empty               | string  | the string preceding the      |
|                       |                     |         | mdengine binary command,      |
|                       |                     |         | e.g. mpirun -np 4 (for MPI    |
|                       |                     |         | programs)                     |
+-----------------------+---------------------+---------+-------------------------------+
| mdengine.postfix      | empty               | string  | the string following the      |
|                       |                     |         | mdengine binary command,      |
|                       |                     |         | e.g. +p2 +isomalloc_sync      |
|                       |                     |         | (for namd multicore)          |
+-----------------------+---------------------+---------+-------------------------------+
| parmchk_version       | 2, 1                | integer | parmchk version               |
+-----------------------+---------------------+---------+-------------------------------+
| mcs.timeout           | 60.0                | float   | timeout in seconds for fmcs,  |
|                       |                     |         |    0 means no timeout         |
+-----------------------+---------------------+---------+-------------------------------+
| mcs.match_by          | none, shapealign    | string  | controls RDKit structure      |
|                       | spatially-closest   |         | matching behaviour            |
+-----------------------+---------------------+---------+-------------------------------+
| remake                | False, True         | bool    | remake already done           |
|                       |                     |         | molecules (excluding morphs)  |
+-----------------------+---------------------+---------+-------------------------------+
| overwrite             | False, True         | bool    | by default no files are ever  |
|                       |                     |         | overwritten in the _          |
|                       |                     |         | directories, use this to      |
|                       |                     |         | change this behaviour         |
+-----------------------+---------------------+---------+-------------------------------+
| user_params           | False, True         | bool    | read user force field         |
|                       |                     |         | parameter files, i.e. all     |
|                       |                     |         | .frmod, .preb and .lib (OFF   |
|                       |                     |         | format) files are read in     |
+-----------------------+---------------------+---------+-------------------------------+
 
[ligand]
 	  	  	 
+-------------------------+-----------------------+---------+------------------------------------------------------+
| Key                     | Values, default       | Type    | Explanation                                          |
|                         | listed first          |         |                                                      |
+=========================+=======================+=========+======================================================+
| basedir                 | none,                 | string  | base directory to find ligands                       |
|                         | must be set by user   |         |                                                      |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| file.name               | ligand.pdb            | string  | ligand input file name                               |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| file.format             | none, determined from | string  | format of file.name, can be used to overwrite if     |
|                         | extension of filename |         | file extension is different from actual file format  |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| ions.conc               | 0.0                   | float   | sets the NaCl concentration in mol/l                 |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| ions.dens               | 1.0                   | float   | density for which the ion concentration is wanted    |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| calc_charge             | False, True           | bool    | Force calculation of molecule's formal charge,       |
|                         |                       |         | required e.g. for mol2 format for which Openbabel    |
|                         |                       |         | computes the charge only in select cases.            |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.conj_econv  | 1e-06                 | float   | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.conj_steps  | 250                   | integer | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.ffield      | mmff94                | string  | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.geomsteps   | 5                     | integer | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.numconf     | 0                     | integer | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.steep_econv | 0.0001                | float   | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| conf_search.steep_steps | 100                   | integer | conformation search option                           |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| add_hydrogens           | False, True           | bool    | add hydrogens, allow for adjusting pH                |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| correct_for_pH          | False, True           | bool    | adjust protonation state /partial charges of ligand, |
|                         |                       |         | only works if add_hydrogens = True                   |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| pH                      | 7.4                   | float   | pH for ligand                                        |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| molecules 	          | none,                 | list of | list of molecules                                    | 
|                         | must be set by user   | strings |                                                      |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| morph.absolute          | False, True           | bool    | write absolute transformation MORPH.pert             |
|                         |                       |         | files for Sire                                       |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| morph_pairs             | none,                 | list of | list of pairs in the form lig1 > lig2, overwrites    |
|                         | must be set by user   | strings | 'molecules', do not use '>' in file names            |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| neutralize              | False, True           | bool    | neutralize the solvation box by adding minimum       |
|                         |                       |         | couterions required                                  |
+-------------------------+-----------------------+---------+------------------------------------------------------+
| skip_param              | False, True           | bool    | skip the parameterisation step, useful in            |
|                         |                       |         | conjunction with user_params or ff_addons            |
|                         |                       |         | see [globals]                                        |
+-------------------------+-----------------------+---------+------------------------------------------------------+
  	  	  	 
[protein] 	  	  	 

+----------------------+--------------------------+---------+------------------------------------------------------+
| Key                  | Values, default          | Type    | Explanation                                          |
|                      | listed first             |         |                                                      |
+======================+==========================+=========+======================================================+
| align_axes           | False, True              | bool    | align protein along principal axes before hydrating  |
+----------------------+--------------------------+---------+------------------------------------------------------+
| basedir              | none,                    | string  | base directory to proteins                           |
|                      | must be set by user      |         |                                                      |
+----------------------+--------------------------+---------+------------------------------------------------------+
| ions.conc            | 0.0                      | float   | sets the NaCl concentration in mol/l                 |
+----------------------+--------------------------+---------+------------------------------------------------------+
| ions.dens            | 1.0                      | float   | density for which the ion concentration is wanted    |
+----------------------+--------------------------+---------+------------------------------------------------------+
| molecules            | none,                    | list of | list of molecules                                    |
|                      | must be set by user      | strings |                                                      |
+----------------------+--------------------------+---------+------------------------------------------------------+
| neutralize           | False, True              | bool    | neutralize the solvation box by adding minimum       |
|                      |                          |         | couterions required                                  |
+----------------------+--------------------------+---------+------------------------------------------------------+
| propka               | False, True              | bool    | use ProPKA to protonate protein                      |
+----------------------+--------------------------+---------+------------------------------------------------------+
| propka.pH            | 7.0                      | float   | pH for ProPKA                                        |
+----------------------+--------------------------+---------+------------------------------------------------------+
  	  	  	 
[complex] 	  	  	 

+----------------------+--------------------------+---------+------------------------------------------------------+
| Key                  | Values, default          | Type    | Explanation                                          |
|                      | listed first             |         |                                                      |
+======================+==========================+=========+======================================================+
| align_axes           | False, True              | bool    | align protein along principal axes before hydrating  |
+----------------------+--------------------------+---------+------------------------------------------------------+
| ions.conc            | 0.0                      | float   | sets the NaCl concentration in mol/l                 |
+----------------------+--------------------------+---------+------------------------------------------------------+
| ions.dens            | 1.0                      | float   | density for which the ion concentration is wanted    |
+----------------------+--------------------------+---------+------------------------------------------------------+
| neutralize           | False, True              | bool    | neutralize the solvation box by adding minimum       |
|                      |                          |         | couterions required                                  |
+----------------------+--------------------------+---------+------------------------------------------------------+
| flatten_rings        | False, True              | bool    | make aromatic rings fully planar, for MC with Sire   |
+----------------------+--------------------------+---------+------------------------------------------------------+
| pairs                | none, must be set by     | list of | list of pairs in the form protein:ligand, do not use |
|                      | user                     | strings |  ':' in file names                                   |
+----------------------+--------------------------+---------+------------------------------------------------------+
  	  	  	 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The minimisation and MD options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following options are the minimsation and MD options for molecule setup common to ligands, proteins and complexes.  To allow minimisation and MD 'box.type' has to be set explicilty which also creates a water box.  If 'box.type' is not set by the user then no box will be created and minimisation or MD will not be carried out.  To actually run a minimsation or simulation  you will need to set any of the '.nsteps' keys to a value larger than 0.  The only difference is relaxation where setting 'md.relax.nrestr' to a value larger than 0 will trigger restraint relaxation.  The order of simulation protocols is fixed as heating (md.heat.*), constant volume and temperature (md.constT.*), pressurising = density adjustment (md.press.*), relaxation at NpT conditions (md.relax.*).  If any of those steps are not needed set '.nsteps' to 0 but be aware that there are no further sanity checks.  The MD protocol can be preceded by a minimisation step (min.*).

 
+-----------------+------------------------------+---------+--------------------------------------------------------+
| Key             | Values, default listed first | Type    | Explanation                                            |
+=================+==============================+=========+========================================================+
| box.type        | empty string = no box        | string  | creates a box of water                                 |
|                 | created, rectangular,        |         |                                                        |
|                 | octahedron (limited support) |         |                                                        |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| box.length      | 10.0                         | float   | the distance in Ångström between solute and the box    |
|                 |                              |         | edges, NOTE: the TIP3P box will create a system of low |
|                 |                              |         | density and thus this distance will decrease on        |
|                 |                              |         | pressuring the sytem.                                  |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| min.ncyc        | 0                            | integer | number of steepest decent steps in minimisation        |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| .nsteps         | 0                            | integer | number of steps; e.g. min.nsteps                       |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| .restr_force    | 10.0                         | float   | restraint force; e.g. md.heat.restr_force              |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| .restraint      | protein, backbone, heavy,    | string  | restraint type, if other string then in the list it is |
|                 | notligand, notsolvent        |         | the restraintmask for sander; e.g. md.constT.restraint |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| .T              | 300.0                        | float   | temperature; e.g. md.press.T                           |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| .p              | 1.0                          | float   | pressure; e.g. md.relax.p                              |
+-----------------+------------------------------+---------+--------------------------------------------------------+
| md.relax.nrestr | 0                            | integer | number of relaxation steps, needed to trigger          |
|                 |                              |         | restraint relaxation                                   |
+-----------------+------------------------------+---------+--------------------------------------------------------+
