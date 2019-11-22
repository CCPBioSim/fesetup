===============
Running FESetup
===============

The script FESetup in the release is the command line tool for the end-user.  This shell script sets a few environment variables and eventually calls dGprep.py.  dGprep.py is the actually code running the routines for AFE setup.

Note: this manual describes options used as of Release 1.2.

------------ 
Command line
------------

Default values for all key–value pairs are written to stdout for each of the four sections when FESetup is called without any command-line parameters.  Calling FESetup with '--help' gives information of all possible command line parameters.  Currently all options are for information purposes only so will not affect the setup in any way.  The option --tracebacklimit is really only of use for debugging.  As noted above FESetup is the front-end script to the Python code dGprep.py. ::

    > FESetup --help

    usage: dGprep.py [-h] [-v] [--tracebacklimit N] [infile]

    positional arguments:
     infile  input file in INI format, if not given then just output defaults 

    optional arguments:
     -h, --help          show this help message and exit                
     -v, --version       full version information                       
     --tracebacklimit N  set the Python traceback limit (for debugging) 

------------
Input format
------------

The input file format for the FESetup script (dGprep.py) is an INI like format, popular in the MS-DOS/Windows world.  It is not exactly the same format but a simplified version of it.

The input file may contain four sections where the section names are delimited with brackets:
 
The four sections of the INI file

#. [globals]    global settings, the section name is optional if it is the first section in the file

#. [ligand] 	settings for the ligand

#. [protein] 	settings for the protein

#. [complex] 	settings for the complex

 
Each section consists of various key-value pairs which are two strings separated with an equal sign ("=").  The key must not contain any whitespace.::

    # a typical key-value pair
    morph_pairs = p-aminophenol > o-cresol

Lines may be continued with an initial whitespace (normal space, TAB) on the following line::

    # multiple continuation lines
    morph_pairs = ethane > methanol, ethane > tbutane, ethane > propane,
                  tbutane > propane, tbutane > acetone,
                  propane > acetone, propane > methane

However, list pairs must always appear on the same line because each line is parsed individually.  So the following will cause an error::

    # this will result in an error
    morph_pairs = ethane > methanol, ethane > tbutane, ethane >
                  propane,  # the string "propane" must appear in the previous line!

 
Comments are either empty lines or lines starting with '#' or ';' (leading whitespace is removed).  Inline comments are allowed too provided the comment character is preceded by a space.   Otherwise the string is part of the preceding string, e.g::

    basedir = smallmols # this is an inline comment
    basedir = smallmol#foo # the valid string and directory name 'smallmol#foo'

-------------------------- 
Explicit tagging mechanism
--------------------------

In some cases it may be necessary to explicitly map certain atoms to one another e.g. to preserve a certain spatial rearrangment as for binding modes. This means that the default mapping algorithm (maximum common substructure search, MCSS) can be given additional hints.  In the following example the benzofuran is oriented in the way outlined by forcing a certain mapping.  As should be evident there are twelve ways of mapping the two molecules but if e.g. the binding mode is known a prior the user has here a chance to preserve it!::

    [ligand]
    basedir = smallmols
    morph_pairs = benzol > benzofuran /1=3/2=2  # indices start from 1

An alternative mechanism is to create a special file in basedir (as set in the [ligand] section).  The name of the file must be in the form ligand1~ligand2.map, e.g. if the input reads::

    [ligand]
    basedir = smallmols
    morph_pairs = benzol > benzofuran

then the map file must be in smallmols i.e. smallmols/benzol~benzofuran.map .

The map file has a two column format indicating which atom index in the initial state maps to which atom index in the final state.  Atom indexes start from 1::

    # example mapping file benzol~benzofuran.map in the basedir smallmols/
    # explicitly map the following atom indexes onto each other
    1 3 # this mapping and...
    2 2 # ...this one will fix the orientation of benzofuran in space

