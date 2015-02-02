#  Copyright (C) 2012-2014  Hannes H Loeffler
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
Utility functions for the amber package.
"""


__revision__ = "$Id$"


import sys, os, re, shlex, string
import subprocess as subp

from FESetup import const, errors, logger



def self_check():
    """
    Check if AMBER is properly set up by checking for bin or exe directory in
    AMBERHOME

    :returns: error message or None for OK
    """

    if 'AMBERHOME' in os.environ:
        const.AMBER_BIN_PATH = os.path.join(os.environ['AMBERHOME'], 'bin')

        if not os.access(const.AMBER_BIN_PATH, os.X_OK):
            const.AMBER_BIN_PATH = os.path.join(os.environ['AMBERHOME'], 'exe')

            if not os.access(const.AMBER_BIN_PATH, os.X_OK):
                return ('AMBERHOME is not set properly: AMBERHOME = %s'
                        % os.environ['AMBERHOME'])
            else:
                return ('AMBERHOME is not set')

    return None


def _cleanup_string(in_str):
    """
    Check if in_str is 'empty'.  Internal function only.

    Returns empty string if in_str contains only white space.
    Returns in_str with beginning and ending new lines stripped and each
    line preceded with four spaces.
    """


    for char in in_str:
        if char not in string.whitespace:
            text = re.sub('^\n+|\n+$', '', in_str)
            text = re.sub('\n', '\n    ', text)

            return text

    return ''


def check_amber(program):
    """
    Check the AMBER executable path for the existence of program.

    :param program: AMBER executable file name
    :type program: string
    :returns: full path of the program
    """


    exe = os.path.join(const.AMBER_BIN_PATH, program)

    if not os.access(exe, os.X_OK):
        sys.exit('Cannot run %s' % exe)

    return exe


def run_amber(program, params):
    """
    Simple wrapper to execute external AMBER programs through subprocess.

    :param program: AMBER program file name
    :type program: string
    :param params: paramters to the AMBER program
    :type params: string
    :raises: SetupError
    :returns: True on failure
    """


    cmd = shlex.split(program)
    cmd.extend(shlex.split(params))

    logger.write('Executing command:\n%s %s\n' % (program, params) )

    proc = subp.Popen(cmd, stdout = subp.PIPE, stderr = subp.PIPE)
    out, err = proc.communicate()

    for stream in out, err:
        text = _cleanup_string(stream)

        if text:
            logger.write('  %s {\n    %s\n  }\n' %
                         ('stdout' if stream is out else 'stderr', text) )

    if proc.returncode:
        return out, err

    return False


def run_leap(top, crd, program = 'tleap', script = 'leap.in'):
    """
    Simple wrapper to execute the AMBER leap program.

    :param top: topology file name, used to check if created
    :type top: string
    :param crd: coordinate file name, used to check if created
    :type crd: string
    :param program: program file name, either tleap or sleap (not recommended)
    :type program: string
    :param script: leap script as string, if 'leap.in' read from respective file
      name
    :type script: string
    :returns: output from leap
    :raises: SetupError
    """


    leap = check_amber(program)
    cmd = [leap, '-f']

    if script == 'leap.in':
        cmd.append(script)
        logger.write('Executing command:\n%s' % ' '.join(cmd) )

        proc = subp.Popen(cmd, stdin = None, stdout = subp.PIPE,
                          stderr = subp.PIPE)
        out = proc.communicate()[0]
    else:
        cmd.append('-')
        logger.write('Executing command:\n%s -f - <<_EOF \n%s\n_EOF\n' %
                     (leap, script) )

        proc = subp.Popen(cmd, stdin = subp.PIPE, stdout = subp.PIPE,
                          stderr = subp.PIPE)
        out = proc.communicate(script)[0]

    if top and crd:
        if os.path.getsize(top) == 0 or os.path.getsize(crd) == 0:
            raise errors.SetupError(
                'Leap did not create the topology and/or coordinate '
                'file(s): %s, %s' % (top, crd)
                )

    return out


def run_ambpdb(params, infile, outfile):
    """
    Simple wrapper to execute the AMBER ambpdb module.

    :param params: input parameters to ambpdb
    :type param: string
    :param infile: input file name
    :type infile: string
    :param outfile: output file name
    :type outfile: string
    :raises: SetupError
    """


    program = check_amber('ambpdb')
    cmd = [program]
    cmd.extend(shlex.split(params))

    logger.write('Executing command: ambpdb %s < %s > %s\n' %
                 (params, infile, outfile) )

    with open(infile, 'r') as inp:
        with open(outfile, 'w') as out:
            proc = subp.Popen(cmd, stdin = inp, stdout = out,
                              stderr = subp.PIPE)
            err =  proc.communicate()[1]

    text = _cleanup_string(err)

    if text:
        logger.write('  stderr {\n    %s\n  }\n' % text)

    if proc.returncode:
        raise errors.SetupError('%s failed' % program)


def run_namd(program, prefix, postfix, params, std_out):
    """
    Simple wrapper to execute the external NAMD program through subprocess.

    :param program: executable name
    :type program: string
    :param std_out: file name for captured NAMD output from stdout
    :type std_out: string
    :raises: SetupError
    """

    exe = os.path.join(os.environ['NAMDHOME'], program)

    if not os.access(exe, os.X_OK):
        sys.exit('Cannot run %s' % exe)

    cmdline = ' '.join( (prefix, exe, postfix, params) ) 
    cmd = shlex.split(cmdline)

    logger.write('Executing command:\n%s\n' % cmdline)

    proc = subp.Popen(cmd, stdout = subp.PIPE, stderr = subp.PIPE)
    out, err =  proc.communicate()

    with open(std_out, 'w') as text:
        text.writelines(out)

    if proc.returncode:
        return err

    return None


def run_gromacs(program, prefix, postfix, params):
    """
    Simple wrapper to execute the external GROMACS program through subprocess.

    :param program: executable name
    :type program: string
    :param params: paramters to the AMBER program
    :type params: string
    :raises: SetupError
    """

    exe = os.path.join(os.environ['GMXHOME'], 'bin', program)

    if not os.access(exe, os.X_OK):
        sys.exit('Cannot run %s' % exe)

    cmdline = ' '.join( (prefix, exe, postfix, params) ) 
    cmd = shlex.split(cmdline)

    logger.write('Executing command:\n%s\n' % cmdline)

    proc = subp.Popen(cmd, stdout = subp.PIPE, stderr = subp.PIPE)
    out, err =  proc.communicate()

    if proc.returncode:
        return 1, err

    return out, err


def run_dlpoly(program, prefix, postfix):
    """
    Simple wrapper to execute the external DL_POLY program through subprocess.

    :param program: executable name
    :type program: string
    :raises: SetupError
    """

    exe = os.path.join(os.environ['DLPOLYHOME'], 'execute', program)

    if not os.access(exe, os.X_OK):
        sys.exit('Cannot run %s' % exe)

    cmdline = ' '.join( (prefix, exe, postfix) ) 
    cmd = shlex.split(cmdline)

    logger.write('Executing command:\n%s\n' % cmdline)

    proc = subp.Popen(cmd, stdout = subp.PIPE, stderr = subp.PIPE)
    out, err =  proc.communicate()

    if proc.returncode:
        return 1, err

    return out, err
