#  Copyright (C) 2012-2014,2016  Hannes H Loeffler
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


import sys
import os
import re
import shlex
import string
import glob
import subprocess as subp

from FESetup import const, errors, logger



def self_check():
    """
    Check if AMBER is properly set up by checking for bin or exe directory in
    AMBERHOME

    :returns: error message or None for OK
    :rtype: str
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


def _setenv():
    """
    Check if we are running on the internal or an external AMBER.  If it
    is external, make sure that it doesn't use LD_LIBRARY_PATH but a copy
    if it exists.  Otherwise, use the current LD_LIBRARY_PATH.
    """
    env = os.environ.copy()

    if 'LOCAL_AMBER' in os.environ and \
           os.environ['LOCAL_AMBER'] != os.environ['AMBERHOME'] and \
           'COPY_LD_LIBRARY_PATH' in os.environ:
        env['LD_LIBRARY_PATH'] = env['COPY_LD_LIBRARY_PATH']

    return env


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

    env = _setenv()
    proc = subp.Popen(cmd, stdout=subp.PIPE, stderr=subp.PIPE, env=env)
    out, err = proc.communicate()

    for stream in out, err:
        text = _cleanup_string(stream)

        if text:
            logger.write('  %s {\n    %s\n  }\n' %
                         ('stdout' if stream is out else 'stderr', text) )

    if proc.returncode:
        return out, err

    return False


def run_leap(top, crd, program='tleap', script=''):
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

    env = _setenv()

    if script == 'leap.in':
        cmd.append(script)
        logger.write('Executing command:\n%s' % ' '.join(cmd) )

        proc = subp.Popen(cmd, stdin=None, stdout=subp.PIPE, stderr=subp.PIPE,
                          env=env)
        out = proc.communicate()[0]
    else:
        cmd.append('-')
        logger.write('Executing command:\n%s -f - <<_EOF \n%s\n_EOF\n' %
                     (leap, script) )

        proc = subp.Popen(cmd, stdin=subp.PIPE, stdout=subp.PIPE,
                          stderr=subp.PIPE, env=env)
        out = proc.communicate(script)[0]

    if top and crd:
        if os.path.getsize(top) == 0 or os.path.getsize(crd) == 0:
            raise errors.SetupError(
                'Leap did not create the topology and/or coordinate '
                'file(s): %s, %s' % (top, crd)
                )

    return out


def run_exe(cmdline):
    """
    Simple wrapper to execute the external programs through subprocess.

    :param cmdline: complete command line as given on a shell prompt
    :type cmdline: str
    """

    logger.write('Executing command:\n%s\n' % cmdline)

    env = os.environ.copy()

    if 'COPY_LD_LIBRARY_PATH' in os.environ:
         env['LD_LIBRARY_PATH'] = os.environ['COPY_LD_LIBRARY_PATH']
    else:
         env['LD_LIBRARY_PATH'] = ''

    proc = subp.Popen(shlex.split(cmdline), stdout=subp.PIPE, stderr=subp.PIPE,
                      env=env)
    out, err =  proc.communicate()

    return proc.returncode, out, err
