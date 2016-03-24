#  Copyright (C) 2012-2015  Hannes H Loeffler
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
The FESetup package.  FESetup is a tool for automatic setup of free energy
simulations like TI or MM-PBSA.

.. moduleauthor: Hannes Loeffler <Hannes.Loeffler@stfc.ac.uk>

"""


__revision__ = "$Id$"

__version__ = "0.4.0"


import sys, os

from datetime import datetime
from cStringIO import StringIO



class Logger(object):
    """
    A simple logger class which redirects output to the 'channels'
    stdout, stderr or a filename.  All output can also be suppressed.
    """


    _instance = None


    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Logger, cls).__new__(cls, *args, **kwargs)

        return cls._instance


    def __init__(self, filename = ''):
        """
        Set up channel to write to.  Add banner with current time and date.

        :param filename: name of file into which log data is written
        :type filename: string
        """
        
        if not filename:
            self.logfile = None
        elif filename == 'stdout':
            self.logfile = sys.stdout
        elif filename == 'stderr':
            self.logfile = sys.stderr
        else:
            try:
                if os.access(filename, os.F_OK):
                    self.logfile = open(filename, 'a')
                else:
                    self.logfile = open(filename, 'w')
            except IOError as why:
                sys.exit('%s' % why)

        if self.logfile:
            self.logfile.write('\n===== START ===== %s ===== START =====\n' %
                               datetime.now() )

        self.filename = filename


    def finalize(self, *args, **kwargs):
        if self.logfile:
            self.logfile.write('\n====== END ====== %s ====== END ======\n\n' %
                               datetime.now() )

            self.logfile.close()


    def write(self, message):
        """
        Write output to selected channel: stdout, stderr, filename; or suppress
        output alltogether.

        :param message: message to be written to the chosen log channel
        :type messgae: string
        """
        
        if self.logfile:
            self.logfile.write('%s\n' % message)


logger = Logger('')

def create_logger(filename):
    global logger
    logger = Logger(filename)


class DirManager(object):
    """
    Context manager for controlled entry and exit of directories.
    """

    def __init__(self, dst):
        self.topdir = os.getcwd()
        self.dst = os.path.join(self.topdir, dst)

    def __enter__(self):
        if not os.access(self.dst, os.F_OK):
            logger.write('Creating directory %s' % self.dst)
            os.makedirs(self.dst)

        logger.write('Entering directory %s' % self.dst)
        os.chdir(self.dst)

        return None

    # FIXME: check typ, value, traceback
    def __exit__(self, typ, value, traceback):
        logger.write('Entering directory %s\n' % self.topdir)
        os.chdir(self.topdir)

        return


class CaptureOutput(object):
  """
  Context manager to capture stdout and stderr.  Both outputs are
  redirected to a StringIO and made available to the caller.
  """

  def __enter__(self):
    self.oldout, self.olderr = sys.stdout, sys.stderr
    self.out = [StringIO(), StringIO()]
    sys.stdout, sys.stderr = self.out

    return self.out

  def __exit__(self, exc_type, exc_value, traceback):
    sys.stdout, sys.stderr = self.oldout, self.olderr
    self.out[0] = self.out[0].getvalue()
    self.out[1] = self.out[1].getvalue()

    return True


def report(func):
    """
    Primitive report decorator which signals start and end of a function.

    .. py:decorator:: report
    
    :param func: decorated function
    :type function: callable
    :rtype: function
    """

    def decorator(self, *args, **kwargs):
        logger.write('\n** Starting %s' % func.__name__)

        ret = func(self, *args, **kwargs)

        logger.write('** Finished with %s\n' % func.__name__)

        return ret

    # play nicely with introspection
    decorator.__name__ = func.__name__
    decorator.__doc__ = func.__doc__
    decorator.__dict__.update(func.__dict__)

    return decorator
