#  Copyright (C) 2014  Hannes H Loeffler
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
MD engine loader
"""

__revision__ = "$Id$"



import sys



def get_MDEngine(mdengine):
    """
    MD engine factory.
    """

    try:
        mde = __import__('FESetup.prepare.mdengines.' + mdengine,
                         fromlist = '*')
    except ImportError as detail:
        print >> sys.stderr, 'Error: %s' % detail
        print >> sys.stderr, 'Ensure all modules are properly installed'
        sys.exit(1)
    except AttributeError as detail:
        print >> sys.stderr, 'Error: %s' % detail
        print >> sys.stderr, 'Failed to properly initialize %s' % ff
        sys.exit(1)

    return mde.MDEngine
