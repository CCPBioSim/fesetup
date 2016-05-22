#  Copyright (C) 2013-2016  Hannes H Loeffler
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
A INI like parser.
"""

__revision__ = "$Id$"


import sys, difflib



COMMENT_CHARS = '#;'
SPACE_CHARS = ' \t'

TRUE_STR = frozenset( ('1', 'yes', 'y', 'true', 't', 'on') )
FALSE_STR = frozenset( ('0', 'no', 'n', 'false', 'f', 'off') )



def strip_eol_comment(s):
    """
    Strip end-of-line comments which are a comment character preceded by a
    space or tab.
    """

    if len(s) < 3:
        return s

    for pos in range(1, len(s) ):
        if s[pos] in COMMENT_CHARS and s[pos-1] in SPACE_CHARS:
            return s[:pos-1]

    return s

    
class IniParserError(Exception):
    pass


class IniParser(dict):

    def __init__(self, *args, **kwargs):

        super(IniParser, self).__init__(*args, **kwargs)

        if not any(self):
            raise TypeError('__init__() takes at least 1 parameter: '
                            'dictionary required for object initialization')

        self.filename = None
        self.lineno = 0
        self.conv = {}

        for section in self:
            self.conv[section] = {}
            
            for key in self[section]:
                self.conv[section][key] = self[section][key][1]
                self[section][key] = self[section][key][0]

        self.section = None
        self.key = None
        self.value = None



    def format(self):
        """
        Return a formatted version of the option dictionary.
        """

        all = []

        for section in self:
            all.append('[' + section + ']')

            for key in sorted(self[section]):
                val = self[section][key]

                if val == '':
                    val = '(empty)'
                else:
                    # FIXME: special cases
                    if key == 'morph_pairs':
                        val = ', '.join('%s > %s' % (a, b) for a, b in val)
                    elif key == 'pairs':
                        val = ', '.join('%s : %s' % (a, b) for a, b in val)
                    elif type(val) in (tuple, list):
                        val = ', '.join(val)

                all.append('%s = %s' % (key, val))

            all.append('')

        return all


    def parse(self, filename, section = 'defaults'):
        """
        Poor man's INI parser.
        """

        self.filename = filename
        self.section = section
        ncont = False

        with open(self.filename, 'r') as infile:
            for line in infile:
                self.lineno += 1
                line_check = line.strip()

                if line_check == '' or line_check[0] in COMMENT_CHARS:
                    continue

                if line_check[0] == '[':
                    pos = line_check.find(']')

                    if pos < 0:
                        raise IniParserError('\ninput file error in %s: missing '
                                             '"]" in section on line %i' %
                                             (self.filename, self.lineno) )

                    self.section = line_check[1:pos]

                    continue

                if line[0] in SPACE_CHARS:
                    if not ncont:
                        raise IniParserError('\ninput file error in %s: file '
                                             'cannot start with continuation '
                                             'lines' % self.filename)

                    self.val = strip_eol_comment(line_check).strip()
                    self[self.section][self.key] += self._convert()

                    continue

                ncont = True

                try:
                    self.key, self.val = line_check.split('=', 1)
                except ValueError:
                    raise IniParserError('\ninput file error in %s: line %i '
                                         'not in "key = value" format' %
                                         (self.filename, self.lineno) )

                self.key = self.key.strip()
                self.val = strip_eol_comment(self.val).strip()

                if not self.key:
                    raise IniParserError('\ninput file error in %s: no key in '
                                         'line %i' %
                                         (self.filename, self.lineno) )

                if not self.val:
                    raise IniParserError('\ninput file error in %s: no value '
                                         'in line %i' %
                                         (self.filename, self.lineno) )

                if self.key not in self[self.section]:
                    # NOTE: fuzzy matching may need refinement
                    maybe_list = difflib.get_close_matches(self.key,
                                            self[self.section].keys(), 5, 0.5)

                    if len(maybe_list) > 1:
                        maybe_text = (', did you mean any of %s...?' %
                                      ', '.join(['"%s"' % s
                                                 for s in maybe_list] ) )
                    elif len(maybe_list) == 1:
                        maybe_text = (', did you mean "%s"?' % maybe_list[0])
                    else:
                        maybe_text = ''

                    raise IniParserError('\ninput file error in %s: unknown '
                                         'key "%s" in line %i%s'
                                         % (self.filename, self.key,
                                            self.lineno, maybe_text) )


                self[self.section][self.key] = self._convert()


    def _convert(self):
        """
        Use function in conversion table to convert value to right type.
        """

        try:
            funct = self.conv[self.section][self.key]
        except KeyError:
            print >>sys.stderr, ('missing key in conversion table: '
                                 '[%s] %s' % (self.section, self.key) )
            raise

        if not funct:
            return self.val

        if type(funct[0]) == str:
            try:
                method = getattr(self, '_str2' + funct[0])
            except AttributeError:
                print >> sys.stderr, ('unknown conversion function %s' %
                                      funct[0])
                raise

            try:
                val = method(self.val, *funct[1:])
            except TypeError:
                print >>sys.stderr, '%s is not callable ' % method
                raise
        elif callable(funct[0]):
            val = funct[0](self.val, *funct[1:])
        else:
            print >>sys.stderr, ('unkown function type for %s' % funct[0])
            raise TypeError

        if val == None:
            raise IniParserError('\ninput file error in %s: data type '
                                 'conversion of \'%s = %s\' failed in '
                                 'line %i' % (self.filename, self.key,
                                              self.val, self.lineno) )

        return val


    def _ltok(self, val, sep):
        """
        Simple tokenizer for lists with quoting for separator sep.
        """

        tok = []
        cs = []
        inquote = False

        for c in val:
            if c == '"':
                inquote = not inquote
                continue

            if c == sep and not inquote:
                tok.append(''.join(cs).strip() )
                cs = []
                continue

            cs.append(c)

        if inquote:
            return None

        tok.append(''.join(cs).strip() )

        return tok


    def _str2list(self, val, sep):
        """
        Convert a string to a list using sep as separator.  
        """

        lst = self._ltok(val, sep)

        if lst == None:
            return None

        temp = []

        for val in lst:
            if val:
                temp.append(val.strip())

        return temp


    def _str2pairlist(self, val, sep, sep2):
        """
        Convert a string to a list of pairs using sep and sep2 as separators.
        sep must be different from sep2!
        """

        lst = self._ltok(val, sep)

        if lst == None:
            return None

        temp = []

        for pair in lst:
            p = pair.strip()

            if p:
                a, b = p.split(sep2)
                temp.append( (a.strip(), b.strip()) )

        return temp


    def _str2bool(self, val):
        """
        Convert a string to a bool.  
        """

        b = val.lower()

        if b in TRUE_STR:
            return True
        elif b in FALSE_STR:
            return False

        return None



if __name__ == '__main__':

    import argparse
    
    LIST_SEP = ','
    MORPH_PAIR_SEP = '>'

    options = {}

    options['globals'] = {
        'logfile': ('dGprep.log', None),
        'forcefield': (['amber', 'ff12SB', 'tip3p'], ('list', LIST_SEP) ),
        'parmchk_version': (1, (int, ) ),
        'mcs.timeout': (600.0, (float, ) ),
        }

    options['ligand'] = {
        'basedir': ('', None),
        'molecules': ('', ('list', LIST_SEP) ),
        'morph_pairs': ('', ('pairlist', LIST_SEP, MORPH_PAIR_SEP) ),
        }

    options['protein'] = {
        'basedir': ('', None),
        'box.length': (10.0, (float, ) )
        }

    options['complex'] = {
        'pairs': ('', None),
        'neutralize': (False, ('bool', ) )
        }


    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        help = 'input file in INI format, if not given then '
                        'just output defaults')
    args = parser.parse_args()


    opts = IniParser(options)
    opts.output()
    opts.parse(args.infile, 'globals')
    opts.output()
