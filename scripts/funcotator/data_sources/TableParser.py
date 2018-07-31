#!/usr/local/bin/python
# coding: UTF-8
# modified code downloaded from:
# http://devwiki.beloblotskiy.com/index.php5/Generic_HTML_Table_parser_(python)
# mods by: Aquil H. Abdullah
from HTMLParser import HTMLParser
# import pdb

# Print debug info
markup_debug_low = not True
markup_debug_med = not True
 
class NestedTableError(Exception):
    """
    Error raised when TableParser finds a nested table.
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)

# Generic HTML table parser
class TableParser(HTMLParser):
    """
    Class to handle extracting a table from an HTML Page.
    NOTE: Does not handle Tables within 
    """
    def __init__(self):
        HTMLParser.__init__(self)
        # Can't use super HTMLParser is an old-style class
        # super(TableParser, self).__init__()
        self._tables = list() # Added to generic class
        self._curr_table = list() # Added to generic class
        self._curr_row = list()# Added to generic class
        self._curr_cell = '' # Added to generic class
        self._in_table = False # Added to generic class
        self._td_cnt = 0
        self._tr_cnt = 0
        self._curr_tag = ''
        self._colspan = 1
    def get_tables(self):
        """
        Return the list of tables scraped from html page
        """
        return self._tables

    def handle_starttag(self, tag, attrs):
        self._curr_tag = tag
        if tag.upper() == 'TABLE' and not self._in_table:
            self._in_table = True
        elif tag.upper() == 'TABLE' and self._in_table:
            raise NestedTableError("Parsing Failed Nested Table Found.")

        if tag == 'td':
            self._td_cnt += 1
            for attr in attrs:
                if attr[0].upper() == 'COLSPAN':
                    self._colspan = int(attr[1])
            self.col_start(self._td_cnt)
            if markup_debug_low: print "<TD> --- %s ---" % self._td_cnt
        elif tag == 'tr':
            self._td_cnt = 0
            self._tr_cnt += 1
            self.row_start(self._tr_cnt)
            if markup_debug_low: print "<TR> === %s ===" % self._tr_cnt
        else:
            if markup_debug_low: print "<%s>" % tag
 
    def handle_endtag(self, tag):
        if tag.upper() == 'TABLE':
            self._in_table = False
            self._tables.append(self._curr_table)
            self._curr_table = list()
        if markup_debug_low: print "</%s>" % tag
        # it's possible to check "start tag - end tag" pair here (see, tag and self._curr_tag)
        if tag == 'tr':
            self.row_finish(self._tr_cnt)
        elif tag == 'td':
            self.col_finish(self._td_cnt)
            self._colspan = 1
        else:
            pass
 
    def handle_data(self, data):
        #if markup_debug_low: print u'[%s,%s] %s: "%s"' % (self._tr_cnt, self._td_cnt, self._curr_tag, unicode(data, 'mbcs'))
        self.process_raw_data(self._tr_cnt, self._td_cnt, self._curr_tag, data)
 
    # Overridable 
    def process_raw_data(self, row, col, tag, data):
        if row > 0 and col > 0:
            self.process_cell_data(row, col, tag, data)
        else:
            pass    # outside the table
 
    # Overridable 
    def process_cell_data(self, row, col, tag, data):
        # pass
        self._curr_cell += data.strip() + ' '
 
    # Overridable 
    def row_start(self, row):
        # pass
        self._curr_row = list()

    # Overridable 
    def row_finish(self, row):
        # pass
        row = self._curr_row[:]
        self._curr_table.append(row)

    # Overridable 
    def col_start(self, col):
        # pass
        self._curr_cell = ''
 
    # Overridable 
    def col_finish(self, col):
        # pass
        self._curr_row.append(self._curr_cell)
        pad = self._colspan - 1
        if pad > 0:
            for i in range(pad):
                self._curr_row.append('')

