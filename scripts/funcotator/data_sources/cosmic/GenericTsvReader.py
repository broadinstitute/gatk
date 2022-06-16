#!/usr/bin/env python


"""
Created on Jul 5, 2012

@author: lichtens
"""
import csv
import os


class GenericTsvReader(object):
    """
    Read a TSV file.  
    
    This class wraps a DictReader, but handles comments, which are not handled gracefully in the python csv library. 
    
    The next() method assumes user is interested in the content, not the comments.  
        Get the comments using getComments or getCommentsAsList.  The latter assumes each comment is a line of text. 
    
    Notes:
    IMPORTANT:  At this time, this class does not support comments below the header line.  
    This class will load all comment lines into RAM at one time.  This could theoretically cause a bottleneck in some files.
    
    TODO: Low priority: Possibly make this inherit from DictReader, since it is exactly the same plus some functionality
    """

    def __init__(self, filename, commentPrepend='#', fieldNames=None, delimiter='\t'):
        """
        Constructor
        """
        self.__dict__.update(locals())
        self.inputContentFP = file(filename, 'r')
        self.commentLines = ''
        self.commentPrepend = commentPrepend
        
        # The comment lines must be loaded before the dict reader is initialized.
        self._loadCommentLines()
        
        self.dictReader = csv.DictReader(self.inputContentFP, delimiter=delimiter, fieldnames=fieldNames)

    def _loadCommentLines(self):
        resetLocation = self.inputContentFP.tell()
        nextChar = self.inputContentFP.read(1)

        # Get rid of blank lines
        while nextChar in ['\n', '\r']:
            resetLocation = self.inputContentFP.tell()
            nextChar = self.inputContentFP.read(1)
            
        while nextChar == self.commentPrepend:
            self.commentLines = self.commentLines + (self.commentPrepend + self.inputContentFP.readline())
            resetLocation = self.inputContentFP.tell()
            nextChar = self.inputContentFP.read(1)
        
        # Go back one character to make sure that we have moved the file pointer to the
        #  beginning of the first non-comment line.
        self.inputContentFP.seek(resetLocation, os.SEEK_SET)

    def next(self):
        return self.dictReader.next()
        
    def getFieldNames(self):
        return self.dictReader.fieldnames
    
    def getComments(self):
        return self.commentLines

    def getCommentsAsList(self):
        """ Return each comment line as an entry in a list """
        return self.commentLines.strip().split('\n')

    def getInputContentFP(self):
        return self.inputContentFP

    def close(self):
        self.inputContentFP.close()

    def __iter__(self):
        return self
