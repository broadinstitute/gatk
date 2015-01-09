/*
* Copyright (c) 2012 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.exceptions;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

import java.io.File;

/**
 * <p/>
 * Class UserException.
 * <p/>
 * This exception is for errors that are due to user mistakes, such as non-existent or malformed files.
 */
public class UserException extends RuntimeException {
    private static final long serialVersionUID = 0L;

    public UserException(final String msg) {
        super("A USER ERROR has occurred: " + msg);
    }

    public UserException(final String message, final Throwable throwable) {
        super("A USER ERROR has occurred: " + message, throwable);
    }

    protected static String getMessage(final Throwable t) {
        final String message = t.getMessage();
        return message != null ? message : t.getClass().getName();
    }

    /**
     * Subtypes of UserException for common kinds of errors
     */

    /**
     * <p/>
     * Class UserException.CouldNotReadInputFile
     * <p/>
     * For generic errors opening/reading from input files
     */
    public static class CouldNotReadInputFile extends UserException {
        private static final long serialVersionUID = 0L;


        public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Couldn't read file because %s caused by %s", message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Couldn't read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(String file, String message) {
            super(String.format("Couldn't read file %s because %s", file, message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, getMessage(e));
        }

        public CouldNotReadInputFile(String message) {
            super(message);
        }
    }

    public static class CommandLineException extends UserException {
        public CommandLineException(String message) {
            super(String.format("Invalid command line: %s", message));
        }
    }

    public static class BadArgumentValue extends CommandLineException {
        public BadArgumentValue(String arg, String message) {
            super(String.format("Argument %s has a bad value: %s", arg, message));
        }
    }

    public static class CannotHandleGzippedRef extends UserException {
        public CannotHandleGzippedRef() {
            super("The GATK cannot process compressed (.gz) reference sequences. Please unzip the file and try again.  Sorry for the inconvenience.");
        }
    }

    /**
     * <p/>
     * Class UserException.CouldNotCreateOutputFile
     * <p/>
     * For generic errors writing to output files
     */
    public static class CouldNotCreateOutputFile extends UserException {
        private static final long serialVersionUID = 0L;

        public CouldNotCreateOutputFile(File file, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, String message) {
            super(String.format("Couldn't write file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", filename, message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Couldn't write file %s because exception %s", file.getAbsolutePath(), getMessage(e)));
        }

        public CouldNotCreateOutputFile(String message, Exception e) {
            super(message, e);
        }
    }

    public static class MalformedBAM extends UserException {
        public MalformedBAM(SAMRecord read, String message) {
            this(read.getFileSource() != null ? read.getFileSource().getReader().toString() : "(none)", message);
        }

        public MalformedBAM(File file, String message) {
            this(file.toString(), message);
        }

        public MalformedBAM(String source, String message) {
            super(String.format("SAM/BAM file %s is malformed: %s", source, message));
        }
    }

    public static class BadInput extends UserException {
        public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }


    public static class BadTmpDir extends UserException {
        public BadTmpDir(String message) {
            super(String.format("Failure working with the tmp directory %s. Override with -Djava.io.tmpdir=X on the command line to a bigger/better file system.  Exact error was %s", System.getProperties().get("java.io.tmpdir"), message));
        }
    }

    public static class MalformedGenomeLoc extends UserException {
        public MalformedGenomeLoc(String message, GenomeLoc loc) {
            super(String.format("Badly formed genome loc: %s: %s", message, loc));
        }

        public MalformedGenomeLoc(String message) {
            super(String.format("Badly formed genome loc: %s", message));
        }
    }

    /**
     * <p/>
     * Class UserException.MalformedFile
     * <p/>
     * For errors parsing files
     */
    public static class MalformedFile extends UserException {
        private static final long serialVersionUID = 0L;

        public MalformedFile(String message) {
            super(String.format("Unknown file is malformed: %s", message));
        }

        public MalformedFile(String message, Exception e) {
            super(String.format("Unknown file is malformed: %s caused by %s", message, getMessage(e)));
        }

        public MalformedFile(File f, String message) {
            super(String.format("File %s is malformed: %s", f.getAbsolutePath(), message));
        }

        public MalformedFile(File f, String message, Exception e) {
            super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, getMessage(e)));
        }

        public MalformedFile(String name, String message) {
            super(String.format("File associated with name %s is malformed: %s", name, message));
        }

        public MalformedFile(String name, String message, Exception e) {
            super(String.format("File associated with name %s is malformed: %s caused by %s", name, message, getMessage(e)));
        }
    }

    /**
     * <p/>
     * Class UserException.CommandLineParseException
     * <p/>
     * For errors parsing the user's command line
     */
    public static class CommandLineParseException extends UserException {
        private static final long serialVersionUID = 0L;
        public CommandLineParseException(final String s) {
            super("Malformed command line: " + s);
        }

        public CommandLineParseException(final String s, final Throwable throwable) {
            super("Malformed command line: " + s, throwable);
        }
    }

    public static class MissingReferenceFaiFile extends UserException {
        public MissingReferenceFaiFile( final File indexFile, final File fastaFile ) {
            super(String.format("Fasta index file %s for reference %s does not exist. Please see %s for help creating it.",
                    indexFile.getAbsolutePath(), fastaFile.getAbsolutePath(),
                    HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }

    public static class MissingReferenceDictFile extends UserException {
        public MissingReferenceDictFile( final File dictFile, final File fastaFile ) {
            super(String.format("Fasta dict file %s for reference %s does not exist. Please see %s for help creating it.",
                    dictFile.getAbsolutePath(), fastaFile.getAbsolutePath(),
                    HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }
}



