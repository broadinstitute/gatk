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

/**
 * <p/>
 * Class GATKException.
 * <p/>
 * This exception is for errors that are beyond the user's control, such as internal pre/post condition failures
 * and "this should never happen" kinds of scenarios.
 */
public class GATKException extends RuntimeException {
    private static final long serialVersionUID = 0L;

    public GATKException( String msg ) {
        super(msg);
    }

    public GATKException( String message, Throwable throwable ) {
        super(message, throwable);
    }

    /**
     * Subtypes of GATKException for common kinds of errors
     */

    /**
     * <p/>
     * Class GATKException.CommandLineParserInternalException
     * <p/>
     * For internal errors in the command line parser not related to syntax errors in the command line itself.
     */
    public static class CommandLineParserInternalException extends GATKException {
        private static final long serialVersionUID = 0L;
        public CommandLineParserInternalException( final String s ) {
            super(s);
        }

        public CommandLineParserInternalException( final String s, final Throwable throwable ) {
            super(s, throwable);
        }
    }

    /**
     * <p/>
     * For wrapping errors that are believed to never be reachable
     */
    public static class ShouldNeverReachHereException extends GATKException {
        private static final long serialVersionUID = 0L;
        public ShouldNeverReachHereException( final String s ) {
            super(s);
        }
        public ShouldNeverReachHereException( final String s, final Throwable throwable ) {
            super(s, throwable);
        }
        public ShouldNeverReachHereException( final Throwable throwable) {this("Should never reach here.", throwable);}
    }
}

