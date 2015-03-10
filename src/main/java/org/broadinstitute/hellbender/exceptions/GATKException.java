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

