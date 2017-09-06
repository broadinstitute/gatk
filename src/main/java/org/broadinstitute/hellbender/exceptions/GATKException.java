package org.broadinstitute.hellbender.exceptions;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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

    /*
      Subtypes of GATKException for common kinds of errors
     */

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


    public static class MissingReadField extends GATKException {
        private static final long serialVersionUID = 0L;

        public MissingReadField( final String fieldName ) {
            super(String.format("Attempted to access field \"%s\" in read, but field is not present", fieldName));
        }

        public MissingReadField( final String fieldName, final String message ) {
            super(String.format("Attempted to access field \"%s\" in read, but field is not present. %s", fieldName, message));
        }

        public MissingReadField( final String fieldName, final GATKRead read ) {
            super(String.format("Attempted to access field \"%s\" in read %s, but field is not present", fieldName, read));
        }
    }

    public static class ReadAttributeTypeMismatch extends GATKException {
        private static final long serialVersionUID = 0L;

        public ReadAttributeTypeMismatch( final String attributeName, final String targetType ) {
            super(String.format("Attribute %s not of (or convertible to) type %s", attributeName, targetType));
        }

        public ReadAttributeTypeMismatch( final String attributeName, final String targetType, final Throwable throwable ) {
            super(String.format("Attribute %s not of (or convertible to) type %s", attributeName, targetType), throwable);
        }

        public ReadAttributeTypeMismatch( final SAMRecord read, final String attributeName, final String targetType) {
            super(String.format("In read %s @ %s attribute %s not of (or convertible to) type %s", read.getReadName(), "" + read.getContig() + ":" + read.getStart(), attributeName, targetType));
        }

        public ReadAttributeTypeMismatch( final SAMRecord read, final String attributeName, final String targetType, final Throwable ex) {
            super(String.format("In read %s @ %s attribute %s not of (or convertible to) type %s", read.getReadName(), "" + read.getContig() + ":" + read.getStart(), attributeName, targetType, ex));
        }
    }
}

