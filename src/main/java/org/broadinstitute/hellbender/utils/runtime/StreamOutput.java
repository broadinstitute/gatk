package org.broadinstitute.hellbender.utils.runtime;

/**
 * The content of stdout or stderr.
 */
public abstract class StreamOutput {
    /**
     * Empty stream output when no output is captured due to an error.
     */
    public static final StreamOutput EMPTY = new StreamOutput() {
        @Override
        public byte[] getBufferBytes() {
            return new byte[0];
        }

        @Override
        public boolean isBufferTruncated() {
            return false;
        }
    };

    /**
     * Returns the content as a string.
     *
     * @return The content as a string.
     */
    public String getBufferString() {
        return new String(getBufferBytes());
    }

    /**
     * Returns the content as a string.
     *
     * @return The content as a string.
     */
    public abstract byte[] getBufferBytes();

    /**
     * Returns true if the buffer was truncated.
     *
     * @return true if the buffer was truncated.
     */
    public abstract boolean isBufferTruncated();

    @Override
    public String toString(){
        return getBufferString();
    }
}
