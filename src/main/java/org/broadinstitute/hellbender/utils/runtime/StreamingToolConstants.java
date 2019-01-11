package org.broadinstitute.hellbender.utils.runtime;

/**
 * Various constants used by StreamingProcessController that require synchronized equivalents in
 * the companion process, i.e., if the streaming process is written in Python, there must be
 * equivalent Python constants for use by the Python code.
 *
 * See the equivalents for Python in toolcontants.py.
 */
public class StreamingToolConstants {
    /**
     * Command acknowledgement messages used to signal positive acknowledgement ('ack'),
     * negative acknowledgement ('nck'), and negative acknowledgement with an accompanying
     * message ('nkm').
     */
    public static String STREAMING_ACK_MESSAGE = "ack";
    public static String STREAMING_NCK_MESSAGE = "nck";
    public static String STREAMING_NCK_WITH_MESSAGE_MESSAGE = "nkm";

    // This is only used by Java, but is kept here since it represents the length of the constant
    // strings defined above.
    protected static int STREAMING_ACK_MESSAGE_SIZE = 3; // "ack", "nck", or "nkm"

    /**
     * Number of characters used to represent the length of the serialized message, fixed at a constant
     * 4 characters to ensure we can deterministically know how much input to wait for when looking for
     * a message length in the incoming stream.
     */
    public static int STREAMING_NCK_WITH_MESSAGE_MESSAGE_LEN_SIZE = 4;
    public static int STREAMING_NCK_WITH_MESSAGE_MAX_MESSAGE_LENGTH = 9999;
}
