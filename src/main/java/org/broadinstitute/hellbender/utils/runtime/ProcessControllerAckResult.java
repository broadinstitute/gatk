package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Command acknowledgements that are returned from a process managed by StreamingProcessController.
 * Ack results can be positive, negative, or negative with a message. Positive acks never have a messsage,
 * negative acks may optionally have a message.
 */
public class ProcessControllerAckResult {

    private final boolean isPositiveAck;
    private final String message;

    // three message types can be used by the remote process
    private static String displayMessageFormat = "%s received\n\n";
    private static String ACK_LOG_MESSAGE               = String.format(displayMessageFormat, StreamingToolConstants.STREAMING_ACK_MESSAGE);
    private static String NCK_LOG_MESSAGE               = String.format(displayMessageFormat, StreamingToolConstants.STREAMING_NCK_MESSAGE);
    private static String NCK_WITH_MESSAGE_LOG_MESSAGE  = String.format(displayMessageFormat, StreamingToolConstants.STREAMING_NCK_WITH_MESSAGE_MESSAGE);

    /**
     * Creates an ack result, for ACK or NCK.
     * @param isPositiveAck true for a positive ack (ACK), false for negative ack (NCK)
     */
    public ProcessControllerAckResult(final boolean isPositiveAck) {
        this.isPositiveAck = isPositiveAck;
        this.message = null;
    }

    /**
     * Creates an (negative NKM) ack result, with a message.
     * @param nckMessage Message detail indicating reason for negative ack (NKM).
     */
    public ProcessControllerAckResult(final String nckMessage) {
        this.isPositiveAck = false;
        this.message = nckMessage;
    }

    /**
     * @return true if this represents a positive ack, otherwise false
     */
    public boolean isPositiveAck() {
        return isPositiveAck;
    }

    /**
     * @return true if this ack is negative and includes a message
     */
    public boolean hasMessage() {
        return !isPositiveAck() && message != null && !message.isEmpty();
    }

    /**
     * @return A message string representing this ack/nck suitable for logging/display to the user.
     */
    public String getDisplayMessage() {
        if (isPositiveAck()) {
            return ACK_LOG_MESSAGE;
        } else if (hasMessage()) {
            return String.format("%s: %s", NCK_WITH_MESSAGE_LOG_MESSAGE, message);
        } else {
            return NCK_LOG_MESSAGE;
        }
    }
}
