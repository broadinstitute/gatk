package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Command acknowledgements that are returned from a process managed by StreamingProcessController.
 * Ack results can be positive, negative, or negative with a message.
 */
public class ProcessControllerAckResult {

    private final boolean isPositiveAck;
    private final String message;

    // three message types can be used by the remote process
    private static String ACK_LOG_MESSAGE               = "Ack received\n\n";
    private static String NCK_LOG_MESSAGE               = "Nck received\n\n";
    private static String NCK_WITH_MESSAGE_LOG_MESSAGE  = "Nkm received\n\n";

    public ProcessControllerAckResult(final boolean isPositiveAck) {
        this(isPositiveAck, "");
    }

    public ProcessControllerAckResult(final boolean isPositiveAck, final String message) {
        this.isPositiveAck = isPositiveAck;
        this.message = message;
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
        return !isPositiveAck() && !message.isEmpty();
    }

    /**
     *
     * @return A (possibly empty) String with any message sent from the remote process.
     * Only defined for negative acknowledgements {@link #hasMessage()}.
     */
    public String getNegativeACKMessage() {
        if (isPositiveAck()) {
            throw new GATKException("Can only retrieve messages for negative acknowledgements");
        }
        return message;
    }

    /**
     * @return A message string representing this ack/nck suitable for logging/display to the user.
     */
    public String getDisplayMessage() {
        if (isPositiveAck()) {
            return ACK_LOG_MESSAGE;
        } else if (hasMessage()) {
            return String.format("%s: %s", NCK_WITH_MESSAGE_LOG_MESSAGE, getNegativeACKMessage());
        } else {
            return NCK_LOG_MESSAGE;
        }
    }
}
