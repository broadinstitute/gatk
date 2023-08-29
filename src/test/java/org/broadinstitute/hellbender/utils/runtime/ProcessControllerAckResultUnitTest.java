package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ProcessControllerAckResultUnitTest extends GATKBaseTest {

    @DataProvider(name="ackResultTests")
    private Object[][] getAckResultTests() {
        return new Object[][] {
                // ack result,
                //      isPositiveAck, hasMessage, getDisplayMessage prefix
                { new ProcessControllerAckResult(true),
                        true, false, StreamingToolConstants.STREAMING_ACK_MESSAGE },
                { new ProcessControllerAckResult(false),
                        false, false, StreamingToolConstants.STREAMING_NCK_MESSAGE },
                { new ProcessControllerAckResult("some message"),
                        false, true, StreamingToolConstants.STREAMING_NCK_WITH_MESSAGE_MESSAGE },
        };
    }

    @Test(dataProvider = "ackResultTests")
    private void testHasMessage(
            final ProcessControllerAckResult ackResult,
            final boolean unusedIsPositiveAck,
            final boolean hasMessage,
            final String unusedGetDisplayMessagePrefix) {
        Assert.assertEquals(ackResult.hasMessage(), hasMessage);
    }

    @Test(dataProvider = "ackResultTests")
    private void testIsPositiveAck(
            final ProcessControllerAckResult ackResult,
            final boolean isPositiveAck,
            final boolean unusedHasMessage,
            final String unusedGetDisplayMessagePrefix) {
        Assert.assertEquals(ackResult.isPositiveAck(), isPositiveAck);
    }

    @Test(dataProvider = "ackResultTests")
    private void testGetDisplayMessage(
            final ProcessControllerAckResult ackResult,
            final boolean unusedIsPositiveAck,
            final boolean unusedHasMessage,
            final String getDisplayMessagePrefix) {
        Assert.assertTrue(ackResult.getDisplayMessage().startsWith(getDisplayMessagePrefix));
    }

}
