package org.broadinstitute.hellbender.utils.runtime;

import htsjdk.samtools.util.BufferedLineReader;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

import org.testng.Assert;
import org.testng.annotations.Test;

public class AsynchronousStreamWriterUnitTest {
    // Use a relatively long timeout when running tests to allow for variation in test environment run times
    final static TimeUnit TIMEOUT_TIMEUNIT = TimeUnit.SECONDS;
    final static int TIMEOUT_TIME = 10;

    @Test
    public void testAsyncWriteInBatches() throws IOException, InterruptedException, ExecutionException {
        final int ITEM_COUNT = 100;
        final int BATCH_SIZE = 12;

        final List<String> readCommandStrings = new ArrayList<>();
        final List<String> expectedReadCommandStrings = new ArrayList<>();
        final List<String> expectedItems = new ArrayList<>();
        List<String> batchItems;
        int batchCount = 0;
        AsynchronousStreamWriter<String> asyncWriteService = null;

        // create an executor service, and write ITEM_COUNT lines in batches of BATCH_SIZE (where ITEM_COUNT is not
        // an integral multiple of BATCH_SIZE) and leave one final odd-sized batch at the end
        final ExecutorService executorService = Executors.newSingleThreadExecutor();
        try (final ByteArrayOutputStream streamWriter = new ByteArrayOutputStream()) {
            asyncWriteService = new AsynchronousStreamWriter<>(executorService, streamWriter, AsynchronousStreamWriter.stringSerializer);
            batchItems = new ArrayList<>(BATCH_SIZE);
            for (int i = 0; i < ITEM_COUNT; i++) {
                if (batchCount == BATCH_SIZE) {
                    dispatchABatch(asyncWriteService, batchItems, batchCount, readCommandStrings, expectedReadCommandStrings);
                    batchItems = new ArrayList<>(BATCH_SIZE);
                    batchCount = 0;
                }
                final String itemToWrite = Integer.toString(i) + "\n";
                batchItems.add(itemToWrite);
                expectedItems.add(itemToWrite);
                batchCount++;
            }

            // write the last, odd-sized batch
            if (batchCount != 0) {
                dispatchABatch(asyncWriteService, batchItems, batchCount, readCommandStrings, expectedReadCommandStrings);
            }

            final Future<Integer> batchResult = asyncWriteService.waitForPreviousBatchCompletion();
            Assert.assertTrue(batchResult.get().equals(batchCount));
            Assert.assertTrue(asyncWriteService.terminate());

            // read the output stream file in and make sure everything was written to the stream by the async writer,
            // and all command strings were executed on the background thread
            try (final ByteArrayInputStream is= new ByteArrayInputStream(streamWriter.toByteArray());
                 final BufferedLineReader br = new BufferedLineReader(is)) {
                expectedItems.forEach(expectedLine -> Assert.assertEquals(br.readLine() + '\n', expectedLine));
            }
            Assert.assertEquals(readCommandStrings, expectedReadCommandStrings);
        } finally {
            if (asyncWriteService != null) {
                asyncWriteService.terminate();
            }
            executorService.shutdown();
        }
    }

    private void dispatchABatch(
            final AsynchronousStreamWriter<String> asyncWriteService,
            final List<String> batchItems,
            int batchSize,
            final List<String> readCommandStrings,
            final List<String> expectedReadCommandStrings)
    {
        final String commandString = "This would usually be a command to tell python to read %d items\n";

        // wait for the last batch to complete before we start a new one
        final Future<Integer> batchResult = asyncWriteService.waitForPreviousBatchCompletion();
        readCommandStrings.add(String.format(commandString, batchSize));
        asyncWriteService.startBatchWrite(batchItems);

        // write to our expected read string list
        expectedReadCommandStrings.add(String.format(commandString, batchSize));
    }

    @Test
    public void testDuplicateWaitForPreviousBatch() throws IOException {
        final int BATCH_SIZE = 5;
        final List<String> batchItems = new ArrayList<>(BATCH_SIZE);
        for (int i = 0; i < BATCH_SIZE; i++) { batchItems.add(Integer.toString(i)); }

        final ExecutorService executorService = Executors.newSingleThreadExecutor();
        AsynchronousStreamWriter<String> asyncWriteService = null;
        try (final ByteArrayOutputStream streamWriter = new ByteArrayOutputStream()) {
            asyncWriteService = new AsynchronousStreamWriter<>(executorService, streamWriter, AsynchronousStreamWriter.stringSerializer);
            asyncWriteService.startBatchWrite(batchItems);
            Future<Integer> batchResult = asyncWriteService.waitForPreviousBatchCompletion();
            Assert.assertNotNull(batchResult);
            batchResult = asyncWriteService.waitForPreviousBatchCompletion();
            Assert.assertNull(batchResult);
        } finally {
            if (asyncWriteService != null) {
                asyncWriteService.terminate();
            }
            executorService.shutdown();
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testWriteBeforePreviousBatchComplete() throws IOException {
       final int BATCH_SIZE = 5;
        final List<String> batchItems = new ArrayList<>(BATCH_SIZE);
        for (int i = 0; i < BATCH_SIZE; i++) { batchItems.add(Integer.toString(i)); }

        AsynchronousStreamWriter<String> asyncWriteService = null;
        final ExecutorService executorService = Executors.newSingleThreadExecutor();
        try (final ByteArrayOutputStream streamWriter = new ByteArrayOutputStream()) {
            asyncWriteService = new AsynchronousStreamWriter<>(executorService, streamWriter, AsynchronousStreamWriter.stringSerializer);
            asyncWriteService.startBatchWrite(batchItems);
            // try to start a new batch without retrieving the completion of the previous one
            asyncWriteService.startBatchWrite(batchItems);
        } finally {
            if (asyncWriteService != null) {
                asyncWriteService.terminate();
            }
            executorService.shutdown();
        }
    }

}
