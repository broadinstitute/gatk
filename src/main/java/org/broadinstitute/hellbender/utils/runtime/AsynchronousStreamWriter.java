package org.broadinstitute.hellbender.utils.runtime;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.concurrent.*;
import java.util.function.Function;

/**
 * A service that can be used to write to a stream using a thread background thread and an executor service. This
 * is typically used to write items to a buffered stream that might block until the stream is consumed by a reader.
 * @param <T> Type of items to be written.
 */
public class AsynchronousStreamWriter<T> {
    private static final Logger logger = LogManager.getLogger(AsynchronousStreamWriter.class);

    final ExecutorService executorService;
    final OutputStream streamWriter;
    final Function<T, ByteArrayOutputStream> itemSerializer;
    Future<Integer> previousBatch;

    /**
     * @param executorService executor service to be used to dispatch background tasks
     * @param streamWriter target stream to which items should be written
     * @param itemSerializer function that converts an item of type {@code T} to a {@code ByteArrayOutputStream} for serialization
     */
    public AsynchronousStreamWriter(
            final ExecutorService executorService,
            final OutputStream streamWriter,
            final Function<T, ByteArrayOutputStream> itemSerializer)
    {
        Utils.nonNull(executorService);
        Utils.nonNull(streamWriter);
        Utils.nonNull(itemSerializer);

        this.streamWriter = streamWriter;
        this.executorService = executorService;
        this.itemSerializer = itemSerializer;
        previousBatch = null;
    }

    /**
     * Request that a batch of items be written to the stream on a background thread. Any previously requested batch
     * must have already been completed and retrieved via {@link #waitForPreviousBatchCompletion}.
     *
     * @param batchList a list of items to be written
     */
    public void startBatchWrite(final List<T> batchList) {
        Utils.nonNull(batchList);
        Utils.nonEmpty(batchList);

        if (previousBatch != null) {
            throw new IllegalStateException("Previous batch not yet complete");
        }

        previousBatch = executorService.submit(() -> {
            try {
                Integer batchSize = batchList.size();
                for (int i = 0; i < batchList.size(); i++) {
                    T element = batchList.get(i);
                    itemSerializer.apply(element).writeTo(streamWriter);
                }
                // this can block, waiting for the stream to be consumed
                streamWriter.flush();
                return batchSize; // return the number of items this batch was asked to write
            } catch (IOException e) {
                throw new GATKException("IOException converting bytes for serialization", e);
            }
        });
    }

    /**
     * Waits for a batch that was previously initiated via {@link #startBatchWrite(List)}}
     * to complete, flushes the target stream and returns the corresponding completed Future. The Future representing
     * a given batch can only be obtained via this method once. If no work is outstanding, and/or the previous batch
     * has already been retrieved, null is returned.
     * @return returns null if no previous work to complete, otherwise a completed Future
     */
    public Future<Integer> waitForPreviousBatchCompletion() {
        final Future<Integer> lastCompleteBatch = previousBatch;
         if (previousBatch != null) {
            try {
                try {
                    previousBatch.get();
                } catch (ExecutionException | InterruptedException e) {
                    throw new GATKException("Interrupted during background stream write", e);
                }
                streamWriter.flush();
            } catch (IOException e) {
                throw new GATKException("IOException waiting for asynchronous batch completion", e);
            }
            previousBatch = null;
        }
        return lastCompleteBatch;
    }

    /**
     * Terminate the async writer, cancelling any outstanding work.
     * @return false if a batch was outstanding and could not be cancelled, true otherwise
     */
    public boolean terminate() {
        boolean isCancelled = true;
        if (previousBatch != null) {
            logger.warn("Cancelling outstanding asynchronous writing");
            isCancelled = previousBatch.cancel(true);
        }
        previousBatch = null;
        return isCancelled;
    }

    /**
     * Convenience function that can be provided to an {@code AsynchronousStreamWriter} to serialize String objects.
     */
    public static Function<String, ByteArrayOutputStream> stringSerializer =
            (String item) -> {
                final ByteArrayOutputStream bos = new ByteArrayOutputStream();
                try {
                    bos.write(item.getBytes());
                } catch (IOException e) {
                    throw new GATKException("IOException converting bytes for serialization", e);
                }
                return bos;
            };

}
