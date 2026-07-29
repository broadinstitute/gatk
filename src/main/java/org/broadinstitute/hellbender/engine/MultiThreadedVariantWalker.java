package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Base class for variant walkers that process each driving variant on a worker thread pool.
 *
 * <p>Subclasses implement {@link #applyParallel} (invoked on worker threads, must be thread-safe) and
 * {@link #consumeOutput} (invoked on a single drainer thread, in input variant order). Driving variants
 * are read sequentially on the main traversal thread, dispatched to a fixed-size worker pool, and their
 * results are drained in submission order so that downstream output preserves the order in which variants
 * were read.</p>
 *
 * <p><b>Thread-safety contract for subclasses:</b>
 * <ul>
 *   <li>{@link #applyParallel} runs on multiple threads concurrently. It must not mutate shared instance
 *       state without external synchronization.</li>
 *   <li>The {@link ReadsContext}, {@link ReferenceContext}, and {@link FeatureContext} objects passed to
 *       {@link #applyParallel} share underlying data sources (the engine's {@code reads}, {@code reference},
 *       and {@code features} fields). These data sources are NOT thread-safe. Subclasses that need to query
 *       them from {@link #applyParallel} must arrange external synchronization, or pre-resolve required data
 *       on the main thread.</li>
 *   <li>{@link #consumeOutput} runs on a single drainer thread; no synchronization is required there.</li>
 * </ul>
 */
public abstract class MultiThreadedVariantWalker extends VariantWalker {

    public static final String THREADS_LONG_NAME = "threads";
    public static final String THREADS_SHORT_NAME = "threads";

    /**
     * Number of worker threads to use. A value &lt; 1 selects the default
     * ({@code Runtime.getRuntime().availableProcessors() - 1}, with a floor of 1).
     */
    @Argument(fullName = THREADS_LONG_NAME,
              shortName = THREADS_SHORT_NAME,
              doc = "Number of worker threads for parallel variant processing. " +
                    "Default is (number of available processors - 1).",
              optional = true)
    @Advanced
    public int numThreads = -1;

    /**
     * Bound on the future queue, as a multiple of {@code numThreads}. Keeps memory footprint controlled
     * while still letting workers stay busy when one variant is slow.
     */
    private static final int FUTURE_QUEUE_DEPTH_PER_THREAD = 4;

    private static final long WORKER_SHUTDOWN_TIMEOUT_SECONDS = 60L;

    /**
     * Pairs a driving variant with the future result of processing it. The queue carries these in
     * submission order so the drainer can call {@link #consumeOutput} with matching input/output.
     */
    private static final class PendingResult {
        final VariantContext input;
        final Future<List<VariantContext>> future;
        PendingResult(final VariantContext input, final Future<List<VariantContext>> future) {
            this.input = input;
            this.future = future;
        }
    }

    private static final PendingResult POISON = new PendingResult(null, null);

    /**
     * Worker-thread variant processing. Returns 0 or more {@link VariantContext} objects to emit for
     * this input variant; an empty list filters the variant out.
     *
     * <p>Must be thread-safe. See class-level Javadoc for the full contract.</p>
     */
    protected abstract List<VariantContext> applyParallel(VariantContext variant,
                                                          ReadsContext readsContext,
                                                          ReferenceContext referenceContext,
                                                          FeatureContext featureContext);

    /**
     * Drainer-thread output consumption. Invoked once per processed input variant, in submission order,
     * with the original input variant and the (possibly empty) list of output variants produced for it.
     * Implementations typically buffer or write {@code outputs} to disk, and may use {@code inputVariant}
     * to maintain ordering invariants whose state depends on the input stream position.
     */
    protected abstract void consumeOutput(VariantContext inputVariant, List<VariantContext> outputs);

    /**
     * Disallowed for subclasses of this walker - override {@link #applyParallel} instead.
     */
    @Override
    public final void apply(final VariantContext variant,
                            final ReadsContext readsContext,
                            final ReferenceContext referenceContext,
                            final FeatureContext featureContext) {
        throw new GATKException(
                "MultiThreadedVariantWalker subclasses must implement applyParallel(), not apply().");
    }

    private int resolveNumThreads() {
        if (numThreads >= 1) {
            return numThreads;
        }
        if (numThreads != -1) {
            throw new UserException.BadInput(
                    "--" + THREADS_LONG_NAME + " must be >= 1 (got " + numThreads + ")");
        }
        final int available = Runtime.getRuntime().availableProcessors();
        return Math.max(1, available - 1);
    }

    @Override
    public void traverse() {
        final int n = resolveNumThreads();
        logger.info(getClass().getSimpleName() + ": dispatching variant processing across "
                    + n + " worker thread" + (n == 1 ? "" : "s"));

        final CountingReadFilter readFilter = makeReadFilter();
        final ExecutorService workers = Executors.newFixedThreadPool(n, r -> {
            final Thread t = new Thread(r, "mt-vw-worker");
            t.setDaemon(true);
            return t;
        });
        final BlockingQueue<PendingResult> futureQueue =
                new ArrayBlockingQueue<>(n * FUTURE_QUEUE_DEPTH_PER_THREAD);
        final AtomicReference<RuntimeException> drainerError = new AtomicReference<>();

        final Thread drainer = new Thread(() -> drainerLoop(futureQueue, drainerError),
                                          "mt-vw-drainer");
        drainer.setDaemon(true);
        drainer.start();

        try {
            getTransformedVariantStream(makeVariantFilter()).forEach(variant -> {
                // If the drainer has already captured a worker-side exception, stop submitting more
                // work — the queue would otherwise back-pressure forever once consumeOutput stops.
                if (drainerError.get() != null) {
                    throw new ProducerAbortException();
                }
                final SimpleInterval interval = new SimpleInterval(variant);
                final ReadsContext rc = new ReadsContext(reads, interval, readFilter);
                final ReferenceContext refC = new ReferenceContext(reference, interval);
                final FeatureContext fc = new FeatureContext(features, interval);

                // htsjdk VariantContext lazy-decodes its genotypes on first access. Force the
                // decode here, on the single producer thread, so concurrent workers never race on
                // the lazy state of an undecoded VC. Cost is paid up-front but the work would be
                // done by the first worker anyway.
                variant.getGenotypes().size();

                final Future<List<VariantContext>> f = workers.submit(
                        () -> applyParallel(variant, rc, refC, fc));
                try {
                    futureQueue.put(new PendingResult(variant, f));
                } catch (final InterruptedException ie) {
                    Thread.currentThread().interrupt();
                    throw new GATKException("Interrupted while enqueuing variant for parallel processing", ie);
                }
                progressMeter.update(interval);
            });
        } catch (final ProducerAbortException pae) {
            // Drainer captured a worker exception; fall through to cleanup and re-throw below.
        } finally {
            try {
                futureQueue.put(POISON);
            } catch (final InterruptedException ie) {
                Thread.currentThread().interrupt();
            }
            try {
                drainer.join();
            } catch (final InterruptedException ie) {
                Thread.currentThread().interrupt();
            }
            workers.shutdown();
            try {
                if (!workers.awaitTermination(WORKER_SHUTDOWN_TIMEOUT_SECONDS, TimeUnit.SECONDS)) {
                    workers.shutdownNow();
                }
            } catch (final InterruptedException ie) {
                Thread.currentThread().interrupt();
                workers.shutdownNow();
            }
        }

        if (drainerError.get() != null) {
            throw drainerError.get();
        }
    }

    private void drainerLoop(final BlockingQueue<PendingResult> futureQueue,
                             final AtomicReference<RuntimeException> error) {
        boolean draining = true;
        try {
            while (draining) {
                final PendingResult pr = futureQueue.take();
                if (pr == POISON) {
                    return;
                }
                if (error.get() != null) {
                    // Already failed. Cancel the future and discard — keep emptying the queue so
                    // the producer can unblock and reach the cleanup path. (Cancellation is a hint;
                    // the worker may still complete, that's fine.)
                    pr.future.cancel(true);
                    continue;
                }
                final List<VariantContext> out;
                try {
                    out = pr.future.get();
                } catch (final ExecutionException ee) {
                    final Throwable cause = ee.getCause() != null ? ee.getCause() : ee;
                    // Preserve UserException / GATKException so test expectations and the CLI's
                    // user-facing error reporting still see the original exception type.
                    if (cause instanceof RuntimeException) {
                        error.compareAndSet(null, (RuntimeException) cause);
                    } else {
                        error.compareAndSet(null, new GATKException(
                                "Worker thread failed during applyParallel()", cause));
                    }
                    // Fall through to drain-and-discard mode.
                    continue;
                }
                // consumeOutput is invoked on every input variant — even when outputs is empty —
                // so subclasses can maintain ordering invariants that depend on input position.
                try {
                    consumeOutput(pr.input, out == null ? Collections.emptyList() : out);
                } catch (final RuntimeException re) {
                    error.compareAndSet(null, re);
                }
            }
        } catch (final InterruptedException ie) {
            Thread.currentThread().interrupt();
        } catch (final RuntimeException re) {
            error.compareAndSet(null, re);
        }
    }

    /** Signals the producer that it should stop submitting work because the drainer has failed. */
    private static final class ProducerAbortException extends RuntimeException {
        private static final long serialVersionUID = 1L;
        ProducerAbortException() { super(null, null, false, false); }
    }
}
