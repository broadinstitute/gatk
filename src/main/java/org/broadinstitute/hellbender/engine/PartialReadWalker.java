package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Spliterator;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import java.util.stream.Stream;

/**
 * A specialized read walker that may be gracefully stopped before the input stream ends
 * 
 * A tool derived from this class should implement {@link PartialReadWalker#shouldExitEarly(GATKRead)}
 * to indicate when to stop. This method is called before {@link ReadWalker#apply(GATKRead, ReferenceContext, FeatureContext)}
 *
 */
abstract public class PartialReadWalker extends ReadWalker {

    /**
     * traverse is overridden to consult the implementation class whether to stop
     *
     * The stoppage is implemented using a custom forEach method to compensate for the
     * lack of .takeWhile() in Java 8
     */

    @Override
    public void traverse() {

        final CountingReadFilter countedFilter = makeReadFilter();
        breakableForEach(getTransformedReadStream(countedFilter), (read, breaker) -> {

            // check if we should stop
            if ( shouldExitEarly(read) ) {
                breaker.set(true);
            } else {
                // this is the body of the iteration
                final SimpleInterval readInterval = getReadInterval(read);
                apply(read,
                        new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                        new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null

                progressMeter.update(readInterval);
            }
        });

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Method to be overridden by the implementation class to determine when to stop the read stream traversal
     * @param read - the read to be processed next (in case it is needed)
     * @return boolean indicator: true means stop!
     */
    protected abstract boolean shouldExitEarly(GATKRead read);

    /**
     * Java 8 does not have a .takeWhile() on streams. The code below implements a custom forEach to allow
     * breaking out of a stream prematurely.
     *
     * code adapted from: https://www.baeldung.com/java-break-stream-foreach
     */
    private static <T> void breakableForEach(Stream<T> stream, BiConsumer<T, AtomicBoolean> consumer) {
        Spliterator<T> spliterator = stream.spliterator();
        boolean hadNext = true;
        AtomicBoolean breaker = new AtomicBoolean();

        while (hadNext && !breaker.get()) {
            hadNext = spliterator.tryAdvance(elem -> {
                consumer.accept(elem, breaker);
            });
        }
    }
}
