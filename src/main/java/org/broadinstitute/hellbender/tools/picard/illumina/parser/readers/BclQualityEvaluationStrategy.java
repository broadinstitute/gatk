package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import htsjdk.samtools.util.CollectionUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaParserException;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Describes a mechanism for revising and evaluating qualities read from a BCL file.  This class accumulates observations about low quality
 * scores that it evaluates, so distinct instances should be used for unrelated sets of BCL readers.
 * 
 * The mechanism for revising qualities is not configurable.  The qualities that are less than 1 are revised to 1, and other qualities are
 * not affected.
 *
 * This class is thread-safe and a single instance can and should be passed to {@link BclReader}s running in separate threads.
 * 
 * To replicate the functionality of {@link BclReader}s prior to the introduction of this class, create a single instance passing 
 * {@link #ILLUMINA_ALLEGED_MINIMUM_QUALITY} to the constructor, and then call {@link #assertMinimumQualities()} once the readers finish
 * their work.
 * 
 * @author mccowan
 */
public class BclQualityEvaluationStrategy {
    public static final int ILLUMINA_ALLEGED_MINIMUM_QUALITY = 2;
    private final int minimumRevisedQuality;
    /** A thread-safe defaulting map that injects an AtomicInteger starting at 0 when a uninitialized key is get-ted. */
    private Map<Byte, AtomicInteger> qualityCountMap = Collections.synchronizedMap(new CollectionUtil.DefaultingMap<Byte, AtomicInteger>(
            new CollectionUtil.DefaultingMap.Factory<AtomicInteger, Byte>() {
                @Override
                public AtomicInteger make(final Byte _) {
                    return new AtomicInteger(0);
                }
            }, true));

    /**
     * @param minimumRevisedQuality The minimum quality that should be seen from revised qualities; controls whether or not an exception
     *                              is thrown when calling {@link #assertMinimumQualities()}
     */
    public BclQualityEvaluationStrategy(final int minimumRevisedQuality) {
        this.minimumRevisedQuality = minimumRevisedQuality;
    }

    /** The rule used to revise quality scores, which is: if it's less than 1, make it 1. */
    private static byte generateRevisedQuality(final byte quality) { return (byte) Math.max(quality, 1); }
    
    /**
     * Accepts a quality read from a BCL file and (1) returns a 1 if the value was 0 and (2) makes a note of the provided quality if it is
     * low.  Because of (2) each record's quality should be passed only once to this method, otherwise it will be observed multiple times.
     *
     * @param quality The quality score read from the BCL
     * @return The revised new quality score
     */
    public byte reviseAndConditionallyLogQuality(final byte quality) {
        final byte revisedQuality = generateRevisedQuality(quality);
        if (quality < ILLUMINA_ALLEGED_MINIMUM_QUALITY) {
            qualityCountMap.get(quality).incrementAndGet();
        }
        return revisedQuality;
    }

    /**
     * Reviews the qualities observed thus far and throws an exception if any are below the minimum quality threshold.
     */
    public void assertMinimumQualities() {
        final Collection<String> errorTokens = new LinkedList<String>();
        for (final Map.Entry<Byte, AtomicInteger> entry : this.qualityCountMap.entrySet()) {
            /**
             * We're comparing revised qualities here, not observed, but the qualities that are logged in qualityCountMap are observed
             * qualities.  So as we iterate through it, convert observed qualities into their revised value. 
             */
            if (generateRevisedQuality(entry.getKey()) < minimumRevisedQuality) { 
                errorTokens.add(String.format("quality %s observed %s times", entry.getKey(), entry.getValue()));
            }
        }
        if (!errorTokens.isEmpty()) {
            throw new IlluminaReaderException(String.format(
                    "Found BCL qualities that fell beneath minimum threshold of %s: %s.",
                    minimumRevisedQuality,
                    CollectionUtil.join(errorTokens, "; ")
            ));
        }
    }

    /**
     * Returns a view of number of qualities that failed, where the key is the quality score and the value is the number of observations.
     */
    public Map<Byte, Integer> getPoorQualityFrequencies() {
        final Map<Byte, Integer> qualityCountMapCopy = new HashMap<Byte, Integer>();
        for (final Map.Entry<Byte, AtomicInteger> entry : qualityCountMap.entrySet()) {
            qualityCountMapCopy.put(entry.getKey(), entry.getValue().intValue());
        }
        return Collections.unmodifiableMap(qualityCountMapCopy);
    }
}
