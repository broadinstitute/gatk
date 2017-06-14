package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongBloomFilter;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

/**
 * PathSeq utilities for kmer libraries
 */
public class PSKmerUtils {

    public static final String HOPSCOTCH_SET_EXTENSION = ".hss";
    public static final String BLOOM_FILTER_EXTENSION = ".bfi";
    private static final Logger logger = LogManager.getLogger(PSKmerUtils.class);

    /**
     * Gets all kmers from a given reference as a Collection of long arrays, while logging handy progress messages.
     */
    protected static Collection<long[]> getMaskedKmersFromLocalReference(final ReferenceFileSource ref, final int kSize,
                                                                         final int kSpace, final SVKmerShort mask) {

        //Load reference records
        final Map<String, ReferenceBases> records;
        try {
            records = ref.getAllReferenceBases();
        } catch (final IOException e) {
            throw new GATKException("Could not get reference bases");
        }

        final Collection<long[]> collection = new ArrayList<>(records.size());

        //Initialize progress state
        final long totalBases = records.values().stream().mapToLong(refBases -> refBases.getBases().length).sum();
        logger.info("Generating kmers from " + totalBases + " bases in " + records.size() + " records...");
        final ProgressCounter counter = new ProgressCounter(totalBases, 1e6, "Mbp", logger);

        //Get kmers from each reference record
        for (final String recName : records.keySet()) {

            //Kmerize the record
            final byte[] bases = records.get(recName).getBases();
            final long[] list = SVKmerizer.stream(bases, kSize, kSpace, new SVKmerShort(kSize))
                    .mapToLong(kmer -> new Long(PSKmerCollection.canonicalizeAndMask((SVKmerShort) kmer, kSize, mask)))
                    .toArray();

            //Add kmers to the result
            collection.add(list);

            counter.update(bases.length);
        }
        logger.info("Finished generating kmers!");
        return collection;
    }

    /**
     * Returns the total number of longs in the Collection
     */
    public static long longArrayCollectionSize(final Collection<long[]> lists) {
        return lists.stream().mapToLong(arr -> arr.length).sum();
    }

    /**
     * Converts a Collection of Lists of Longs's into a Hopscotch set
     */
    protected static LargeLongHopscotchSet longArrayCollectionToSet(final Collection<long[]> longs, final long numLongs) {
        final LargeLongHopscotchSet kmerHopscotchSet = new LargeLongHopscotchSet(numLongs);
        for (final long[] array : longs) {
            for (final long val : array) {
                kmerHopscotchSet.add(val);
            }
        }
        return kmerHopscotchSet;
    }

    /**
     * Converts a Collection of Lists of Longs's into a Bloom filter
     */
    protected static LongBloomFilter longArrayCollectionToBloomFilter(final Collection<long[]> longs, final long numLongs, final double bloomFpp) {
        final LongBloomFilter bloomFilter = new LongBloomFilter(numLongs, bloomFpp);
        final ProgressCounter counter = new ProgressCounter(numLongs, 1e6, "million kmers", logger);
        for (final long[] array : longs) {
            bloomFilter.addAll(array);
            counter.update(array.length);
        }
        return bloomFilter;
    }

    /**
     * Writes an object to a URI using Kryo serialization.
     */
    public static void writeKryoObject(final Object obj, String uri) {
        final Output output = new Output(BucketUtils.createFile(uri));
        final Kryo kryo = new Kryo();
        kryo.writeObject(output, obj);
        output.close();
    }

    public static void writeKmerSet(final String uri, final PSKmerSet set) {
        String filePath = uri;
        if (HOPSCOTCH_SET_EXTENSION != null && !uri.toLowerCase().endsWith(HOPSCOTCH_SET_EXTENSION.toLowerCase())) {
            filePath = filePath + HOPSCOTCH_SET_EXTENSION;
        }
        writeKryoObject(set, filePath);
    }

    public static void writeKmerBloomFilter(final String uri, final PSKmerBloomFilter bloomFilter) {
        String filePath = uri;
        if (BLOOM_FILTER_EXTENSION != null && !uri.toLowerCase().endsWith(BLOOM_FILTER_EXTENSION.toLowerCase())) {
            filePath = filePath + BLOOM_FILTER_EXTENSION;
        }
        writeKryoObject(bloomFilter, filePath);
    }

    public static PSKmerCollection readKmerFilter(final String uri) {
        final Input input = new Input(BucketUtils.openFile(uri));
        final Kryo kryo = new Kryo();
        if (uri.endsWith(HOPSCOTCH_SET_EXTENSION)) {
            return kryo.readObject(input, PSKmerSet.class);
        } else if (uri.endsWith(BLOOM_FILTER_EXTENSION)) {
            return kryo.readObject(input, PSKmerBloomFilter.class);
        }
        throw new UserException.BadInput("Unknown kmer set extension in file name " + uri);
    }

    private final static class ProgressCounter {
        long processedItems, processedItemsSinceLast;
        final long initialTime, totalItems, itemsInterval;
        final double normFactor;
        final Logger logger;
        final String name;

        public ProgressCounter (final long totalItems, final double normFactor, final String name, final Logger logger) {
            this.processedItems = 0;
            this.processedItemsSinceLast = 0;
            this.totalItems = totalItems;
            this.itemsInterval = totalItems / 100;
            this.initialTime = System.currentTimeMillis();
            this.normFactor = normFactor;
            this.logger = logger;
            this.name = name;
        }

        public void update(final long numItems) {
            processedItems += numItems;
            processedItemsSinceLast += numItems;
            if (processedItemsSinceLast > itemsInterval) {
                final long currentTime = System.currentTimeMillis();
                final double percentComplete = 100.0 * processedItems / (double) totalItems;
                final double averageRate = processedItems / (normFactor * (currentTime - initialTime) / (1000.0 * 60.0));
                final double estimatedTime = ((totalItems - processedItems) / normFactor) / averageRate;
                final String percentStr = String.format("%1$.1f", percentComplete);
                final String totalStr = String.format("%1$.1f", processedItems / normFactor);
                final String rateStr = String.format("%1$.1f", averageRate);
                final String timeStr = String.format("%1$.2f", estimatedTime);
                logger.info(percentStr + "% complete - " + totalStr + " " + name + " at " + rateStr + " " + name + "/min, " + timeStr + " min remaining");
                processedItemsSinceLast = 0;
            }
        }
    }
}
