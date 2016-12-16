package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongBloomFilter;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LargeQueryableLongSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.util.Collection;
import java.util.Random;

/**
 * PathSeq utilities for kmer libraries
 */
public class PSKmerLibUtils {

    public static final String HOPSCOTCH_SET_EXTENSION = ".hss";
    public static final String BLOOM_FILTER_EXTENSION = ".bfi";
    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024 * 1024 / REF_RECORD_LEN;

    /**
     * Get kmers in the reference sequence
     */
    protected static LargeLongHopscotchSet getKmersFromReference(final JavaSparkContext ctx,
                                                                 final long maxPartitionBytes,
                                                                 final int kSize,
                                                                 final int kSpace,
                                                                 final byte[] mask,
                                                                 final ReferenceMultiSource ref,
                                                                 final PipelineOptions options,
                                                                 final SAMSequenceDictionary refDict) {
        // Generate reference sequence RDD
        final JavaRDD<byte[]> refRDD = SVUtils.getRefRDD(ctx, kSize, ref, options, refDict, REF_RECORD_LEN, REF_RECORDS_PER_PARTITION);

        // Kmerize, mapping map to the collection of all canonicalized kmers and returning result to the driver. Note we do
        // not remove duplicates at this stage because the shuffle is expensive. Assuming there are relatively few
        // duplicates (~10-20% for hg38), we will collect all kmers and rely on HopscotchSet to do the work later.
        final Collection<long[]> kmerLists = refRDD.map(seq ->
                SVKmerizer.stream(seq, kSize, kSpace, new SVKmerShort(kSize))
                        .map(kmer -> kmer.mask(mask, kSize).canonical(kSize - mask.length))
                        .mapToLong(SVKmer::getLong)
                        .toArray())
                .collect();

        return longArrayCollectionToKmerSet(kmerLists, maxPartitionBytes);
    }

    /**
     * Converts a Collection of long[]'s into a Hopscotch set
     */
    protected static LargeLongHopscotchSet longArrayCollectionToKmerSet(final Collection<long[]> kmerLists, final long maxPartitionBytes) {

        //Count the number of kmers so we can initialize the HopscotchSet with an appropriate capacity
        long numKmers = 0;
        for (final long[] arr : kmerLists) {
            numKmers += arr.length;
        }

        //Construct the Hopscotch set
        final LargeLongHopscotchSet kmerSet = new LargeLongHopscotchSet(maxPartitionBytes, numKmers);
        for (final long[] arr : kmerLists) {
            kmerSet.addAll(arr);
        }

        return kmerSet;
    }

    /**
     * Applies a mask to a set of kmers. If the kmers were loaded from a file and are already masked, their lengths
     * will be the original kmer size minus the number of masked bases.
     * For example, if the kmer size was 31 and bases 0 and 15 were masked, then kmerSize should be 29
     * when loaded in here. The argument kmerMask is then applied on top of the current masked kmer. The same rule
     * applies to any subsequent iterations of writing/loading masked kmers.
     */
    protected static LargeLongHopscotchSet maskKmers(LargeLongHopscotchSet kmerSet,
                                                     final int kmerSize,
                                                     final byte[] kmerMask,
                                                     final long maxPartitionBytes) {
        if (kmerMask.length > 0) {
            final LargeLongHopscotchSet fullKmerSet = kmerSet;
            kmerSet = new LargeLongHopscotchSet(maxPartitionBytes, fullKmerSet.size());
            final LongIterator itr = fullKmerSet.iterator();
            while (itr.hasNext()) {
                final SVKmer kmer = new SVKmerShort(itr.next()).mask(kmerMask, kmerSize).canonical(kmerSize - kmerMask.length);
                kmerSet.add(kmer.getLong());
            }
        }
        return kmerSet;
    }

    /**
     * Returns new Hopscotch derived from an input set, where each item is retained with a given probability.
     * Note the number of actual reads in the downsampled set is binomially distributed.
     */
    public static LargeLongHopscotchSet downsampleKmerSet(LargeLongHopscotchSet kmerSet, final long maxPartitionBytes, final double downsampleProbability, final long downsampleSeed) {

        if (downsampleProbability < 1) {
            final LargeLongHopscotchSet fullKmerSet = kmerSet;
            final long newCapacity = (long) Math.ceil(fullKmerSet.size() * downsampleProbability);
            kmerSet = new LargeLongHopscotchSet(maxPartitionBytes, newCapacity);

            final Random rand = new Random(downsampleSeed);
            final LongIterator itr = fullKmerSet.iterator();
            while (itr.hasNext()) {
                final long val = itr.next();
                if (rand.nextDouble() < downsampleProbability) {
                    kmerSet.add(val);
                }
            }
        }
        return kmerSet;
    }

    /**
     * Converts a Hopscotch set into a Bloom filter
     */
    public static LargeLongBloomFilter createBloomFilterFromSet(final LargeLongHopscotchSet kmerSet, final long maxPartitionBytes, final double bloomFilterFPP) {

        final LargeLongBloomFilter bloomFilter = new LargeLongBloomFilter(maxPartitionBytes, kmerSet.size(), bloomFilterFPP);
        final LongIterator itr = kmerSet.iterator();
        while (itr.hasNext()) {
            bloomFilter.add(itr.next());
        }
        return bloomFilter;
    }

    /**
     * Writes an object to a URI using Kryo serialization. Appends mandatoryExtension to the URI if specified
     */
    public static void writeKryoObject(final Object obj, String uri, final PipelineOptions options, final String mandatoryExtension) {

        if (mandatoryExtension != null && !uri.toLowerCase().endsWith(mandatoryExtension.toLowerCase())) {
            uri = uri + mandatoryExtension;
        }
        final Output output = new Output(BucketUtils.createFile(uri, options));
        final Kryo kryo = new Kryo();
        kryo.setReferences(false);
        kryo.writeClassAndObject(output, obj);
        output.close();
    }

    public static void writeLargeLongHopscotchSet(final LargeLongHopscotchSet set, final String uri, final PipelineOptions options) {
        writeKryoObject(set, uri, options, HOPSCOTCH_SET_EXTENSION);
    }

    public static void writeLargeLongBloomFilter(final LargeLongBloomFilter set, final String uri, final PipelineOptions options) {
        writeKryoObject(set, uri, options, BLOOM_FILTER_EXTENSION);
    }

    public static LargeQueryableLongSet readLargeQueryableSet(final String uri, final PipelineOptions options) {

        final Input input = new Input(BucketUtils.openFile(uri, options));
        final Kryo kryo = new Kryo();

        //There is an odd bug where a ClassNotFound exception gets thrown in Spark mode when Kryo-deserializing the library
        //This has something to do with the default system ClassLoader not loading GATK classes when in Spark mode
        //Working solution: explicitly replace the Kryo system class loader with GATKTool's
        final ClassLoader loader = GATKTool.class.getClassLoader();
        kryo.setClassLoader(loader);

        kryo.setReferences(false);
        return (LargeQueryableLongSet) kryo.readClassAndObject(input);
    }

    public static LargeLongHopscotchSet readLargeLongHopscotchSet(final String uri, final PipelineOptions options) {
        final Input input = new Input(BucketUtils.openFile(uri, options));
        final Kryo kryo = new Kryo();

        //See PSKmerLibUtils.readLargeQueryableSet
        final ClassLoader loader = GATKTool.class.getClassLoader();
        kryo.setClassLoader(loader);

        kryo.setReferences(false);
        return (LargeLongHopscotchSet) kryo.readClassAndObject(input);
    }
}
