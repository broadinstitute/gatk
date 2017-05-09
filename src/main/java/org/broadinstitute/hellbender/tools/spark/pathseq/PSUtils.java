package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.apache.logging.log4j.Logger;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.*;

/**
 * Common functions for PathSeq
 */
public final class PSUtils {

    public static JavaRDD<GATKRead> primaryReads(final JavaRDD<GATKRead> reads) {
        return reads.filter(read -> !(read.isSecondaryAlignment() || read.isSupplementaryAlignment()));
    }

    /**
     * Gets RDDs of the paired and unpaired reads
     */
    public static Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> splitPairedAndUnpairedReads(final JavaRDD<GATKRead> reads, final int readsPerPartitionGuess) {

        //Shuffle reads into partitions by read name hash code, then map each partition to a pair of lists, one
        //  containing the paired reads and the other the unpaired reads
        final JavaRDD<Tuple2<List<GATKRead>, List<GATKRead>>> repartitionedReads = reads.mapToPair(read -> new Tuple2<>(read.getName(), read))
                .partitionBy(new HashPartitioner(reads.getNumPartitions()))
                .map(Tuple2::_2)
                .mapPartitions(iter -> mapPartitionsToPairedAndUnpairedLists(iter, readsPerPartitionGuess));

        //Split the RDD into paired and unpaired reads
        repartitionedReads.cache();
        final JavaRDD<GATKRead> pairedReads = repartitionedReads.flatMap(tuple -> tuple._1.iterator());
        final JavaRDD<GATKRead> unpairedReads = repartitionedReads.flatMap(tuple -> tuple._2.iterator());
        repartitionedReads.unpersist();

        return new Tuple2<>(pairedReads, unpairedReads);
    }

    /**
     * Maps each partition to a Tuple of two Lists, the first containing the paired reads, the second containing unpaired
     */
    private static Iterator<Tuple2<List<GATKRead>, List<GATKRead>>> mapPartitionsToPairedAndUnpairedLists(final Iterator<GATKRead> iter, final int readsPerPartitionGuess) {
        //Find the paired and unpaired reads by scanning the partition for repeated names
        final List<GATKRead> pairedReadsList = new ArrayList<>(readsPerPartitionGuess);
        final Map<String, GATKRead> unpairedReads = new HashMap<>(readsPerPartitionGuess);
        while (iter.hasNext()) {
            final GATKRead read = iter.next();
            final String readName = read.getName();
            //If read's mate is already in unpairedReads then we have a pair, which gets added to the ordered List
            if (unpairedReads.containsKey(readName)) {
                pairedReadsList.add(read);
                pairedReadsList.add(unpairedReads.remove(readName));
            } else {
                unpairedReads.put(readName, read);
            }
        }
        //Get the unpaired reads out of the hashmap
        final List<GATKRead> unpairedReadsList = new ArrayList<>(unpairedReads.values());

        //Minimize unpairedReads memory footprint (don't rely on readsPerPartitionGuess)
        final List<GATKRead> pairedReadsListResized = new ArrayList<>(pairedReadsList.size());
        pairedReadsListResized.addAll(pairedReadsList);

        //Wrap the result
        final Tuple2<List<GATKRead>, List<GATKRead>> lists = new Tuple2<>(pairedReadsListResized, unpairedReadsList);
        final List<Tuple2<List<GATKRead>, List<GATKRead>>> listOfTuple = new ArrayList<>(1);
        listOfTuple.add(lists);
        return listOfTuple.iterator();
    }

    public static String[] parseCommaDelimitedArgList(final String arg) {
        if (arg == null || arg.isEmpty()) {
            return new String[0];
        }
        return arg.split(",");
    }

    /**
     * Parses command-line option for specifying kmer spacing masks
     */
    public static byte[] parseMask(final String maskArg, final int kSize) {

        final String[] kmerMaskString = parseCommaDelimitedArgList(maskArg);
        final byte[] kmerMask = new byte[kmerMaskString.length];
        for (int i = 0; i < kmerMaskString.length; i++) {
            kmerMask[i] = (byte) Integer.parseInt(kmerMaskString[i]);
            Utils.validateArg(kmerMask[i] >= 0 && kmerMask[i] < kSize, "Invalid kmer mask index: " + kmerMaskString[i]);
        }
        return kmerMask;
    }

    /**
     * Prints warning message followed by a list of relevant items
     */
    public static void logItemizedWarning(final Logger logger, final Collection<String> items, final String warning) {
        if (!items.isEmpty()) {
            String str = "";
            for (final String acc : items) str += acc + ", ";
            str = str.substring(0, str.length() - 2);
            logger.warn(warning + " : " + str);
        }
    }

    /**
     * Writes two objects using Kryo to specified local file path.
     * NOTE: using setReferences(false), which must also be set when reading the file. Does not work with nested
     * objects that reference its parent.
     */
    public static void writeKryoTwo(final String filePath, final Object obj1, final Object obj2) {
        try {
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Output output = new Output(new FileOutputStream(filePath));
            kryo.writeClassAndObject(output, obj1);
            kryo.writeClassAndObject(output, obj2);
            output.close();
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not serialize objects to file", e);
        }
    }

    /**
     * Same as GATKSparkTool's getRecommendedNumReducers(), but can specify input BAM path (for when --input is not used)
     */
    public static int pathseqGetRecommendedNumReducers(final String inputPath, final int numReducers,
                                                       final PipelineOptions options, final int targetPartitionSize) {
        if (numReducers != 0) {
            return numReducers;
        }
        return 1 + (int) (BucketUtils.dirSize(inputPath, options) / targetPartitionSize);
    }
}
