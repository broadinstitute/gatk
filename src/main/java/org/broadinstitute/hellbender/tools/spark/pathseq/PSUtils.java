package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.Collection;

/**
 * Common functions for PathSeq
 */
public final class PSUtils {

    public static JavaRDD<GATKRead> primaryReads(final JavaRDD<GATKRead> reads) {
        return reads.filter(read -> !(read.isSecondaryAlignment() || read.isSupplementaryAlignment()));
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
