package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;

/**
 * Wrapper to load a reference sequence from the Google Genomics API, or a file stored on HDFS or locally.
 *
 * This class needs to subclassed by test code, so it cannot be declared final.
 */
public class ReferenceMultiSparkSource implements ReferenceSparkSource, Serializable {
    private static final long serialVersionUID = 1L;

    private ReferenceSparkSource referenceSource;
    private SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    @VisibleForTesting
    protected ReferenceMultiSparkSource() {};

    /**
     * @param referenceURL the name of the reference (if using the Google Genomics API), or a path to the reference file
     * @param referenceWindowFunction the custom reference window function used to map reads to desired reference bases
     */
    public ReferenceMultiSparkSource( final String referenceURL,
                                      final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction) {
        Utils.nonNull(referenceWindowFunction);
        if ( ReferenceTwoBitSparkSource.isTwoBit(referenceURL)) {
            try {
                referenceSource = new ReferenceTwoBitSparkSource(referenceURL);
            } catch (IOException e) {
                throw new UserException("Failed to create a ReferenceTwoBitSource object" + e.getMessage());
            }
        } else if (isFasta(referenceURL)) {
            if (BucketUtils.isHadoopUrl(referenceURL)) {
                referenceSource = new ReferenceHadoopSparkSource(referenceURL);
            } else {
                referenceSource = new ReferenceFileSparkSource(referenceURL);
            }
        } else {
            throw new UserException.CouldNotReadInputFile("Couldn't read the given reference, reference must be a .fasta or .2bit file.\n" +
                                                                  " Reference provided was: " + referenceURL);
        }
        this.referenceWindowFunction = referenceWindowFunction;
    }

    private static boolean isFasta(String reference) {
        for (final String ext : ReferenceSequenceFileFactory.FASTA_EXTENSIONS) {
            if (reference.endsWith(ext)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether this reference source can be used with Spark broadcast.
     */
    @Override
    public boolean isCompatibleWithSparkBroadcast(){
        return referenceSource.isCompatibleWithSparkBroadcast();
    }

    /**
     * @return the custom reference window function used to map reads to desired reference bases
     */
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return referenceWindowFunction;
    }

    /**
     * Return reference bases for the given interval.
     * @param interval the interval to return reference bases for
     * @return reference bases for the given interval
     */
    @Override
    public ReferenceBases getReferenceBases(final SimpleInterval interval) throws IOException {
        return referenceSource.getReferenceBases(interval);
    }

    /**
     * Return a sequence dictionary for the reference.
     * @param optReadSequenceDictionaryToMatch - (optional) the sequence dictionary of the reads, we'll match its order if possible.
     * @return sequence dictionary for the reference
     */
    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        try {
            return referenceSource.getReferenceSequenceDictionary(optReadSequenceDictionaryToMatch);
        }
        catch ( IOException e ) {
            throw new GATKException("Error getting reference sequence dictionary");
        }
    }
}
