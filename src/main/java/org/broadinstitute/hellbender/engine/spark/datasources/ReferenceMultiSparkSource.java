package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.FileExtensions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;

/**
 * Wrapper to load a reference sequence from a file stored on HDFS, GCS, or locally.
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
     * @param referencePathSpecifier local path or URL to the reference file
     * @param referenceWindowFunction the custom reference window function used to map reads to desired reference bases
     */
    public ReferenceMultiSparkSource( final GATKPath referencePathSpecifier,
                                      final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction) {
        Utils.nonNull(referenceWindowFunction);
        if ( ReferenceTwoBitSparkSource.isTwoBit(referencePathSpecifier)) {
            try {
                referenceSource = new ReferenceTwoBitSparkSource(referencePathSpecifier);
            } catch (IOException e) {
                throw new UserException("Failed to create a ReferenceTwoBitSource object" + e.getMessage());
            }
        } else if (referencePathSpecifier.isFasta()) {
            if (referencePathSpecifier.isHadoopURL()) {
                referenceSource = new ReferenceHadoopSparkSource(referencePathSpecifier);
            } else {
                referenceSource = new ReferenceFileSparkSource(referencePathSpecifier);
            }
        } else {
            throw new UserException.CouldNotReadInputFile("Couldn't read the given reference, reference must be a .fasta or .2bit file.\n" +
                    " Reference provided was: " + referencePathSpecifier);
        }
        this.referenceWindowFunction = referenceWindowFunction;
    }

    static boolean isFasta(final GATKPath referencePathSpecifier) {
        final String referencePathString = referencePathSpecifier.getURI().getPath();
        for (final String ext : FileExtensions.FASTA) {
            if (referencePathString.endsWith(ext)) {
                return true;
            }
        }
        return false;
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
