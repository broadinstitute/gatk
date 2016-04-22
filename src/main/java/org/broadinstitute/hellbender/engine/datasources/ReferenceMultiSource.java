package org.broadinstitute.hellbender.engine.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceTwoBitSource;
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
 * This class needs to be mocked, so it cannot be declared final.
 */
public class ReferenceMultiSource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    private ReferenceSource referenceSource;
    private SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    /**
     * @param pipelineOptions the pipeline options; must be GCSOptions if using the Google Genomics API
     * @param referenceURL the name of the reference (if using the Google Genomics API), or a path to the reference file
     * @param referenceWindowFunction the custom reference window function used to map reads to desired reference bases
     */
    public ReferenceMultiSource( final PipelineOptions pipelineOptions, final String referenceURL,
                                 final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction ) {
        Utils.nonNull(referenceWindowFunction);
        if (ReferenceTwoBitSource.isTwoBit(referenceURL)) {
            try {
                referenceSource = new ReferenceTwoBitSource(pipelineOptions, referenceURL);
            } catch (IOException e) {
                throw new UserException("Failed to create a ReferenceTwoBitSource object" + e.getMessage());
            }
        } else if (isFasta(referenceURL)) {
            if (BucketUtils.isHadoopUrl(referenceURL)) {
                referenceSource = new ReferenceHadoopSource(referenceURL);
            } else {
                referenceSource = new ReferenceFileSource(referenceURL);
            }
        } else { // use the Google Genomics API
            referenceSource = new ReferenceAPISource(pipelineOptions, referenceURL);
        }
        this.referenceWindowFunction = referenceWindowFunction;
    }

    /**
     * @param auth authentication information
     * @param referenceURL the name of the reference (if using the Google Genomics API), or a path to the reference file
     * @param referenceWindowFunction the custom reference window function used to map reads to desired reference bases
     */
    public ReferenceMultiSource(final AuthHolder auth, final String referenceURL,
                                   final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction) {
        this(auth.asPipelineOptionsDeprecated(), referenceURL, referenceWindowFunction);
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
     * @param pipelineOptions the pipeline options; must be GCSOptions if using the Google Genomics API
     * @param interval the interval to return reference bases for
     * @return reference bases for the given interval
     */
    @Override
    public ReferenceBases getReferenceBases(final PipelineOptions pipelineOptions, final SimpleInterval interval) throws IOException {
        return referenceSource.getReferenceBases(pipelineOptions, interval);
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
