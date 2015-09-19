package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.bdgenomics.adam.util.TwoBitFile;
import org.bdgenomics.utils.io.ByteAccess;
import org.bdgenomics.utils.io.ByteArrayByteAccess;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;

/**
 * Created by laserson on 9/18/15.
 */
public class ReferenceTwoBitSource implements ReferenceSource, Serializable {
    private String referenceURL;
    private TwoBitFile twoBitFile;

    public ReferenceTwoBitSource(String referenceURL) {
        if (!isTwoBit(referenceURL)) {
            throw new IllegalArgumentException("ReferenceTwoBitSource can only take .2bit files");
        }
        byte[] bytes;
        // TODO: obtain the bytes in the different cases (local file, HDFS, S3, GCS)
        if(BucketUtils.isHadoopUrl(referenceURL)) {
            bytes = ;
        }
        ByteAccess byteAccess = new DirectFullByteArrayByteAccess(bytes);
        this.twoBitFile = new TwoBitFile(byteAccess);
    }

    @Override
    public ReferenceBases getReferenceBases(PipelineOptions pipelineOptions, SimpleInterval interval) throws IOException {
        return null;
    }

    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException {
        return null;
    }

    private static boolean isTwoBit(String reference) {
        if (reference.endsWith(".2bit")) {
            return true;
        }
        return false;
    }
}
