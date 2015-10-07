package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.io.ByteStreams;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.bdgenomics.adam.models.ReferenceRegion;
import org.bdgenomics.adam.util.TwoBitFile;
import org.bdgenomics.adam.util.TwoBitRecord;
import org.bdgenomics.utils.io.ByteAccess;
import org.bdgenomics.utils.io.ByteArrayByteAccess;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.collection.JavaConversions;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * A ReferenceSource impl that is backed by a .2bit representation of a reference genome.  This loads an entire .2bit
 * file into a byte array that is encapsulated by this object.  This is particularly useful for fast reference queries
 * if the entire reference can fit into memory.
 */
public class ReferenceTwoBitSource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    public static final String TWO_BIT_EXTENSION = ".2bit";

    private String referenceURL;
    private TwoBitFile twoBitFile;

    public ReferenceTwoBitSource(PipelineOptions popts, String referenceURL) throws IOException {
        this.referenceURL = referenceURL;
        if (!isTwoBit(this.referenceURL)) {
            throw new IllegalArgumentException("ReferenceTwoBitSource can only take .2bit files");
        }
        byte[] bytes = ByteStreams.toByteArray(BucketUtils.openFile(this.referenceURL, popts));
        ByteAccess byteAccess = new DirectFullByteArrayByteAccess(bytes);
        this.twoBitFile = new TwoBitFile(byteAccess);
    }

    @Override
    public ReferenceBases getReferenceBases(PipelineOptions pipelineOptions, SimpleInterval interval) throws IOException {
        String bases = twoBitFile.extract(simpleIntervalToReferenceRegion(interval));
        return new ReferenceBases(bases.getBytes(), interval);
    }

    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException {
        Map<String, TwoBitRecord> twoBitEntries = JavaConversions.mapAsJavaMap(twoBitFile.seqRecords());
        List<SAMSequenceRecord> records = twoBitEntries.entrySet().stream()
                .map(pair -> new SAMSequenceRecord(pair.getKey(), pair.getValue().dnaSize()))
                .collect(Collectors.toList());
        return new SAMSequenceDictionary(records);
    }

    public static boolean isTwoBit(String file) {
        return file.endsWith(TWO_BIT_EXTENSION);
    }

    private static ReferenceRegion simpleIntervalToReferenceRegion(SimpleInterval interval) {
        // ADAM and GA4GH both use zero-based half-open intervals for genome coordinates
        String contig = interval.getContig();
        long start = interval.getGA4GHStart();
        long end = interval.getGA4GHEnd();
        return new ReferenceRegion(contig, start, end, null);
    }
}
