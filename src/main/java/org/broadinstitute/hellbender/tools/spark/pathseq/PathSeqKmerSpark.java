package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.*;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.util.List;

/**
 * SparkTool to build kmer library from a given host reference. This library is required for PathSeqFilterSpark.
 */
@CommandLineProgramProperties(summary="Builds library of reference kmers used in the PathSeq subtraction phase.",
        oneLineSummary="Builds library of reference kmers",
        programGroup = SparkProgramGroup.class)
public final class PathSeqKmerSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;
    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    @Argument(doc = "file for kmer output",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String OUTPUT_FILE;

    @Argument(doc = "Kmer size",
            fullName = "kSize",
            shortName = "kSize",
            optional = true)
    private int KMER_SIZE = 31;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /** Get the list of distinct kmers in the reference, and write them to a file as a HopScotchSet. */
    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        final PipelineOptions options = getAuthenticatedGCSOptions();
        final ReferenceMultiSource referenceMultiSource = getReference();

        final List<SVKmer> kmerList = findKmers(ctx, KMER_SIZE, referenceMultiSource, options, dict);
        final HopscotchSet<SVKmer> kmerSet = new HopscotchSet<>(kmerList);

        final Output output = new Output(BucketUtils.createFile(OUTPUT_FILE));
        final Kryo kryo=new Kryo();
        kryo.setReferences(false);
        kryo.writeClassAndObject(output, kmerSet);
        output.close();

    }

    /** Get kmers in the reference sequence */
    private static List<SVKmer> findKmers(final JavaSparkContext ctx,
                                               final int kSize,
                                               final ReferenceMultiSource ref,
                                               final PipelineOptions options,
                                               final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final JavaRDD<byte[]> refRDD = SVUtils.getRefRDD(ctx, kSize, ref, options, readsDict, REF_RECORD_LEN, REF_RECORDS_PER_PARTITION);

        // List of kmers
        return processRefRDD(kSize, refRDD);
    }

    /**
     * Turn a text file of overlapping records from a reference sequence into an RDD, and do a classic map/reduce:
     * Kmerize, mapping map to the collection of all canonicalized kmers and returning result to the driver. Note we do
     * not remove duplicates at this stage because the shuffle is expensive. Assuming there are relatively few
     * duplicates, we will collect all kmers and rely on HopscotchSet to do the work later.
     */
    private static List<SVKmer> processRefRDD(final int kSize, final JavaRDD<byte[]> refRDD ) {
        return refRDD.flatMap(seq ->
                    SVKmerizer.stream(seq, kSize, new SVKmerShort(kSize))
                        .map(kmer -> kmer.canonical(kSize))
                        .collect(SVUtils.arrayListCollector(Math.max(0,seq.length-kSize+1)))
                        .iterator())
                .collect();
    }
}
