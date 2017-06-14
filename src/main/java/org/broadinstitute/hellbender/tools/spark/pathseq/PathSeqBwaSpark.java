package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Aligns reads to a reference using Bwa-mem. This is a specialized version of " +
        "BwaSpark designed for the PathSeq pipeline. It takes in paired read and/or unpaired read bams (do not use --input) " +
        "and a BWA-MEM reference image file. For small input, it is recommended that the user reduce --bamPartitionSize in " +
        "order to increase parallelism.",
        oneLineSummary = "Aligns reads to reference using Bwa-mem",
        programGroup = SparkProgramGroup.class)
public final class PathSeqBwaSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Base uri for the output file(s). Paired and unpaired reads will be written to uri appended with" +
            " '.paired.bam' and '.unpaired.bam'",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String ouputUri;

    @Argument(doc = "URI to paired reads BAM",
            fullName = "pairedInput",
            optional = true)
    public String inputPaired = null;

    @Argument(doc = "URI to unpaired reads BAM",
            fullName = "unpairedInput",
            optional = true)
    public String inputUnpaired = null;

    @ArgumentCollection
    public PSBwaArgumentCollection bwaArgs = new PSBwaArgumentCollection();

    /**
     * Reads bam from path and returns tuple of the header and reads RDD
     */
    private Tuple2<SAMFileHeader, JavaRDD<GATKRead>> loadBam(final String path, final PipelineOptions options,
                                                             final ReadsSparkSource readsSource) {
        if (path == null) return null;
        if (BucketUtils.fileExists(path)) {
            final SAMFileHeader header = readsSource.getHeader(path, null);
            if (header.getSequenceDictionary() != null && !header.getSequenceDictionary().isEmpty()) {
                throw new UserException.BadInput("Input BAM should be unaligned, but found one or more sequences in the header.");
            }
            final JavaRDD<GATKRead> reads = readsSource.getParallelReads(path, null, null, bamPartitionSplitSize);
            PSBwaUtils.addReferenceSequencesToHeader(header, bwaArgs.referencePath, getReferenceWindowFunction(), options);
            return new Tuple2<>(header, reads);
        }
        logger.warn("Could not find file " + path + ". Skipping...");
        return null;
    }

    /**
     * Writes RDD of reads to path. Note writeReads() is not used because there are separate paired/unpaired outputs.
     * Header sequence dictionary is reduced to only those that were aligned to.
     */
    private void writeBam(final JavaRDD<GATKRead> reads, final String inputBamPath, final boolean isPaired,
                          final JavaSparkContext ctx, final SAMFileHeader header) {

        //Only retain header sequences that were aligned to
        reads.cache();
        final List<String> usedSequences = PSBwaUtils.getAlignedSequenceNames(reads);
        final List<SAMSequenceRecord> usedSequenceRecords = usedSequences.stream()
                .map(seqName -> header.getSequence(seqName))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        header.setSequenceDictionary(new SAMSequenceDictionary(usedSequenceRecords));

        try {
            final String outputSuffix = isPaired ? ".paired.bam" : ".unpaired.bam";
            ReadsSparkSink.writeReads(ctx, ouputUri + outputSuffix, bwaArgs.referencePath, reads, header,
                    shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                    PSUtils.pathseqGetRecommendedNumReducers(inputBamPath, numReducers, getTargetPartitionSize()));
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(ouputUri, "writing failed", e);
        }
        reads.unpersist();
    }

    /**
     * Loads a bam, aligns using the given aligner, and writes to a new bam. Returns false if the input bam could not
     * be read.
     */
    private boolean alignBam(final String inputBamPath, final PSBwaAlignerSpark aligner, final boolean isPaired,
                             final JavaSparkContext ctx, final PipelineOptions options, final ReadsSparkSource readsSource) {

        final Tuple2<SAMFileHeader, JavaRDD<GATKRead>> loadedBam = loadBam(inputBamPath, options, readsSource);
        if (loadedBam == null) return false;
        final SAMFileHeader header = loadedBam._1;
        final JavaRDD<GATKRead> reads = loadedBam._2;
        Utils.nonNull(header);
        Utils.nonNull(reads);
        if (isPaired && !header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
            throw new UserException.BadInput("Paired input BAM must be sorted by queryname");
        }

        final JavaRDD<GATKRead> alignedReads = aligner.doBwaAlignment(reads, isPaired, ctx.broadcast(header));
        writeBam(alignedReads, inputBamPath, isPaired, ctx, header);
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        if (!readArguments.getReadFiles().isEmpty()) {
            throw new UserException.BadInput("Please use --pairedInput or --unpairedInput instead of --input");
        }
        final ReadsSparkSource readsSource = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());
        final PipelineOptions options = getAuthenticatedGCSOptions();
        try (final PSBwaAlignerSpark aligner = new PSBwaAlignerSpark(ctx, bwaArgs)) {
            boolean bPairedSuccess = alignBam(inputPaired, aligner, true, ctx, options, readsSource);
            boolean bUnpairedSuccess = alignBam(inputUnpaired, aligner, false, ctx, options, readsSource);
            if (!bPairedSuccess && !bUnpairedSuccess) {
                throw new UserException.BadInput("No reads were loaded. Ensure --pairedInput and/or --unpairedInput are set and valid.");
            }
        }
    }

}
