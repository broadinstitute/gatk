package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.MetagenomicsProgramGroup;
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

@DocumentedFeature
@CommandLineProgramProperties(summary = "Aligns reads to a reference using Bwa-mem. This is a specialized version of " +
        "BwaSpark designed for the PathSeq pipeline. User must supply unaligned paired read and/or unpaired reads " +
        "(do not use --input) and a BWA reference image file created using BwaMemIndexImageCreator. For small " +
        "input, it is recommended that the user reduce --bam-partition-size in order to increase parallelism.",
        oneLineSummary = "Step 2: Aligns reads to the pathogen reference",
        programGroup = MetagenomicsProgramGroup.class)
@BetaFeature
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
    private Tuple2<SAMFileHeader, JavaRDD<GATKRead>> loadBam(final String path,
                                                             final ReadsSparkSource readsSource) {
        if (path == null) return null;
        if (BucketUtils.fileExists(path)) {
            final SAMFileHeader header = readsSource.getHeader(path, null);
            if (header.getSequenceDictionary() != null && !header.getSequenceDictionary().isEmpty()) {
                throw new UserException.BadInput("Input BAM should be unaligned, but found one or more sequences in the header.");
            }
            PSBwaUtils.addReferenceSequencesToHeader(header, bwaArgs.referencePath, getReferenceWindowFunction());
            final JavaRDD<GATKRead> reads = readsSource.getParallelReads(path, null, null, bamPartitionSplitSize);
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
                          final JavaSparkContext ctx, SAMFileHeader header) {


        //Only retain header sequences that were aligned to.
        //This invokes an action and therefore the reads must be cached.
        reads.persist(StorageLevel.MEMORY_AND_DISK_SER());
        header = PSBwaUtils.removeUnmappedHeaderSequences(header, reads, logger);

        try {
            final String outputSuffix = isPaired ? ".paired.bam" : ".unpaired.bam";
            ReadsSparkSink.writeReads(ctx, ouputUri + outputSuffix, bwaArgs.referencePath, reads, header,
                    shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                    PSUtils.pathseqGetRecommendedNumReducers(inputBamPath, numReducers, getTargetPartitionSize()));
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(ouputUri, "writing failed", e);
        }
    }

    /**
     * Loads a bam, aligns using the given aligner, and writes to a new bam. Returns false if the input bam could not
     * be read.
     */
    private boolean alignBam(final String inputBamPath, final PSBwaAlignerSpark aligner, final boolean isPaired,
                             final JavaSparkContext ctx, final ReadsSparkSource readsSource) {

        final Tuple2<SAMFileHeader, JavaRDD<GATKRead>> loadedBam = loadBam(inputBamPath, readsSource);
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

        final PSBwaAlignerSpark aligner = new PSBwaAlignerSpark(ctx, bwaArgs);
        boolean bPairedSuccess = alignBam(inputPaired, aligner, true, ctx, readsSource);
        boolean bUnpairedSuccess = alignBam(inputUnpaired, aligner, false, ctx, readsSource);
        if (!bPairedSuccess && !bUnpairedSuccess) {
            throw new UserException.BadInput("No reads were loaded. Ensure --pairedInput and/or --unpairedInput are set and valid.");
        }
        aligner.close();
    }

}
