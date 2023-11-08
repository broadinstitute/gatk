package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSPairedUnpairedSplitterSpark;
import org.broadinstitute.hellbender.tools.spark.utils.ReadFilterSparkifier;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Realigns soft-clipped reads. Intended for use with short-read Dragen v3.7.8 BAMs/CRAMs.
 */
@CommandLineProgramProperties(
        summary = "Realigns soft-clipped reads to a given reference.",
        oneLineSummary = "Realigns soft-clipped reads to a given reference.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public final class RealignSoftClippedReads extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc="Output bam file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public GATKPath output;

    @Argument(doc="Minimum length of soft clips for realignment.",
            fullName = "min-clipped-length",
            optional = true,
            minValue = 0)
    public int minSoftClipLength = 0;

    @Argument(doc = "Retain duplicate flags on realigned reads.",
            fullName = "reset-duplicate-flag",
            optional = true)
    public boolean resetDuplicateFlag = false;

    @Advanced
    @Argument(doc = "Number of reads per partition to use for intermediate steps.",
            fullName = "alignment-reads-per-partition",
            optional = true,
            minValue = 100)
    public int intermediateReadsPerPartition = 5000;

    @Advanced
    @Argument(doc = "Estimated number of reads per partition in the input.",
            fullName = "input-reads-per-partition",
            optional = true,
            minValue = 1)
    public int inputReadsPerPartition = 200000;

    @Argument(doc = "The BWA-MEM index image file name that you've distributed to each executor",
            fullName = BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_FULL_NAME,
            shortName = BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_SHORT_NAME,
            optional = true)
    public String indexImageFile;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection(){
        return new SequenceDictionaryValidationArgumentCollection.NoValidationCollection();
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        // First traversal to find reads that are soft-clipped
        final Set<String> splitReadIds = new HashSet<>(getReads().filter(new ReadFilterSparkifier(ReadFilterLibrary.PRIMARY_LINE))
                .filter(new ReadFilterSparkifier(new ValidSoftClipReadFilter(minSoftClipLength))).map(GATKRead::getName).collect());
        final Broadcast<Set<String>> splitReadIdsBroadcast = ctx.broadcast(splitReadIds);
        final Broadcast<SAMFileHeader> headerBroadCast = ctx.broadcast(getHeaderForReads());

        // Second traversal to gather soft-clipped read pairs
        final JavaRDD<GATKRead> primarySplitReadsRdd = getReads()
                .filter(r -> splitReadIdsBroadcast.getValue().contains(r.getName()))
                .filter(r -> ReadFilterLibrary.PRIMARY_LINE.test(r))
                .map(r -> unmapRead(r, headerBroadCast.getValue()));
        // Cache the result before running the splitter, since usually only <5% of the alignments are soft-clipped
        primarySplitReadsRdd.persist(StorageLevel.DISK_ONLY());
        final long numPrimarySplitReads = primarySplitReadsRdd.count();
        final int numSplitterParitions = 1 + (int) (numPrimarySplitReads / intermediateReadsPerPartition);

        // Split into paired/unpaired soft-clipped reads
        final PSPairedUnpairedSplitterSpark splitter = new PSPairedUnpairedSplitterSpark(primarySplitReadsRdd, inputReadsPerPartition, numSplitterParitions, false);

        // Rebalance partitions using the counts
        JavaRDD<GATKRead> pairedReads = splitter.getPairedReads();
        pairedReads.persist(StorageLevel.DISK_ONLY());
        final long numPairedReads = pairedReads.count();
        final int numPairedPartitions = 1 + (int) (numPairedReads / intermediateReadsPerPartition);
        pairedReads = SparkUtils.repartitionPairedReads(pairedReads, numPairedPartitions, numPairedReads);

        JavaRDD<GATKRead> unpairedReads = splitter.getUnpairedReads();
        unpairedReads.persist(StorageLevel.DISK_ONLY());
        final long numUnpairedReads = unpairedReads.count();
        final int numUnpairedPartitions = 1 + (int) (numUnpairedReads / intermediateReadsPerPartition);
        unpairedReads = unpairedReads.repartition(numUnpairedPartitions);

        // Realign
        final SAMFileHeader alignmentHeader = getHeaderForReads().clone();
        alignmentHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final BwaSparkEngine bwa = new BwaSparkEngine(ctx, referenceArguments.getReferenceFileName(), indexImageFile, alignmentHeader, getReferenceSequenceDictionary(), !resetDuplicateFlag);

        pairedReads = bwa.alignPaired(pairedReads).map(r -> setRealigned(r));
        pairedReads = SparkUtils.sortReadsAccordingToHeader(pairedReads, getHeaderForReads(), pairedReads.getNumPartitions());
        pairedReads.persist(StorageLevel.DISK_ONLY());

        unpairedReads = bwa.alignUnpaired(unpairedReads).map(r -> setRealigned(r));
        unpairedReads = SparkUtils.sortReadsAccordingToHeader(unpairedReads, getHeaderForReads(), unpairedReads.getNumPartitions());
        unpairedReads.persist(StorageLevel.DISK_ONLY());

        // Third traversal to gather unclipped read pairs
        final JavaRDD<GATKRead> notSplitReadsRdd = getReads()
                .filter(r -> !splitReadIdsBroadcast.getValue().contains(r.getName()));

        // Merge reads
        final JavaRDD<GATKRead> mergedReads = notSplitReadsRdd.union(pairedReads.union(unpairedReads));

        try {
            // Reduce number of partitions since we previously went to ~5K reads per partition, which
            // is far too small for sharded output.
            ReadsSparkSink.writeReads(ctx, output.toString(),
                    hasReference() ? referenceArguments.getReferenceSpecifier() : null, mergedReads, getHeaderForReads(),
                    shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, mergedReads.getNumPartitions(),
                    shardedPartsDir, false, splittingIndexGranularity);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, "writing failed", e);
        }
    }

    static GATKRead unmapRead(final GATKRead read, final SAMFileHeader header) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(read.getName());
        record.setReadBases(read.getBases());
        record.setBaseQualities(read.getBaseQualities());
        record.setAttribute(SAMTag.RG, read.getReadGroup());
        record.setDuplicateReadFlag(read.isDuplicate());
        if (read.isReverseStrand()) {
            record.reverseComplement();
        }
        return SAMRecordToGATKReadAdapter.headerlessReadAdapter(record);
    }

    static GATKRead setRealigned(GATKRead read) {
        final GATKRead copy = read.copy();
        copy.setAttribute("RA", 1);
        return copy;
    }

    private static final class ValidSoftClipReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        private final int minSoftClipLength;

        public ValidSoftClipReadFilter(final int minSoftClipLength) {
            this.minSoftClipLength = minSoftClipLength;
        }

        @Override
        public boolean test(final GATKRead read) {
            return ReadFilterLibrary.PRIMARY_LINE.test(read)
                    && read.getCigarElements().stream().anyMatch(c -> c.getOperator() == CigarOperator.SOFT_CLIP
                    && c.getLength() >= this.minSoftClipLength);
        }

    }
}
