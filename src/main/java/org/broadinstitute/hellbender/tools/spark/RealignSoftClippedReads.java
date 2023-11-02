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
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSPairedUnpairedSplitterSpark;
import org.broadinstitute.hellbender.tools.spark.utils.ReadFilterSparkifier;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.IOException;
import java.util.*;

/**
 * Realigns soft-clipped reads.
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
    @Argument(doc = "Number of reads per partition to use for realignment.",
            fullName = "alignment-reads-per-partition",
            optional = true,
            minValue = 100)
    public int alignmentReadsPerPartition = 5000;

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
        primarySplitReadsRdd.persist(StorageLevel.MEMORY_AND_DISK());

        // Third traversal to gather unclipped read pairs
        final JavaRDD<GATKRead> notSplitReadsRdd = getReads()
                .filter(r -> !splitReadIdsBroadcast.getValue().contains(r.getName()));

        // Split into paired/unpaired soft-clipped reads
        final PSPairedUnpairedSplitterSpark splitter = new PSPairedUnpairedSplitterSpark(primarySplitReadsRdd, inputReadsPerPartition, false);
        JavaRDD<GATKRead> pairedReads = splitter.getPairedReads();
        JavaRDD<GATKRead> unpairedReads = splitter.getUnpairedReads();

        // Counting forces an action on the RDDs to guarantee we're done with the Bwa image and kmer filter
        final long numPairedReads = pairedReads.count();
        final long numUnpairedReads = unpairedReads.count();

        // Rebalance partitions using the counts
        final int numPairedPartitions = 1 + (int) (numPairedReads / alignmentReadsPerPartition);
        final int numUnpairedPartitions = 1 + (int) (numUnpairedReads / alignmentReadsPerPartition);
        pairedReads = repartitionPairedReads(pairedReads, numPairedPartitions, numPairedReads);
        unpairedReads = unpairedReads.repartition(numUnpairedPartitions);

        // Cache again after the shuffle
        pairedReads.persist(StorageLevel.MEMORY_AND_DISK());
        unpairedReads.persist(StorageLevel.MEMORY_AND_DISK());

        // Realign
        final SAMFileHeader alignmentHeader = getHeaderForReads().clone();
        alignmentHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final BwaSparkEngine bwa = new BwaSparkEngine(ctx, referenceArguments.getReferenceFileName(), indexImageFile, alignmentHeader, getReferenceSequenceDictionary(), !resetDuplicateFlag);
        pairedReads = bwa.alignPaired(pairedReads).map(r -> setRealigned(r));
        unpairedReads = bwa.alignUnpaired(unpairedReads).map(r -> setRealigned(r));

        // Cache this since it's expensive
        pairedReads.persist(StorageLevel.MEMORY_AND_DISK());
        unpairedReads.persist(StorageLevel.MEMORY_AND_DISK());

        // Merge reads
        final JavaRDD<GATKRead> mergedReads = notSplitReadsRdd.union(pairedReads.union(unpairedReads));

        try {
            // Reduce number of partitions since we previously went to ~5K reads per partition, which
            // is far too small for sharded output.
            ReadsSparkSink.writeReads(ctx, output.toString(),
                    hasReference() ? referenceArguments.getReferenceSpecifier() : null, mergedReads, getHeaderForReads(),
                    shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE, mergedReads.getNumPartitions(),
                    shardedPartsDir, true, splittingIndexGranularity);
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

    /**
     * Reduces number of partitions of paired reads, keeping pairs together.
     */
    private static JavaRDD<GATKRead> repartitionPairedReads(final JavaRDD<GATKRead> pairedReads, final int alignmentPartitions, final long numReads) {
        final int readsPerPartition = 1 + (int) (numReads / alignmentPartitions);
        return pairedReads.mapPartitions(iter -> pairPartitionReads(iter, readsPerPartition))
                .repartition(alignmentPartitions)
                .flatMap(List::iterator);
    }
    /**
     * Maps partition of paired reads to a partition of Lists containing each pair. Assumes pairs are adjacent.
     */
    private static Iterator<List<GATKRead>> pairPartitionReads(final Iterator<GATKRead> iter, final int readsPerPartition) {
        final ArrayList<List<GATKRead>> readPairs = new ArrayList<>(readsPerPartition / 2);
        while (iter.hasNext()) {
            final List<GATKRead> list = new ArrayList<>(2);
            list.add(iter.next());
            if (!iter.hasNext()) throw new GATKException("Odd number of read pairs in paired reads partition");
            list.add(iter.next());
            if (!list.get(0).getName().equals(list.get(1).getName())) throw new GATKException("Pair did not have the same name in a paired reads partition");
            readPairs.add(list);
        }
        return readPairs.iterator();
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
