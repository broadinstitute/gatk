package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Streams;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Realigns soft-clipped reads using BWA.
 *
 * <p>
 * This tool consumes coordinate-sorted alignments and realigns all soft-clipped primary reads and their mates.
 * The new alignments are then efficiently merged with the non-clipped reads. The tool works in the following steps:
 * </p>
 *
 * <ol>
 *     <li>Traverse input to identify soft-clipped primary alignments</li>
 *     <li>Traverse input a second time to pull out soft-clipped primary alignments and their mates</li>
 *     <li>Realign the soft-clipped reads/read pairs with BWA</li>
 *     <li>Sort the realignments and merge them with the original unclipped read pairs</li>
 * </ol>
 *
 * <p>
 * Input reads may be paired or unpaired or a mixture. For paired reads, their names must be identical
 * (no "/1" or "/2" suffixes; SAM flags are ignored). Also, the BWA index image must match the reference of the input
 * alignments. In other words, the realignment must be to the same reference as the original reads.
 * </p>
 *
 * <p>
 * Default BWA alignment parameters are used, with the exception of the addition of the -Y command line argument,
 * which enables soft-clipping of supplementary alignments.
 * </p>
 *
 * <p>
 * Unclipped reads (whose mates are also not soft-clipped) are emitted in the output exactly as they appear in the
 * input. That is, SAM flags and tags are not changed. For the realigned reads, the only tags retained are read
 * groups (RG). The "--keep-duplicate-flag" option may be enabled to retain the duplicate read SAM flag
 * for realigned reads. Reads that were realigned can be identified by an "RA" tag set to "1" in the output.
 * </p>
 *
 * <p>
 * In order to minimize the memory footprint, this tool caches reads to disk in BAM format. Therefore, it is
 * recommended to have free disk space equal to 2-3x of the input BAM size on the java.io.tmpdir volume.
 * </p>
 *
 * <p>
 * It is highly recommended to utilize BWA multi-threading with the "--bwa-threads" parameter. Runtime is also
 * substantially reduced (~2x) using BAM input/output compared to CRAM.
 * </p>
 *
 * <p>
 * Use of interval subsetting (with -L/-XL) is supported. However, use caution because mates outside of the selected
 * intervals will NOT be included. Reads with mates lying outside the intervals will be aligned in unpaired mode.
 * </p>
 *
 * <p>
 * The reference (provided with -R) is strictly required when handling CRAM files.
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 *     <li> BWA index image file of the reference (created with BwaMemIndexImageCreator) </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * Realign reads with soft-clips of at least 10 bases, utilizing 4 threads for alignment.
 * <pre>
 * gatk RealignSoftClippedReads \
 *   -I input.bam \
 *   --bwa-mem-index-image ref.fasta.img \
 *   --bwa-threads 4 \
 *   --min-clipped-length 10 \
 *   -O output.bam
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Realigns soft-clipped reads to a given reference using BWA.",
        oneLineSummary = "Realigns soft-clipped reads to a given reference using BWA.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public final class RealignSoftClippedReads extends MultiplePassReadWalker {

    public static final String MIN_SOFT_CLIP_LENGTH_LONG_NAME = "min-clipped-length";

    @Argument(doc="Output bam file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public GATKPath output;

    @Argument(doc="Minimum length of soft clips required for realignment. Setting to 0 will realign all reads.",
            fullName = MIN_SOFT_CLIP_LENGTH_LONG_NAME,
            optional = true,
            minValue = 0)
    public int minSoftClipLength = 1;

    @ArgumentCollection
    public SubsettingRealignmentArgumentCollection args = new SubsettingRealignmentArgumentCollection();

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public void traverseReads() {
        logger.info("Identifying soft-clipped reads...");
        final Set<String> softclippedReadNames = new HashSet<>();
        final Set<String> pairedClippedReadsWithDistantMates = new HashSet<>();
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) -> {
            checkIfClipped(read, softclippedReadNames, minSoftClipLength, pairedClippedReadsWithDistantMates);
        });
        logger.info("Found " + softclippedReadNames.size() + " soft-clipped reads / read pairs.");
        logger.info("Found " + pairedClippedReadsWithDistantMates.size() + " distant mates.");

        try (final SubsettingRealignmentEngine engine = new SubsettingRealignmentEngine(args.indexImage.toString(),
                getHeaderForReads(), args.bufferSize, args.bwaThreads, args.keepDuplicateFlag);
             final SAMFileGATKReadWriter writer = createSAMWriter(output, true)) {
            final OverlapDetector<SimpleInterval> traversalIntervals = intervalArgumentCollection.intervalsSpecified() ? OverlapDetector.create(intervalArgumentCollection.getIntervals(dictionary)) : null;
            final List<MateInfo> distantMates = new ArrayList<>();
            logger.info("Subsetting soft-clipped reads and close mates...");
            // Add reads
            forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) -> {
                // Clear realignment tags in case the input was previously realigned
                read.clearAttribute(SubsettingRealignmentEngine.REALIGNED_READ_TAG);
                // Submit the read to the engine, selecting soft-clipped reads and their mates
                engine.addRead(read, r -> softclippedReadNames.contains(r.getName()));
                if (pairedClippedReadsWithDistantMates.contains(read.getName())) {
                    addMateLocus(read, distantMates, traversalIntervals);
                }
            });

            // Add distant mates
            SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());
            if (hasReference()) { // pass in reference if available, because CRAM files need it
                factory = factory.referenceSequence(referenceArguments.getReferencePath());
            }
            final ReadsDataSource readSource = new ReadsPathDataSource(readArguments.getReadPaths(), readArguments.getReadIndexPaths(), factory, 1,
                    (cloudIndexPrefetchBuffer < 0 ? 1 : cloudIndexPrefetchBuffer));
            // Sort so we traverse in order
            Collections.sort(distantMates, IntervalUtils.getDictionaryOrderComparator(dictionary));
            logger.info("Adding " + distantMates.size() + " distant mates...");
            for (final MateInfo mate : distantMates) {
                logger.info(mate.name + " " + mate.contig + ":" + mate.start);
                final List<GATKRead> mates = Streams.stream(readSource.query(new SimpleInterval(mate))).filter(r -> validMate(r, mate)).collect(Collectors.toList());
                if (mates.size() != 1) {
                    throw new GATKException("Valid mate of " + mate.name + " not found at " + mate.contig + ":" + mate.start);
                }
                engine.addDistantMate(mates.get(0));
            }

            logger.info("Found " + engine.getSelectedReadsCount() + " soft-clipped reads and mates.");
            logger.info("Found " + engine.getNonselectedReadsCount() + " non-clipped reads.");
            logger.info("Realigning and merging reads...");
            try (final MergingIterator<GATKRead> iter = engine.alignAndMerge()) {
                while (iter.hasNext()) {
                    writer.addRead(iter.next());
                }
            }
            logger.info("Realigned " + engine.getPairedAlignmentReadsCount() + " paired and " +
                    engine.getUnpairedAlignmentReadsCount() + " unpaired reads.");
        }
    }

    /**
     * Gather names of soft-clipped reads
     */
    @VisibleForTesting
    static void checkIfClipped(final GATKRead read, final Set<String> softclippedReadNames, final int minSoftClipLength, final Set<String> pairedClippedReadsWithDistantMates) {
        if (isValidSoftClip(read, minSoftClipLength)) {
            softclippedReadNames.add(read.getName());
            evaluateDistantMate(read, pairedClippedReadsWithDistantMates);
        }
    }

    static void evaluateDistantMate(final GATKRead read, final Set<String> pairedClippedReadsWithDistantMates) {
        if (ReadFilterLibrary.PRIMARY_LINE.test(read) && read.isPaired()) {
            final String name = read.getName();
            if (pairedClippedReadsWithDistantMates.contains(name)) {
                pairedClippedReadsWithDistantMates.remove(name);
            } else {
                pairedClippedReadsWithDistantMates.add(name);
            }
        }
    }

    static void addMateLocus(final GATKRead read, final List<MateInfo> loci, final OverlapDetector<SimpleInterval> traversalIntervals) {
        // Note we exclude unmapped mates, since they are impossible to pull out efficiently
        if (ReadFilterLibrary.PRIMARY_LINE.test(read) && read.isPaired() && read.getMateContig() != null) {
            final SimpleInterval mateLocus = new SimpleInterval(read.getMateContig(), read.getMateStart(), read.getMateStart() + 1);
            if (traversalIntervals == null || !traversalIntervals.overlapsAny(mateLocus)) {
                loci.add(new MateInfo(read.getMateContig(), read.getMateStart(), !read.isFirstOfPair(), read.getName()));
            }
        }
    }
    static boolean validMate(final GATKRead read, final MateInfo mate) {
        if (!(ReadFilterLibrary.PRIMARY_LINE.test(read) && read.isPaired() && read.getName().equals(mate.name) && read.isFirstOfPair() == mate.firstRead)) {
            return false;
        }
        if (read.getContig() == null) {
            return mate.contig == null;
        } else {
            return read.getContig().equals(mate.contig) && read.getStart() == mate.start;
        }
    }

    /**
     * Returns whether the read is non-secondary/non-supplementary and sufficiently soft-clipped
     */
    @VisibleForTesting
    static boolean isValidSoftClip(final GATKRead read, final int minSoftClipLength) {
        return ReadFilterLibrary.PRIMARY_LINE.test(read) && hasMinSoftClip(read, minSoftClipLength);
    }

    private static boolean hasMinSoftClip(final GATKRead read, final int minSoftClipLength) {
        if (minSoftClipLength == 0) {
            return true;
        }
        for (final CigarElement e : read.getCigarElements()) {
            if (e.getOperator() == CigarOperator.SOFT_CLIP && e.getLength() >= minSoftClipLength) {
                return true;
            }
        }
        return false;
    }

    private static final class MateInfo implements Locatable {
        final String contig;
        final int start;
        final boolean firstRead;
        final String name;

        public MateInfo(final String contig, final int start, final boolean firstRead, final String name) {
            this.contig = contig;
            this.start = start;
            this.firstRead = firstRead;
            this.name = name;
        }

        @Override
        public String getContig() {
            return contig;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return start;
        }
    }
}
