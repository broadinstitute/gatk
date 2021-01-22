package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Predicate;

/**
 * Creates discordant read pair and split read evidence files for use in the GATK-SV pipeline.
 *
 * This tool emulates the functionality of the "svtk collect-pesr" used in v1 of the GATK-SV pipeline.
 * The first output file is a tab-delimited file containing information on discordant read pairs in the
 * input cram, with the following columns:
 *
 * <ul>
 *     <li>read contig</li>
 *     <li>read start</li>
 *     <li>read strand</li>
 *     <li>mate contig</li>
 *     <li>mate start</li>
 *     <li>mate strand</li>
 *     <li>sample name</li>
 * </ul>
 *
 * Only one record is emitted for each discordant read pair, at the read in the pair with the "upstream" start
 * position according to the sequence dictionary contig ordering and coordinate.
 *
 * The second file contains the locations of all split read clippings in the input bam or cram, with the
 * following columns:
 *
 * <ul>
 *     <li>contig</li>
 *     <li>clipping position</li>
 *     <li>direction: side of the read that was clipped (either "left" or "right")</li>
 *     <li>count: the number of reads clipped at this location in this direction</li>
 *     <li>sample name</li>
 * </ul>
 */
@BetaFeature
@CommandLineProgramProperties(
        summary = "Gathers paired-end and split read evidence files for use in the GATK-SV pipeline. Output files " +
                "are a file containing the location of and orientation of read pairs marked as discordant, and a " +
                "file containing the clipping location of all soft clipped reads and the orientation of the clipping.",
        oneLineSummary = "Gathers paired-end and split read evidence files for use in the GATK-SV pipeline.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class PairedEndAndSplitReadEvidenceCollection extends ReadWalker {

    public static final String PAIRED_END_FILE_ARGUMENT_SHORT_NAME = "PE";
    public static final String PAIRED_END_FILE_ARGUMENT_LONG_NAME = "pe-file";
    public static final String SPLIT_READ_FILE_ARGUMENT_SHORT_NAME = "SR";
    public static final String SPLIT_READ_FILE_ARGUMENT_LONG_NAME = "sr-file";
    public static final String SAMPLE_NAME_ARGUMENT_LONG_NAME = "sample-name";
    public static final String COMPRESSION_LEVEL_ARGUMENT_LONG_NAME = "compression-level";

    @Argument(shortName = PAIRED_END_FILE_ARGUMENT_SHORT_NAME, fullName = PAIRED_END_FILE_ARGUMENT_LONG_NAME, doc = "Output file for paired end evidence", optional=false)
    public GATKPath peFile;

    @Argument(shortName = SPLIT_READ_FILE_ARGUMENT_SHORT_NAME, fullName = SPLIT_READ_FILE_ARGUMENT_LONG_NAME, doc = "Output file for split read evidence", optional=false)
    public GATKPath srFile;

    @Argument(fullName = SAMPLE_NAME_ARGUMENT_LONG_NAME, doc = "Sample name")
    String sampleName = null;

    @Argument(fullName = COMPRESSION_LEVEL_ARGUMENT_LONG_NAME, doc = "Output compression level")
    int compressionLevel = 4;

    final Set<String> observedDiscordantNames = new HashSet<>();
    final PriorityQueue<SplitPos> splitPosBuffer = new PriorityQueue<>(new SplitPosComparator());
    final List<DiscordantRead> discordantPairs = new ArrayList<>();

    int currentDiscordantPosition = -1;
    String currentChrom = null;

    private FeatureOutputStream<DiscordantPairEvidence> peWriter;
    private FeatureOutputStream<SplitReadEvidence> srWriter;

    private SAMSequenceDictionary sequenceDictionary;

    @Override
    public boolean requiresReads() {
        return true;
    }


    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        sequenceDictionary = getBestAvailableSequenceDictionary();
        peWriter = new FeatureOutputStream<>(peFile, new DiscordantPairEvidenceCodec(), DiscordantPairEvidenceCodec::encode, sequenceDictionary, compressionLevel);
        srWriter = new FeatureOutputStream<>(srFile, new SplitReadEvidenceCodec(), SplitReadEvidenceCodec::encode, sequenceDictionary, compressionLevel);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.MATE_UNMAPPED_AND_UNMAPPED_READ_FILTER);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);
        return readFilters;
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (isSoftClipped(read)) {
            countSplitRead(read, splitPosBuffer, srWriter);
        }

        if (! read.isProperlyPaired()) {
            reportDiscordantReadPair(read);
        }
    }

    private void reportDiscordantReadPair(final GATKRead read) {
        if (read.getStart() != currentDiscordantPosition) {
            flushDiscordantReadPairs();
            currentDiscordantPosition = read.getStart();
            observedDiscordantNames.clear();
        }

        final DiscordantRead reportableDiscordantReadPair = getReportableDiscordantReadPair(read, observedDiscordantNames,
                sequenceDictionary);
        if (reportableDiscordantReadPair != null) {
            discordantPairs.add(reportableDiscordantReadPair);
        }
    }

    @VisibleForTesting
    public DiscordantRead getReportableDiscordantReadPair(final GATKRead read, final Set<String> observedDiscordantNamesAtThisLocus,
                                                          final SAMSequenceDictionary samSequenceDictionary) {
        final int readSeqId = samSequenceDictionary.getSequenceIndex(read.getContig());
        final int mateSeqId = samSequenceDictionary.getSequenceIndex(read.getMateContig());
        if (readSeqId < mateSeqId) {
            return new DiscordantRead(read);
        } else if (readSeqId == mateSeqId) {
            if (read.getStart() < read.getMateStart()) {
                return new DiscordantRead(read);
            } else if (read.getStart() == read.getMateStart()) {
                final boolean seenBefore = observedDiscordantNamesAtThisLocus.remove(read.getName());
                if (! seenBefore) {
                    final DiscordantRead discordantRead = new DiscordantRead(read);
                    observedDiscordantNamesAtThisLocus.add(read.getName());
                    return discordantRead;
                }
            }
        }
        return null;
    }

    private void flushDiscordantReadPairs() {
        final Comparator<DiscordantRead> discReadComparator = new DiscordantReadComparator(sequenceDictionary);

        discordantPairs.sort(discReadComparator);
        discordantPairs.forEach(this::writeDiscordantPair);
        discordantPairs.clear();
    }

    private void writeDiscordantPair(final DiscordantRead r) {
        final boolean strandA = !r.isReadReverseStrand();
        final boolean strandB = !r.isMateReverseStrand();

        final DiscordantPairEvidence discordantPair = new DiscordantPairEvidence(sampleName, r.getContig(), r.getStart(), strandA, r.getMateContig(), r.getMateStart(), strandB);
        peWriter.add(discordantPair);
    }

    /**
     * Adds split read information about the current read to the counts in splitCounts. Flushes split read counts to
     * srWriter if necessary.
     */
    @VisibleForTesting
    public void countSplitRead(final GATKRead read, final PriorityQueue<SplitPos> splitCounts, final FeatureOutputStream<SplitReadEvidence> srWriter) {
        final SplitPos splitPosition = getSplitPosition(read);
        final int readStart = read.getStart();
        if (splitPosition.direction == POSITION.MIDDLE) {
            return;
        }
        if (currentChrom == null) {
            currentChrom = read.getContig();
        } else if (!currentChrom.equals(read.getContig())) {
            flushSplitCounts(splitPos -> true, srWriter, splitCounts);
            currentChrom = read.getContig();
        } else {
            flushSplitCounts(sp -> (sp.pos < readStart - 1), srWriter, splitCounts);
        }

        splitCounts.add(splitPosition);
    }

    private void flushSplitCounts(final Predicate<SplitPos> flushablePosition, final FeatureOutputStream<SplitReadEvidence> srWriter, final PriorityQueue<SplitPos> splitCounts) {

        while (splitCounts.size() > 0 && flushablePosition.test(splitCounts.peek())) {
            SplitPos pos = splitCounts.poll();
            int countAtPos = 1;
            while (splitCounts.size() > 0 && splitCounts.peek().equals(pos)) {
                countAtPos++;
                splitCounts.poll();
            }
            final SplitReadEvidence splitRead = new SplitReadEvidence(sampleName, currentChrom, pos.pos, countAtPos, pos.direction.equals(POSITION.RIGHT));
            srWriter.add(splitRead);
        }
    }

    private SplitPos getSplitPosition(GATKRead read) {
        if (read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.M) {
            final int matchLength = read.getCigar().getCigarElements().stream().filter(e -> e.getOperator().consumesReferenceBases()).mapToInt(CigarElement::getLength).sum();
            return new SplitPos(read.getStart() + matchLength, POSITION.RIGHT);
        } else if (read.getCigar().getLastCigarElement().getOperator() == CigarOperator.M) {
            return new SplitPos(read.getStart(), POSITION.LEFT);
        }

        return new SplitPos(-1, POSITION.MIDDLE);
    }

    private boolean isSoftClipped(final GATKRead read) {
        final CigarOperator firstOperator = read.getCigar().getFirstCigarElement().getOperator();
        final CigarOperator lastOperator = read.getCigar().getLastCigarElement().getOperator();
        return (firstOperator == CigarOperator.SOFT_CLIP && lastOperator != CigarOperator.SOFT_CLIP) ||
                (firstOperator != CigarOperator.SOFT_CLIP && lastOperator == CigarOperator.SOFT_CLIP);
    }

    @Override
    public Object onTraversalSuccess() {
        flushSplitCounts(splitPos -> true, srWriter, splitPosBuffer);
        flushDiscordantReadPairs();
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (peWriter != null) {
            peWriter.close();
        }
        if (srWriter != null) {
            srWriter.close();
        }
    }

    enum POSITION {
        LEFT ("left"),
        MIDDLE ("middle"),
        RIGHT ("right");

        private String description;

        POSITION(final String description) {
            this.description = description;
        }

        public String getDescription() {
            return description;
        }
    }

    @VisibleForTesting final static class DiscordantRead {
        private boolean readReverseStrand;
        private boolean mateReverseStrand;
        private String contig;
        private int start;
        private String mateContig;
        private int mateStart;
        private String name;

        public DiscordantRead(final GATKRead read) {
            this.readReverseStrand = read.isReverseStrand();
            this.mateReverseStrand = read.mateIsReverseStrand();
            this.contig = read.getContig();
            this.start = read.getStart();
            this.mateContig = read.getMateContig();
            this.mateStart = read.getMateStart();
            this.name = read.getName();
        }

        public boolean isReadReverseStrand() {
            return readReverseStrand;
        }

        public void setReadReverseStrand(final boolean readReverseStrand) {
            this.readReverseStrand = readReverseStrand;
        }

        public boolean isMateReverseStrand() {
            return mateReverseStrand;
        }

        public void setMateReverseStrand(final boolean mateReverseStrand) {
            this.mateReverseStrand = mateReverseStrand;
        }

        public String getContig() {
            return contig;
        }

        public void setContig(final String contig) {
            this.contig = contig;
        }

        public int getStart() {
            return start;
        }

        public void setStart(final int start) {
            this.start = start;
        }

        public String getMateContig() {
            return mateContig;
        }

        public void setMateContig(final String mateContig) {
            this.mateContig = mateContig;
        }

        public int getMateStart() {
            return mateStart;
        }

        public void setMateStart(final int mateStart) {
            this.mateStart = mateStart;
        }

        public String getName() {
            return name;
        }

        public void setName(final String name) {
            this.name = name;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final DiscordantRead that = (DiscordantRead) o;

            if (readReverseStrand != that.readReverseStrand) return false;
            if (mateReverseStrand != that.mateReverseStrand) return false;
            if (start != that.start) return false;
            if (mateStart != that.mateStart) return false;
            if (contig != null ? !contig.equals(that.contig) : that.contig != null) return false;
            if (mateContig != null ? !mateContig.equals(that.mateContig) : that.mateContig != null) return false;
            return name != null ? name.equals(that.name) : that.name == null;
        }

        @Override
        public int hashCode() {
            int result = (readReverseStrand ? 1 : 0);
            result = 31 * result + (mateReverseStrand ? 1 : 0);
            result = 31 * result + (contig != null ? contig.hashCode() : 0);
            result = 31 * result + start;
            result = 31 * result + (mateContig != null ? mateContig.hashCode() : 0);
            result = 31 * result + mateStart;
            result = 31 * result + (name != null ? name.hashCode() : 0);
            return result;
        }
    }

    @VisibleForTesting final static class SplitPos {
        public POSITION direction;
        public int pos;

        public SplitPos(final int start, final POSITION direction) {
            this.pos = start;
            this.direction = direction;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final SplitPos splitPos = (SplitPos) o;

            if (pos != splitPos.pos) return false;
            return direction.ordinal() == splitPos.direction.ordinal();
        }

        @Override
        public int hashCode() {
            int result = direction != null ? direction.ordinal() : 0;
            result = 31 * result + pos;
            return result;
        }
    }

    @VisibleForTesting final static class SplitPosComparator implements Comparator<SplitPos> {
        @Override
        public int compare(final SplitPos o1, final SplitPos o2) {
            if (o1.pos != o2.pos) {
                return Integer.compare(o1.pos, o2.pos);
            } else {
                return o1.direction.compareTo(o2.direction);
            }
        }
    }

    @VisibleForTesting final static class DiscordantReadComparator implements Comparator<DiscordantRead> {

        private final Comparator<DiscordantRead> internalComparator;

        public DiscordantReadComparator(final SAMSequenceDictionary sequenceDictionary) {
            internalComparator = Comparator.comparing((DiscordantRead r) -> sequenceDictionary.getSequenceIndex(r.getContig()))
                    .thenComparing(DiscordantRead::getStart)
                    .thenComparing(DiscordantRead::isReadReverseStrand)
                    .thenComparing((DiscordantRead r) -> sequenceDictionary.getSequenceIndex(r.getMateContig()))
                    .thenComparing(DiscordantRead::getMateStart)
                    .thenComparing(DiscordantRead::isMateReverseStrand);

        }

        @Override
        public int compare(final DiscordantRead o1, final DiscordantRead o2) {
            return internalComparator.compare(o1, o2);
        }
    }
}
