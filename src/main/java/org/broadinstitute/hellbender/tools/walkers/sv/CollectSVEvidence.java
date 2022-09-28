package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SiteDepth;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.function.Predicate;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.isBaseInsideAdaptor;

/**
 * Creates discordant read pair, split read evidence, site depth, and read depth files for use in the GATK-SV pipeline.
 * This tool emulates the functionality of the "svtk collect-pesr" used in v1 of the GATK-SV pipeline.
 *
 * The first output file, which should be named "*.pe.txt" or "*.pe.txt.gz" is a tab-delimited file
 * containing information on discordant read pairs in the input cram, with the following columns:
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
 * The second output file, which should be named "*.sr.txt" or "*.sr.txt.gz" contains the locations
 * of all split read clippings in the input bam or cram, with the following columns:
 *
 * <ul>
 *     <li>contig</li>
 *     <li>clipping position</li>
 *     <li>direction: side of the read that was clipped (either "left" or "right")</li>
 *     <li>count: the number of reads clipped at this location in this direction</li>
 *     <li>sample name</li>
 * </ul>
 *
 * The third output file, which should be named "*.sd.txt" or "*.sd.txt.gz" specifies site depth counts:
 * For each locus specified in an input VCF as a simple, biallelic SNP, it gives a count, for each
 * base call, of the number of reads that cover that locus.
 * It has the following columns:
 *
 * <ul>
 *     <li>contig</li>
 *     <li>position</li>
 *     <li>sampleName</li>
 *     <li>A observations</li>
 *     <li>C observations</li>
 *     <li>G observations</li>
 *     <li>T observations</li>
 * </ul>
 *
 * The fourth output file, which should be named "*.rd.txt" or "*.rd.txt.gz" specifies read depths:
 * For each interval specified by a 3-column, tab delimited input file, the number of reads that
 * start in that interval are reported.
 * It has the following columns:
 *
 * <ul>
 *     <li>contig</li>
 *     <li>starting position</li>
 *     <li>ending position</li>
 *     <li>read count</li>
 * </ul>
 *
 * Note: when only collecting RD evidence, users should consider providing the same interval list
 * with -L as --depth-evidence-intervals in order to avoid processing unused reads outside the intervals.
 *
 * Each of these output files may also be written as a block-compressed interval file, rather than
 * as a tab-delimited text file by specifying an output file name that ends with ".bci" rather than
 * ".txt".  These files are self-indexing, and contain complete header information including sample
 * name(s) and a dictionary for the contigs.
 */
@BetaFeature
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Gathers paired-end and split read evidence files for use in the GATK-SV pipeline. Output files " +
                "are a file containing the location of and orientation of read pairs marked as discordant, and a " +
                "file containing the clipping location of all soft clipped reads and the orientation of the clipping.",
        oneLineSummary = "Gathers paired-end and split read evidence files for use in the GATK-SV pipeline.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class CollectSVEvidence extends ReadWalker {

    public static final String PAIRED_END_FILE_ARGUMENT_SHORT_NAME = "PE";
    public static final String PAIRED_END_FILE_ARGUMENT_LONG_NAME = "pe-file";
    public static final String SPLIT_READ_FILE_ARGUMENT_SHORT_NAME = "SR";
    public static final String SPLIT_READ_FILE_ARGUMENT_LONG_NAME = "sr-file";
    public static final String SITE_DEPTH_OUTPUT_ARGUMENT_SHORT_NAME = "SD";
    public static final String SITE_DEPTH_OUTPUT_ARGUMENT_LONG_NAME = "sd-file";
    public static final String SITE_DEPTH_INPUT_ARGUMENT_SHORT_NAME = "F";
    public static final String SITE_DEPTH_INPUT_ARGUMENT_LONG_NAME = "site-depth-locs-vcf";
    public static final String DEPTH_EVIDENCE_OUTPUT_FILE_ARGUMENT_SHORT_NAME = "RD";
    public static final String DEPTH_EVIDENCE_OUTPUT_FILE_ARGUMENT_LONG_NAME = "depth-evidence-file";
    public static final String DEPTH_EVIDENCE_SUMMARY_FILE_ARGUMENT_SHORT_NAME = "DS";
    public static final String DEPTH_EVIDENCE_SUMMARY_FILE_ARGUMENT_LONG_NAME = "depth-summary-file";
    public static final String DEPTH_EVIDENCE_INTERVALS_INPUT_FILE_ARGUMENT_SHORT_NAME = "DI";
    public static final String DEPTH_EVIDENCE_INTERVALS_INPUT_FILE_ARGUMENT_LONG_NAME = "depth-evidence-intervals";
    public static final String MIN_DEPTH_EVIDENCE_MAPQ_ARGUMENT_NAME = "depth-evidence-min-mapq";
    public static final String MIN_SITE_DEPTH_MAPQ_ARGUMENT_NAME = "site-depth-min-mapq";
    public static final String MIN_SITE_DEPTH_BASEQ_ARGUMENT_NAME = "site-depth-min-baseq";
    public static final String SAMPLE_NAME_ARGUMENT_LONG_NAME = "sample-name";
    public static final String COMPRESSION_LEVEL_ARGUMENT_LONG_NAME = "compression-level";

    @Argument(shortName = PAIRED_END_FILE_ARGUMENT_SHORT_NAME,
            fullName = PAIRED_END_FILE_ARGUMENT_LONG_NAME, doc = "Output file for paired end evidence",
            optional=true)
    public GATKPath peFile;

    @Argument(shortName = SPLIT_READ_FILE_ARGUMENT_SHORT_NAME,
            fullName = SPLIT_READ_FILE_ARGUMENT_LONG_NAME, doc = "Output file for split read evidence",
            optional=true)
    public GATKPath srFile;

    @Argument(shortName = SITE_DEPTH_OUTPUT_ARGUMENT_SHORT_NAME,
            fullName = SITE_DEPTH_OUTPUT_ARGUMENT_LONG_NAME,
            doc = "Output file for site depth counts",
            optional = true)
    public GATKPath siteDepthOutputFilename;

    @Argument(shortName = SITE_DEPTH_INPUT_ARGUMENT_SHORT_NAME,
            fullName = SITE_DEPTH_INPUT_ARGUMENT_LONG_NAME,
            doc = "Input VCF of SNPs marking loci for site depth counts",
            optional = true)
    public GATKPath siteDepthInputFilename;

    @Argument(shortName = DEPTH_EVIDENCE_OUTPUT_FILE_ARGUMENT_SHORT_NAME,
            fullName = DEPTH_EVIDENCE_OUTPUT_FILE_ARGUMENT_LONG_NAME,
            doc = "Output file for depth evidence",
            optional = true)
    public GATKPath depthEvidenceOutputFilename;

    @Argument(shortName = DEPTH_EVIDENCE_SUMMARY_FILE_ARGUMENT_SHORT_NAME,
            fullName = DEPTH_EVIDENCE_SUMMARY_FILE_ARGUMENT_LONG_NAME,
            doc = "Output file for depth evidence summary statistics",
            optional = true)
    public GATKPath depthEvidenceSummaryFilename;

    @Argument(shortName = DEPTH_EVIDENCE_INTERVALS_INPUT_FILE_ARGUMENT_SHORT_NAME,
            fullName = DEPTH_EVIDENCE_INTERVALS_INPUT_FILE_ARGUMENT_LONG_NAME,
            doc = "Input feature file specifying intervals where depth evidence will be gathered",
            optional = true)
    public GATKPath depthEvidenceInputFilename;

    @Argument(fullName = MIN_DEPTH_EVIDENCE_MAPQ_ARGUMENT_NAME,
            doc = "minimum mapping quality for read to be counted as depth evidence",
            optional = true)
    public int minDepthEvidenceMapQ = 0;

    @Argument(fullName = MIN_SITE_DEPTH_MAPQ_ARGUMENT_NAME,
            doc = "minimum mapping quality for read to be counted toward site depth",
            optional = true)
    public int minMapQ = 30;

    @Argument(fullName = MIN_SITE_DEPTH_BASEQ_ARGUMENT_NAME,
            doc = "minimum base call quality for SNP to be counted toward site depth",
            optional = true)
    public int minQ = 20;

    @Argument(fullName = SAMPLE_NAME_ARGUMENT_LONG_NAME, doc = "Sample name")
    String sampleName = null;

    @Argument(fullName = COMPRESSION_LEVEL_ARGUMENT_LONG_NAME, doc = "Output compression level")
    int compressionLevel = 4;

    final Set<String> observedDiscordantNames = new HashSet<>();
    final PriorityQueue<SplitPos> splitPosBuffer = new PriorityQueue<>(new SplitPosComparator());
    final List<DiscordantRead> discordantPairs = new ArrayList<>();

    int currentDiscordantPosition = -1;
    String currentChrom = null;

    private FeatureSink<DiscordantPairEvidence> peWriter;
    private FeatureSink<SplitReadEvidence> srWriter;
    private SiteDepthCounter siteDepthCounter;
    private DepthEvidenceCollector depthEvidenceCollector;

    private SAMSequenceDictionary sequenceDictionary;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        sequenceDictionary = getBestAvailableSequenceDictionary();
        peWriter = createPEWriter();
        srWriter = createSRWriter();
        siteDepthCounter = createSiteDepthCounter();
        depthEvidenceCollector = createDepthEvidenceCollector();
        if ( peWriter == null && srWriter == null &&
                siteDepthCounter == null && depthEvidenceCollector == null ) {
            throw new UserException("You must supply at least one output file: PE, SR, SD, or RD");
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        return readFilters;
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if ( !(read.isPaired() && read.mateIsUnmapped()) &&
                !read.isSupplementaryAlignment() &&
                !read.isSecondaryAlignment() ) {
            if ( srWriter != null && isSoftClipped(read) ) {
                countSplitRead(read, splitPosBuffer, srWriter);
            }

            if ( peWriter != null && !read.isProperlyPaired() ) {
                reportDiscordantReadPair(read);
            }
        }

        if ( siteDepthCounter != null ) {
            siteDepthCounter.apply(read);
        }
        if ( depthEvidenceCollector != null ) {
            depthEvidenceCollector.apply(read);
        }
    }

    private FeatureSink<DiscordantPairEvidence> createPEWriter() {
        if ( peFile == null ) {
            return null;
        }
        final String peFilename = peFile.toPath().toString();
        final List<String> sampleNames = Collections.singletonList(sampleName);
        final DiscordantPairEvidenceCodec peCodec = new DiscordantPairEvidenceCodec();
        final DiscordantPairEvidenceBCICodec peBCICodec = new DiscordantPairEvidenceBCICodec();
        if ( peBCICodec.canDecode(peFilename) ) {
            return peBCICodec.makeSink(peFile, sequenceDictionary, sampleNames, compressionLevel);
        }
        if ( !peCodec.canDecode(peFilename) ) {
            throw new UserException("Attempting to write discordant pair evidence to a file that " +
                    "can't be read as discordant pair evidence: " + peFilename + ".  The file " +
                    "name should end with \".pe.txt\", \".pe.txt.gz\", or \".pe.bci\".");
        }
        return peCodec.makeSink(peFile, sequenceDictionary, sampleNames, compressionLevel);
    }

    private FeatureSink<SplitReadEvidence> createSRWriter() {
        if ( srFile == null ) {
            return null;
        }
        final String srFilename = srFile.toPath().toString();
        final List<String> sampleNames = Collections.singletonList(sampleName);
        final SplitReadEvidenceCodec srCodec = new SplitReadEvidenceCodec();
        final SplitReadEvidenceBCICodec srBCICodec = new SplitReadEvidenceBCICodec();
        if ( srBCICodec.canDecode(srFilename) ) {
            return srBCICodec.makeSink(srFile, sequenceDictionary, sampleNames, compressionLevel);
        }
        if ( !srCodec.canDecode(srFilename) ) {
            throw new UserException("Attempting to write split read evidence to a file that " +
                    "can't be read as split read evidence: " + srFilename + ".  The file " +
                    "name should end with \".sr.txt\", \".sr.txt.gz\", or \".sr.bci\".");
        }
        return srCodec.makeSink(srFile, sequenceDictionary, sampleNames, compressionLevel);
    }

    private SiteDepthCounter createSiteDepthCounter() {
        if ( siteDepthInputFilename != null && siteDepthOutputFilename != null ) {
            return new SiteDepthCounter(sequenceDictionary, sampleName, compressionLevel,
                    siteDepthInputFilename, siteDepthOutputFilename,
                                        minMapQ, minQ);
        }
        if ( siteDepthInputFilename != null ) {
            throw new UserException("Having specified a " + SITE_DEPTH_INPUT_ARGUMENT_LONG_NAME +
                    " input, you must also supply an " + SITE_DEPTH_OUTPUT_ARGUMENT_LONG_NAME +
                    " for output.");
        }
        if ( siteDepthOutputFilename != null ) {
            throw new UserException("Having specified an " + SITE_DEPTH_OUTPUT_ARGUMENT_LONG_NAME +
                    " for output, you must also supply a " + SITE_DEPTH_INPUT_ARGUMENT_LONG_NAME +
                    " as input.");
        }
        return null;
    }

    private DepthEvidenceCollector createDepthEvidenceCollector() {
        if ( depthEvidenceInputFilename != null && depthEvidenceOutputFilename != null ) {
            return new DepthEvidenceCollector(sequenceDictionary, sampleName, compressionLevel,
                                            depthEvidenceInputFilename, depthEvidenceOutputFilename,
                                            minDepthEvidenceMapQ);
        }
        if ( depthEvidenceInputFilename != null ) {
            throw new UserException("Having specified an depth-evidence-intervals input, " +
                    "you must also supply a depth-evidence-file for output.");
        }
        if ( depthEvidenceOutputFilename != null ) {
            throw new UserException("Having specified a depth-evidence-file for output, " +
                    "you must also supply depth-evidence-intervals as input.");
        }
        return null;
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
        peWriter.write(new DiscordantPairEvidence(sampleName,
                            r.getContig(), r.getStart(), !r.isReadReverseStrand(),
                            r.getMateContig(), r.getMateStart(), !r.isMateReverseStrand()));
    }

    /**
     * Adds split read information about the current read to the counts in splitCounts. Flushes split read counts to
     * srWriter if necessary.
     */
    @VisibleForTesting
    public void countSplitRead(final GATKRead read,
                               final PriorityQueue<SplitPos> splitCounts,
                               final FeatureSink<SplitReadEvidence> srWriter ) {
        final SplitPos splitPosition = getSplitPosition(read);
        final int readStart = read.getStart();
        if (splitPosition.direction == POSITION.MIDDLE) {
            return;
        }
        if (currentChrom == null) {
            currentChrom = read.getContig();
        } else if (!currentChrom.equals(read.getContig())) {
            flushSplitCounts(splitPos -> true, splitCounts, srWriter);
            currentChrom = read.getContig();
        } else {
            flushSplitCounts(sp -> (sp.pos < readStart - 1), splitCounts, srWriter);
        }

        splitCounts.add(splitPosition);
    }

    private void flushSplitCounts(final Predicate<SplitPos> flushablePosition,
                                  final PriorityQueue<SplitPos> splitCounts,
                                  final FeatureSink<SplitReadEvidence> srWriter) {

        while (splitCounts.size() > 0 && flushablePosition.test(splitCounts.peek())) {
            final SplitPos pos = splitCounts.poll();
            int countAtPos = 1;
            while (splitCounts.size() > 0 && splitCounts.peek().equals(pos)) {
                countAtPos++;
                splitCounts.poll();
            }
            final SplitReadEvidence splitRead = new SplitReadEvidence(sampleName, currentChrom, pos.pos, countAtPos, pos.direction.equals(POSITION.RIGHT));
            srWriter.write(splitRead);
        }
    }

    private SplitPos getSplitPosition( final GATKRead read ) {
        if (read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.M) {
            final int matchLength = read.getCigar().getCigarElements().stream().filter(e -> e.getOperator().consumesReferenceBases()).mapToInt(CigarElement::getLength).sum();
            return new SplitPos(read.getStart() + matchLength, POSITION.RIGHT);
        } else if (read.getCigar().getLastCigarElement().getOperator() == CigarOperator.M) {
            return new SplitPos(read.getStart(), POSITION.LEFT);
        }

        return new SplitPos(-1, POSITION.MIDDLE);
    }

    private boolean isSoftClipped( final GATKRead read ) {
        final CigarOperator firstOperator = read.getCigar().getFirstCigarElement().getOperator();
        final CigarOperator lastOperator = read.getCigar().getLastCigarElement().getOperator();
        return (firstOperator == CigarOperator.SOFT_CLIP && lastOperator != CigarOperator.SOFT_CLIP) ||
                (firstOperator != CigarOperator.SOFT_CLIP && lastOperator == CigarOperator.SOFT_CLIP);
    }

    @Override
    public Object onTraversalSuccess() {
        flushSplitCounts(splitPos -> true, splitPosBuffer, srWriter);
        flushDiscordantReadPairs();
        if ( siteDepthCounter != null ) {
            siteDepthCounter.close();
        }
        if ( depthEvidenceCollector != null ) {
            depthEvidenceCollector.close();
            if ( depthEvidenceSummaryFilename != null ) {
                depthEvidenceCollector.reportSummaryStats(depthEvidenceSummaryFilename, sampleName);
            }
        }
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if ( peWriter != null ) {
            peWriter.close();
        }
        if ( srWriter != null ) {
            srWriter.close();
        }
    }

    enum POSITION {
        LEFT ("left"),
        MIDDLE ("middle"),
        RIGHT ("right");

        private final String description;

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
            if ( !Objects.equals(contig, that.contig) ) return false;
            if ( !Objects.equals(mateContig, that.mateContig) ) return false;
            return Objects.equals(name, that.name);
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

    /**
     * Compare a locus to an interval and indicates whether the locus is
     * upstream (returns -1), within (returns 0), or downstream (returns 1) of an interval
     */
    @VisibleForTesting
    final static class LocusComparator {
        private final SAMSequenceDictionary dict;

        public LocusComparator( final SAMSequenceDictionary dict ) {
            this.dict = dict;
        }

        public int compareLocus( final String contig,
                                 final int position,
                                 final Locatable loc ) {
            int cmp = Integer.compare(dict.getSequenceIndex(contig),
                    dict.getSequenceIndex(loc.getContig()));
            if ( cmp == 0 ) {
                if ( position < loc.getStart() ) {
                    cmp = -1;
                } else if ( position > loc.getEnd() ) {
                    cmp = 1;
                }
            }
            return cmp;
        }
    }

    @VisibleForTesting
    final static class SiteDepthCounter {
        private final LocusComparator lComp;
        private final String sampleName;
        private final FeatureSink<SiteDepth> writer;
        private final int minMapQ;
        private final int minQ;
        private final Iterator<VariantContext> snpSourceItr;
        private final Deque<SiteDepth> siteDepthQueue;

        public SiteDepthCounter( final SAMSequenceDictionary dict,
                                 final String sampleName,
                                 final int compressionLevel,
                                 final GATKPath inputPath,
                                 final GATKPath outputPath,
                                 final int minMapQ,
                                 final int minQ ) {
            this.lComp = new LocusComparator(dict);
            this.sampleName = sampleName;
            final String outputFilename = outputPath.toPath().toString();
            final SiteDepthBCICodec bciCodec = new SiteDepthBCICodec();
            final List<String> sampleNames = Collections.singletonList(sampleName);
            if ( bciCodec.canDecode(outputFilename) ) {
                this.writer = bciCodec.makeSink(outputPath, dict, sampleNames, compressionLevel);
            } else {
                final SiteDepthCodec codec = new SiteDepthCodec();
                if ( !codec.canDecode(outputFilename) ) {
                    throw new UserException("Attempting to write site depth evidence to a file that " +
                            "can't be read as site depth evidence: " + outputFilename + ".  The file " +
                            "name should end with \".sd.txt\", \".sd.txt.gz\", or \".sd.bci\".");
                }
                this.writer = codec.makeSink(outputPath, dict, sampleNames, compressionLevel);
            }
            this.minMapQ = minMapQ;
            this.minQ = minQ;
            final FeatureDataSource<VariantContext> snpSource =
                    new FeatureDataSource<>(inputPath.toPath().toString(),
                                        null,
                                            FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                                            VariantContext.class);
            dict.assertSameDictionary(snpSource.getSequenceDictionary());
            this.snpSourceItr = new BAFSiteIterator(snpSource.iterator());
            this.siteDepthQueue = new ArrayDeque<>(100);
            readNextLocus();
        }

        public void apply( final GATKRead read ) {
            if ( read.getMappingQuality() < minMapQ || siteDepthQueue.isEmpty() ) {
                return;
            }

            // clean queue of SiteDepths that precede the current read
            final SimpleInterval readLoc =
                    new SimpleInterval(read.getContig(), read.getStart(), read.getEnd());
            while ( true ) {
                final SiteDepth siteDepth = siteDepthQueue.getFirst();
                if ( lComp.compareLocus(siteDepth.getContig(), siteDepth.getStart(), readLoc) >= 0 ) {
                    break;
                }
                writer.write(siteDepthQueue.removeFirst());
                if ( siteDepthQueue.isEmpty() ) {
                    if ( !readNextLocus() ) {
                        return;
                    }
                }
            }

            // make sure that the last SiteDepth in the queue occurs after the current read
            //  if such a SiteDepth is available
            while ( true ) {
                final SiteDepth siteDepth = siteDepthQueue.getLast();
                if ( lComp.compareLocus(siteDepth.getContig(), siteDepth.getStart(), readLoc) > 0 ||
                        !readNextLocus() ) {
                    break;
                }
            }

            walkReadMatches(read);
        }

        private void walkReadMatches( final GATKRead read ) {
            walkReadMatches(read, minQ, siteDepthQueue, lComp);
        }

        static void walkReadMatches( final GATKRead read,
                                     final int minQ,
                                     final Iterable<SiteDepth> sites,
                                     final LocusComparator locusComparator ) {
            int opStart = read.getStart();
            int readIdx = 0;
            final byte[] calls = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            for ( final CigarElement cigEle : read.getCigar().getCigarElements() ) {
                final int eleLen = cigEle.getLength();
                final CigarOperator cigOp = cigEle.getOperator();
                if ( cigOp.isAlignment() ) {
                    final int opEnd = opStart + eleLen - 1;
                    final SimpleInterval opLoc =
                            new SimpleInterval(read.getContig(), opStart, opEnd);
                    for ( final SiteDepth siteDepth : sites ) {
                        final String siteContig = siteDepth.getContig();
                        final int sitePos = siteDepth.getStart();
                        final int cmp = locusComparator.compareLocus(siteContig, sitePos, opLoc);
                        if ( cmp > 0 ) {
                            break;
                        }
                        // don't count base calls that aren't really part of the template
                        // (if the template is shorter than the read, we can call into adaptor sequence)
                        if ( cmp == 0 && !isBaseInsideAdaptor(read, sitePos) ) {
                            final int callIdx = readIdx + sitePos - opStart;
                            if ( quals[callIdx] < minQ ) {
                                continue;
                            }
                            final Nucleotide call = Nucleotide.decode(calls[callIdx]);
                            if ( call.isStandard() ) {
                                siteDepth.observe(call.ordinal());
                            }
                        }
                    }
                }
                if ( cigOp.consumesReadBases() ) {
                    readIdx += eleLen;
                }
                if ( cigOp.consumesReferenceBases() ) {
                    opStart += eleLen;
                }
            }
        }

        public void close() {
            while ( !siteDepthQueue.isEmpty() ) {
                writer.write(siteDepthQueue.removeFirst());
            }
            writer.close();
        }

        private boolean readNextLocus() {
            if ( !snpSourceItr.hasNext() ) {
                return false;
            }
            final SiteDepth siteDepth = new SiteDepth(snpSourceItr.next(), sampleName);
            siteDepthQueue.add(siteDepth);
            return true;
        }
    }

    public final static class BAFSiteIterator implements Iterator<VariantContext> {
        final Iterator<VariantContext> vcIterator;
        VariantContext last;
        VariantContext next;

        public BAFSiteIterator( final Iterator<VariantContext> vcIterator ) {
            this.vcIterator = vcIterator;
            advance();
        }

        public boolean hasNext() {
            if ( next != null ) {
                return true;
            }
            advance();
            return next != null;
        }

        public VariantContext next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException("baf sites iterator is exhausted");
            }
            final VariantContext result = next;
            last = next;
            next = null;
            return result;
        }

        private void advance() {
            while ( vcIterator.hasNext() ) {
                final VariantContext vc = vcIterator.next();
                // if it's a SNP, it's biallelic, and it occurs at a new locus
                if ( vc.isSNP() && vc.isBiallelic() &&
                        (last == null || !last.getContig().equals(vc.getContig()) ||
                                last.getStart() < vc.getStart()) ) {
                    next = vc;
                    break;
                }
            }
        }
    }

    @VisibleForTesting
    final static class CountCounter {
        private final int[] lowCounts;
        private final SortedMap<Integer, Integer> highCounts;
        private long nCounts;
        private long totalCounts;

        public CountCounter() {
            lowCounts = new int[3000];
            highCounts = new TreeMap<>();
        }

        public void addCount( final int count ) {
            nCounts += 1;
            totalCounts += count;
            if ( count < lowCounts.length ) {
                lowCounts[count] += 1;
            } else {
                highCounts.compute(count, (k,v) -> v==null ? 1 : v+1);
            }
        }

        public int getNZeroCounts() {
            return lowCounts[0];
        }

        public double getMeanCount() {
            return (double)totalCounts/nCounts;
        }

        public int[] getQuartiles() {
            final int[] quartiles = new int[4];
            int counts = 0;
            int quartile = 0;
            long targetCounts = (nCounts + 3) / 4;
            for ( int idx = 0; idx != lowCounts.length; ++idx ) {
                counts += lowCounts[idx];
                while ( counts >= targetCounts ) {
                    quartiles[quartile++] = idx;
                    targetCounts = ((quartile + 1) * nCounts + 3) / 4;
                }
            }
            for ( final Map.Entry<Integer,Integer> entry : highCounts.entrySet() ) {
                counts += entry.getValue();
                while ( counts >= targetCounts ) {
                    quartiles[quartile++] = entry.getKey();
                    targetCounts = ((quartile + 1) * nCounts + 3) / 4;
                }
            }
            return quartiles;
        }
    }

    @VisibleForTesting
    final static class DepthEvidenceCollector {
        private final LocusComparator lComp;
        private final CountCounter countCounter;
        private final FeatureSink<DepthEvidence> writer;
        private final Iterator<Feature> intervalIterator;
        private final int minMapQ;
        private DepthEvidence depthEvidence;
        private long summmedIntervalLengths;
        private long nIntervals;

        public DepthEvidenceCollector( final SAMSequenceDictionary dict,
                                       final String sampleName,
                                       final int cmprLevel,
                                       final GATKPath inputIntervalsPath,
                                       final GATKPath outputDepthEvidencePath,
                                       final int minMapQ ) {
            lComp = new LocusComparator(dict);
            countCounter = new CountCounter();
            final String outputFilename = outputDepthEvidencePath.toPath().toString();
            final DepthEvidenceBCICodec bciCodec = new DepthEvidenceBCICodec();
            final List<String> sampleNames = Collections.singletonList(sampleName);
            if ( bciCodec.canDecode(outputFilename) ) {
                writer = bciCodec.makeSink(outputDepthEvidencePath, dict, sampleNames, cmprLevel);
            } else {
                final DepthEvidenceCodec codec = new DepthEvidenceCodec();
                if ( !codec.canDecode(outputFilename) ) {
                    throw new UserException("Attempting to write depth evidence to a file that " +
                            "can't be read as depth evidence: " + outputFilename + ".  The file " +
                            "name should end with \".rd.txt\", \".rd.txt.gz\", or \".rd.bci\".");
                }
                writer = codec.makeSink(outputDepthEvidencePath, dict, sampleNames, cmprLevel);
            }

            final FeatureDataSource<Feature> intervalSource =
                    new FeatureDataSource<>(inputIntervalsPath.toPath().toString());
            this.intervalIterator = intervalSource.iterator();
            if ( !intervalIterator.hasNext() ) {
                throw new UserException(inputIntervalsPath + " contains no intervals.");
            }
            this.minMapQ = minMapQ;
            depthEvidence = new DepthEvidence(intervalIterator.next(), new int[1]);
            summmedIntervalLengths += depthEvidence.getLengthOnReference();
            nIntervals += 1;
        }

        void apply( final GATKRead read ) {
            if ( read.getMappingQuality() < minMapQ ) {
                return;
            }
            while ( depthEvidence != null ) {
                final int cmp = lComp.compareLocus(read.getContig(), read.getStart(), depthEvidence);
                if ( cmp < 0 ) { // if read is upstream of interval of interest, nothing to do
                    return;
                }
                if ( cmp == 0 ) { // if read is in the interval of interest, count it
                    depthEvidence.getCounts()[0] += 1;
                    return;
                }

                // read is downstream of interval
                writer.write(depthEvidence);
                countCounter.addCount(depthEvidence.getCounts()[0]);
                if ( !intervalIterator.hasNext() ) {
                    depthEvidence = null;
                    return;
                }
                depthEvidence = new DepthEvidence(intervalIterator.next(), new int[1]);
                summmedIntervalLengths += depthEvidence.getLengthOnReference();
                nIntervals += 1;
            }
        }

        void close() {
            if ( depthEvidence != null ) {
                writer.write(depthEvidence);
                countCounter.addCount(depthEvidence.getCounts()[0]);
            }
            writer.close();
        }

        void reportSummaryStats( final GATKPath summaryPath, final String sampleName ) {
            try ( final BufferedWriter writer
                     = new BufferedWriter(new OutputStreamWriter(summaryPath.getOutputStream())) ) {
                final int[] quartiles = countCounter.getQuartiles();
                writer.write("rd_q25\t" + sampleName + "\t" + quartiles[0]);
                writer.newLine();
                writer.write("rd_q50\t" + sampleName + "\t" + quartiles[1]);
                writer.newLine();
                writer.write("rd_q75\t" + sampleName + "\t" + quartiles[2]);
                writer.newLine();
                final String val = String.format("%.2f",countCounter.getMeanCount());
                writer.write("rd_mean\t" + sampleName + "\t" + val);
                writer.newLine();
                writer.write("rd_num_intervals\t" + sampleName + "\t" + nIntervals);
                writer.newLine();
                writer.write("rd_intervals_size\t" + sampleName + "\t" + summmedIntervalLengths);
                writer.newLine();
            } catch ( final IOException ioe ) {
                throw new UserException("Can't write depth evidence summary statistics to " + summaryPath);
            }
        }
    }
}
