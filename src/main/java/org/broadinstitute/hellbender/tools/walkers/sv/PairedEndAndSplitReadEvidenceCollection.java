package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.*;
import java.util.function.Predicate;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.isBaseInsideAdaptor;

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
    public static final String ALLELE_COUNT_FILE_ARGUMENT_SHORT_NAME = "AC";
    public static final String ALLELE_COUNT_FILE_ARGUMENT_LONG_NAME = "allele-count-file";
    public static final String ALLELE_COUNT_VCF_ARGUMENT_SHORT_NAME = "F";
    public static final String ALLELE_COUNT_VCF_ARGUMENT_LONG_NAME = "allele-count-vcf";
    public static final String SAMPLE_NAME_ARGUMENT_LONG_NAME = "sample-name";

    @Argument(shortName = PAIRED_END_FILE_ARGUMENT_SHORT_NAME, fullName = PAIRED_END_FILE_ARGUMENT_LONG_NAME, doc = "Output file for paired end evidence", optional=false)
    public String peFile;

    @Argument(shortName = SPLIT_READ_FILE_ARGUMENT_SHORT_NAME, fullName = SPLIT_READ_FILE_ARGUMENT_LONG_NAME, doc = "Output file for split read evidence", optional=false)
    public String srFile;

    @Argument(shortName = ALLELE_COUNT_FILE_ARGUMENT_SHORT_NAME,
            fullName = ALLELE_COUNT_FILE_ARGUMENT_LONG_NAME,
            doc = "Output file for allele counts",
            optional = true)
    public String alleleCountFilename;

    @Argument(shortName = ALLELE_COUNT_VCF_ARGUMENT_SHORT_NAME,
            fullName = ALLELE_COUNT_VCF_ARGUMENT_LONG_NAME,
            doc = "Input VCF of SNPs marking loci for allele counts",
            optional = true)
    public String vcfFilename;

    @Argument(fullName = "sv-pipe-counts",
            doc = "True for SV pipeline style allele-count output (designated alt count reported), " +
                    "otherwise CollectAllelicCounts style (max alt attracts all non-ref counts)",
            optional = true)
    public boolean svPipeStyleAlleleCounts = false;

    @Argument(fullName = "allele-count-min-mapq",
            doc = "minimum mapping quality for read to be allele-counted",
            optional = true)
    public int minMapQ = 30;

    @Argument(fullName = "allele-count-min-baseq",
            doc = "minimum base call quality for SNP to be allele-counted",
            optional = true)
    public int minQ = 20;

    @Argument(fullName = SAMPLE_NAME_ARGUMENT_LONG_NAME, doc = "Sample name")
    String sampleName = null;

    final Set<String> observedDiscordantNames = new HashSet<>();
    final PriorityQueue<SplitPos> splitPosBuffer = new PriorityQueue<>(new SplitPosComparator());
    final List<DiscordantRead> discordantPairs = new ArrayList<>();

    int currentDiscordantPosition = -1;
    String currentChrom = null;

    private BufferedWriter peWriter;
    private BufferedWriter srWriter;
    private AlleleCounter alleleCounter;

    private SAMSequenceDictionary sequenceDictionary;

    @Override
    public boolean requiresReads() {
        return true;
    }


    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        peWriter = createOutputFile(peFile);
        srWriter = createOutputFile(srFile);

        sequenceDictionary = getBestAvailableSequenceDictionary();
        if ( vcfFilename != null && alleleCountFilename != null ) {
            final LocusCountsWriter writer = svPipeStyleAlleleCounts ?
                    new LocusCountsVCFWriter(sequenceDictionary, sampleName,
                            createVCFWriter(new GATKPath(alleleCountFilename))) :
                    new LocusCountsFileWriter(sequenceDictionary, alleleCountFilename);
            alleleCounter = new AlleleCounter(sequenceDictionary, vcfFilename, minMapQ, minQ, writer);
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
        if ( read.isPaired() &&
                !read.mateIsUnmapped() && !read.isProperlyPaired() &&
                !read.isSupplementaryAlignment() && !read.isSecondaryAlignment() ) {
            if ( isSoftClipped(read) ) {
                countSplitRead(read, splitPosBuffer, srWriter);
            }

            reportDiscordantReadPair(read);
        }

        if ( alleleCounter != null ) {
            alleleCounter.apply(read);
        }
    }

    private static BufferedWriter createOutputFile( final String fileName ) {
        final GATKPath path = new GATKPath(fileName);
        if ( fileName.endsWith(".gz") ) {
            return new BufferedWriter(
                    new OutputStreamWriter(
                            new BlockCompressedOutputStream(path.getOutputStream(), path.toPath(), 6)));
        }
        return new BufferedWriter(new OutputStreamWriter(path.getOutputStream()));
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
        final String strandA = r.isReadReverseStrand() ? "-" : "+";
        final String strandB = r.isMateReverseStrand() ? "-" : "+";

        try {
            // subtract 1 from positions to match pysam output
            peWriter.write(r.getContig() + "\t" + (r.getStart() - 1) + "\t" + strandA + "\t" + r.getMateContig() + "\t" + (r.getMateStart() - 1) + "\t" + strandB + "\t" + sampleName + "\n");
        } catch (IOException e) {
            throw new GATKException("Could not write to PE file", e);
        }
    }

    /**
     * Adds split read information about the current read to the counts in splitCounts. Flushes split read counts to
     * srWriter if necessary.
     */
    @VisibleForTesting
    public void countSplitRead(final GATKRead read, final PriorityQueue<SplitPos> splitCounts, final BufferedWriter srWriter) {
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

    private void flushSplitCounts(final Predicate<SplitPos> flushablePosition, final BufferedWriter srWriter, final PriorityQueue<SplitPos> splitCounts) {

        while (splitCounts.size() > 0 && flushablePosition.test(splitCounts.peek())) {
            SplitPos pos = splitCounts.poll();
            int countAtPos = 1;
            while (splitCounts.size() > 0 && splitCounts.peek().equals(pos)) {
                countAtPos++;
                splitCounts.poll();
            }
            try {
                srWriter.write(currentChrom + "\t" + (pos.pos - 1) + "\t" + pos.direction.getDescription() + "\t" + countAtPos + "\t" + sampleName + "\n");
            } catch (IOException e) {
                throw new GATKException("Could not write to sr file", e);
            }
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
        if ( alleleCounter != null ) {
            alleleCounter.close();
        }
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        try {
            peWriter.close();
            srWriter.close();
        } catch (IOException e) {
            throw new GATKException("error closing output file", e);
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

    @VisibleForTesting
    final static class AlleleCounter {
        private final SAMSequenceDictionary dict;
        private final int minMapQ;
        private final int minQ;
        private final Iterator<VariantContext> snpSourceItr;
        private final LocusCountsWriter writer;
        private final Deque<LocusCounts> locusCountsQueue;

        public AlleleCounter( final SAMSequenceDictionary dict,
                              final String vcfFilename,
                              final int minMapQ,
                              final int minQ,
                              final LocusCountsWriter writer ) {
            this.dict = dict;
            this.minMapQ = minMapQ;
            this.minQ = minQ;
            final FeatureDataSource<VariantContext> snpSource = new FeatureDataSource<>(vcfFilename);
            dict.assertSameDictionary(snpSource.getSequenceDictionary());
            this.snpSourceItr = snpSource.iterator();
            this.writer = writer;
            this.locusCountsQueue = new ArrayDeque<>(100);
            readLocus();
        }

        public void apply( final GATKRead read ) {
            if ( read.getMappingQuality() < minMapQ || locusCountsQueue.isEmpty() ) {
                return;
            }

            // clean queue of LocusCounts that precede the current read
            final int contigId = dict.getSequenceIndex(read.getContig());
            final Locus readStart = new Locus(contigId, read.getStart());
            while ( locusCountsQueue.getFirst().compareTo(readStart) < 0 ) {
                writer.write(locusCountsQueue.removeFirst());
                if ( locusCountsQueue.isEmpty() ) {
                    if ( !readLocus() ) {
                        return;
                    }
                }
            }

            // make sure that the last LocusCount in the queue happens after the current read
            //  if such a LocusCount is available
            final Locus readEnd = new Locus(contigId, read.getEnd());
            while ( readEnd.compareTo(locusCountsQueue.getLast()) > 0 ) {
                if ( !readLocus() ) {
                    break;
                }
            }

            walkReadMatches(read, contigId);
        }

        public void walkReadMatches( final GATKRead read, final int contigId ) {
            Locus opStart = new Locus(contigId, read.getStart());
            int readIdx = 0;
            final byte[] calls = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            for ( final CigarElement cigEle : read.getCigar().getCigarElements() ) {
                final int eleLen = cigEle.getLength();
                final CigarOperator cigOp = cigEle.getOperator();
                if ( cigOp.isAlignment() ) {
                    final Locus opEnd =
                            new Locus(opStart.getContigId(), opStart.getPosition() + eleLen);
                    for ( final LocusCounts locusCounts : locusCountsQueue ) {
                        if ( locusCounts.compareTo(opEnd) >= 0 ) {
                            break;
                        }
                        if ( opStart.compareTo(locusCounts) <= 0 &&
                                !isBaseInsideAdaptor(read, locusCounts.getPosition() )) {
                            final int callIdx =
                                    readIdx + locusCounts.getPosition() - opStart.getPosition();
                            if ( quals[callIdx] < minQ ) {
                                continue;
                            }
                            final Nucleotide call = Nucleotide.decode(calls[callIdx]);
                            if ( call.isStandard() ) {
                                locusCounts.observe(call.ordinal());
                            }
                        }
                    }
                }
                if ( cigOp.consumesReadBases() ) {
                    readIdx += cigEle.getLength();
                }
                if ( cigOp.consumesReferenceBases() ) {
                    opStart = new Locus(opStart.getContigId(), opStart.getPosition() + eleLen);
                }
            }
        }

        public void close() {
            while ( !locusCountsQueue.isEmpty() ) {
                writer.write(locusCountsQueue.removeFirst());
            }
            writer.close();
        }

        private boolean readLocus() {
            if ( !snpSourceItr.hasNext() ) {
                return false;
            }
            VariantContext snp = snpSourceItr.next();
            while ( !snp.isSNP() ) {
                if ( !snpSourceItr.hasNext() ) {
                    return false;
                }
                snp = snpSourceItr.next();
            }
            final int contigId = dict.getSequenceIndex(snp.getContig());
            final int position = snp.getStart();
            final byte[] refSeq = snp.getReference().getBases();
            final Nucleotide refCall = Nucleotide.decode(refSeq[0]);
            if ( !refCall.isStandard() ) {
                throw new UserException("vcf contains a SNP with a non-standard reference base " +
                        refCall + " at locus " + snp.getContig() + ":" + position);
            }
            final byte[] altSeq = snp.getAlternateAllele(0).getBases();
            final Nucleotide altCall = Nucleotide.decode(altSeq[0]);
            if ( !altCall.isStandard() ) {
                throw new UserException("vcf contains a SNP with a non-standard alt base" +
                        altCall + " at locus " + snp.getContig() + ":" + position);
            }
            final LocusCounts locusCounts =
                    new LocusCounts(contigId, position, refCall.ordinal(), altCall.ordinal());
            locusCountsQueue.add(locusCounts);
            return true;
        }
    }

    public interface LocusCountsWriter {
        void write( LocusCounts locusCounts );
        void close();
    }

    @VisibleForTesting
    static class LocusCountsFileWriter implements LocusCountsWriter {
        private final SAMSequenceDictionary dict;
        private final String outputFilename;
        private final BufferedWriter writer;

        public LocusCountsFileWriter( final SAMSequenceDictionary dict,
                                      final String outputFilename ) {
            this.dict = dict;
            this.outputFilename = outputFilename;
            this.writer = createOutputFile(outputFilename);
        }

        public void write( final LocusCounts locusCounts ) {
            try {
                writer.write(dict.getSequence(locusCounts.getContigId()).getContig());
                writer.write('\t');
                writer.write(Integer.toString(locusCounts.getPosition()));
                writer.write('\t');
                final int refIdx = locusCounts.getRefIdx();
                final int refCount = locusCounts.getCount(refIdx);
                writer.write(Integer.toString(refCount));
                writer.write('\t');
                writer.write(Integer.toString(locusCounts.getTotalCounts() - refCount));
                writer.write('\t');
                writer.write("ACGT".charAt(refIdx));
                writer.write('\t');
                final int altIdx = locusCounts.getMaxAltIdx();
                writer.write(altIdx == LocusCounts.NO_ALT_IDX ? 'N' : "ACGT".charAt(altIdx));
                writer.newLine();
            } catch ( final IOException ioe ) {
                throw new UserException("can't write to allele count output file " + outputFilename, ioe);
            }
        }

        public void close() {
            try {
                writer.close();
            } catch ( final IOException ioe ) {
                throw new UserException("Can't close " + outputFilename);
            }
        }
    }

    @VisibleForTesting
    static class LocusCountsVCFWriter implements LocusCountsWriter {
        final SAMSequenceDictionary dict;
        final String sampleName;
        final VariantContextWriter writer;

        final static String DEPTH_KEY = "DP";
        final static String ALLELE_DEPTH_KEY = "AD";
        final static Allele[][][] allelePairs = {
                {       {},
                        {Allele.REF_A, Allele.ALT_C},
                        {Allele.REF_A, Allele.ALT_G},
                        {Allele.REF_A, Allele.ALT_T} },
                {       {Allele.REF_C, Allele.ALT_A},
                        {},
                        {Allele.REF_C, Allele.ALT_G},
                        {Allele.REF_C, Allele.ALT_T} },
                {       {Allele.REF_G, Allele.ALT_A},
                        {Allele.REF_G, Allele.ALT_C},
                        {},
                        {Allele.REF_G, Allele.ALT_T} },
                {       {Allele.REF_T, Allele.ALT_A},
                        {Allele.REF_T, Allele.ALT_C},
                        {Allele.REF_T, Allele.ALT_G},
                        {} } };

        public LocusCountsVCFWriter( final SAMSequenceDictionary dict,
                                     final String sampleName,
                                     final VariantContextWriter writer ) {
            this.dict = dict;
            this.sampleName = sampleName;
            this.writer = writer;
            final VCFHeader header = new VCFHeader(VCFHeaderVersion.VCF4_2, Collections.emptySet(),
                                                    Collections.singleton(sampleName));
            header.setSequenceDictionary(dict);
            header.addMetaDataLine(
                    new VCFFormatHeaderLine(DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Read depth"));
            header.addMetaDataLine(
                    new VCFFormatHeaderLine(ALLELE_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Alt allele depth"));
            writer.writeHeader(header);
        }

        public void write( final LocusCounts locusCounts ) {
            writer.add(new VariantContextBuilder(
                    "AlleleCounter",
                    dict.getSequence(locusCounts.getContigId()).getSequenceName(),
                    locusCounts.getPosition(),
                    locusCounts.getPosition(),
                    Arrays.asList(allelePairs[locusCounts.getRefIdx()][locusCounts.getAltIdx()]))
                    .genotypes(new GenotypeBuilder(sampleName)
                            .DP(locusCounts.getTotalCounts())
                            .AD(new int[]{locusCounts.getCount(locusCounts.getAltIdx())})
                            .make())
                    .make());
        }

        public void close() {
            writer.close();
        }
    }

    @VisibleForTesting
    static class Locus implements Comparable<Locus> {
        final int contigId;
        final int position;

        public Locus( final int contigId, final int position ) {
            this.contigId = contigId;
            this.position = position;
        }

        public int getContigId() { return contigId; }
        public int getPosition() { return position; }

        @Override public int compareTo( final Locus locus ) {
            int result = Integer.compare(contigId, locus.contigId);
            if ( result == 0 ) {
                result = Integer.compare(position, locus.position);
            }
            return result;
        }
    }

    @VisibleForTesting
    final static class LocusCounts extends Locus {
        public static final int NO_ALT_IDX = -1;
        private final int[] counts;
        private final int refIdx;
        private final int altIdx;

        public LocusCounts( final int contigId, final int position,
                            final int refIdx, final int altIdx ) {
            super(contigId, position);
            this.counts = new int[4];
            this.refIdx = refIdx;
            this.altIdx = altIdx;
        }

        public void observe( final int idx ) { counts[idx] += 1; }

        public int getRefIdx() { return refIdx; }
        public int getAltIdx() { return altIdx; }

        public int getMaxAltIdx() {
            int altIdx = NO_ALT_IDX;
            int altCount = 0;
            for ( int idx = 0; idx != 4; ++idx ) {
                if ( idx != refIdx ) {
                    final int count = counts[idx];
                    if ( count > altCount ) {
                        altIdx = idx;
                        altCount = count;
                    }
                }
            }
            return altIdx;
        }

        public int getCount( final int idx ) {
            return counts[idx];
        }
        public int getTotalCounts() { return counts[0] + counts[1] + counts[2] + counts[3]; }
    }
}
