package org.broadinstitute.hellbender.tools.walkers.indels;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * Define intervals to target for local realignment
 * <p>
 * <p>
 * The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases
 * is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion
 * or deletion (indels) in the individual's genome with respect to the reference genome.  Such alignment artifacts result in many bases mismatching
 * the reference near the misalignment, which are easily mistaken as SNPs.  Moreover, since read mapping algorithms operate on each read independently,
 * it is impossible to place reads on the reference genome such that mismatches are minimized across all reads.  Consequently, even when some reads are
 * correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel,
 * also requiring realignment.  Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus
 * indel suitable for standard variant discovery approaches.
 * </p>
 * <p><b>Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly
 * step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect.</b></p>
 * <p>There are 2 steps to the realignment process:</p>
 * <ol>
 * <li>Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)</li>
 * <li>Running the realigner over those intervals (see the IndelRealigner tool)</li>
 * </ol>
 * <p>
 * <p>
 * For more details, see <a href="http://www.broadinstitute.org/gatk/guide/article?id=38">the indel realignment method documentation</a>.
 * </p>
 * <p>
 * <h3>Inputs</h3>
 * <p>
 * One or more aligned BAM files and optionally, one or more lists of known indels.
 * </p>
 * <p>
 * <h3>Output</h3>
 * <p>
 * A list of target intervals to pass to the IndelRealigner.
 * </p>
 * <p>
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T RealignerTargetCreator \
 *   -R reference.fasta \
 *   -I input.bam \
 *   --known indels.vcf \
 *   -o forIndelRealigner.intervals
 * </pre>
 * <p>
 * <h3>Notes</h3>
 * <ul>
 * <li>The input BAM(s), reference, and known indel file(s) should be the same ones to be used for the IndelRealigner step.</li>
 * <li>When multiple potential indels are found by the tool in the same general region, the tool will choose the most likely
 * one for realignment to the exclusion of the others.  This is a known limitation of the tool.</li>
 * <li>Because reads produced from the 454 technology inherently contain false indels, the realigner will not work with them
 * (or with reads from similar technologies).</li>
 * <li>This tool also ignores MQ0 reads and reads with consecutive indel operators in the CIGAR string.</li>
 * </ul>
 */
@CommandLineProgramProperties(oneLineSummary = "Define intervals to target for local realignment",
        summary = "The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases " +
                "is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion " +
                "or deletion (indels) in the individual's genome with respect to the reference genome.  Such alignment artifacts result in many bases mismatching " +
                "the reference near the misalignment, which are easily mistaken as SNPs.  Moreover, since read mapping algorithms operate on each read independently, " +
                "it is impossible to place reads on the reference genome such at mismatches are minimized across all reads.  Consequently, even when some reads are " +
                "correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, " +
                "also requiring realignment.  Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus " +
                "indel suitable for standard variant discovery approaches.",
        programGroup = ReadProgramGroup.class

)
public class RealignerTargetCreator extends LocusWalker {

    /**
     * The target intervals for realignment.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    protected File out;

    /**
     * Any number of VCF files representing known SNPs and/or indels.  Could be e.g. dbSNP and/or official 1000 Genomes indel calls.
     * SNPs in these files will be ignored unless the --mismatchFraction argument is used.
     */
    @Argument(fullName = "known", shortName = "known", doc = "Input VCF file with known indels", optional = true)
    public List<FeatureInput<VariantContext>> known = new ArrayList<>();

    /**
     * Any two SNP calls and/or high entropy positions are considered clustered when they occur no more than this many basepairs apart. Must be > 1.
     */
    @Argument(fullName = "windowSize", shortName = "window", doc = "window size for calculating entropy or SNP clusters", optional = true)
    protected int windowSize = 10;

    /**
     * To disable this behavior, set this value to <= 0 or > 1.  This feature is really only necessary when using an ungapped aligner
     * (e.g. MAQ in the case of single-end read data) and should be used in conjunction with '--model USE_SW' in the IndelRealigner.
     */
    @Argument(fullName = "mismatchFraction", shortName = "mismatch", doc = "fraction of base qualities needing to mismatch for a position to have high entropy", optional = true)
    protected double mismatchThreshold = 0.0;

    @Argument(fullName = "minReadsAtLocus", shortName = "minReads", doc = "minimum reads at a locus to enable using the entropy calculation", optional = true)
    protected int minReadsAtLocus = 4;

    /**
     * Because the realignment algorithm is N^2, allowing too large an interval might take too long to completely realign.
     */
    @Argument(fullName = "maxIntervalSize", shortName = "maxInterval", doc = "maximum interval size; any intervals larger than this value will be dropped", optional = true)
    protected int maxIntervalSize = 500;

    // genome location parser and sam file header
    private GenomeLocParser genomeLocParser;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean includeDeletions() {
        return true;
    }

    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(9);
        // read filters from Walker (original implementation)
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(ReadFilterLibrary.GOOD_CIGAR);
        // read filters from LocusWalker (original implementation)
        defaultFilters.add(ReadFilterLibrary.MAPPED);
        defaultFilters.add(ReadFilterLibrary.PRIMARY_ALIGNMENT);
        defaultFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        defaultFilters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        // read filters from RealignerTargetCreator (original implementation)
        defaultFilters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        defaultFilters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        defaultFilters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        // TODO: original had Platform454Filter

        return defaultFilters;
    }

    private boolean lookForMismatchEntropy;

    @Override
    public String[] customCommandLineValidation() {
        if (windowSize < 2)
            throw new CommandLineException.BadArgumentValue("windowSize", "Window Size must be an integer greater than 1");

        lookForMismatchEntropy = mismatchThreshold > 0.0 && mismatchThreshold <= 1.0;
        return super.customCommandLineValidation();
    }

    @Override
    public void onTraversalStart() {
        genomeLocParser = new GenomeLocParser(getBestAvailableSequenceDictionary());
    }

    @Override
    public void apply(AlignmentContext context, ReferenceContext ref, FeatureContext tracker) {
        // TODO: this comes from the annotation @Reference(window=@Window(start=-1,stop=50))
        ref.setWindow(0, 50);

        boolean hasIndel = false;
        boolean hasInsertion = false;
        boolean hasPointEvent = false;

        int furthestStopPos = -1;

        // look at the rods for indels or SNPs
        if (!known.isEmpty()) {
            for (VariantContext vc : tracker.getValues(known)) {
                switch (vc.getType()) {
                    case INDEL:
                        hasIndel = true;
                        if (vc.isSimpleInsertion())
                            hasInsertion = true;
                        break;
                    case SNP:
                        hasPointEvent = true;
                        break;
                    case MIXED:
                        hasPointEvent = true;
                        hasIndel = true;
                        if (vc.isSimpleInsertion())
                            hasInsertion = true;
                        break;
                    default:
                        break;
                }
                if (hasIndel)
                    furthestStopPos = vc.getEnd();
            }
        }

        // look at the normal context to get deletions and positions with high entropy
        final ReadPileup pileup = context.getBasePileup();

        int mismatchQualities = 0, totalQualities = 0;
        final byte refBase = ref.getBase();
        for (PileupElement p : pileup) {

            // check the ends of the reads to see how far they extend
            furthestStopPos = Math.max(furthestStopPos, p.getRead().getEnd());

            // is it a deletion or insertion?
            if (p.isDeletion() || p.isBeforeInsertion()) {
                hasIndel = true;
                if (p.isBeforeInsertion())
                    hasInsertion = true;
            }

            // look for mismatches
            else if (lookForMismatchEntropy) {
                if (p.getBase() != refBase)
                    mismatchQualities += p.getQual();
                totalQualities += p.getQual();
            }
        }

        // make sure we're supposed to look for high entropy
        if (lookForMismatchEntropy &&
                pileup.size() >= minReadsAtLocus &&
                (double) mismatchQualities / (double) totalQualities >= mismatchThreshold)
            hasPointEvent = true;

        // return null if no event occurred
        if (!hasIndel && !hasPointEvent) {
            return;
        }
        // return null if we didn't find any usable reads/rods associated with the event
        if (furthestStopPos == -1) {
            return;
        }

        GenomeLoc eventLoc = genomeLocParser.createGenomeLoc(context.getLocation());
        if (hasInsertion)
            eventLoc = genomeLocParser.createGenomeLoc(eventLoc.getContig(), eventLoc.getStart(), eventLoc.getStart() + 1);

        EVENT_TYPE eventType = (hasIndel ? (hasPointEvent ? EVENT_TYPE.BOTH : EVENT_TYPE.INDEL_EVENT) : EVENT_TYPE.POINT_EVENT);

        reduce(new Event(eventLoc, furthestStopPos, eventType), sum);
    }

    // accumulator
    private final EventPair sum = new EventPair(null, null);

    @Override
    public Object onTraversalSuccess() {
        if (sum.left != null && sum.left.isReportableEvent())
            sum.intervals.add(sum.left.getLoc());
        if (sum.right != null && sum.right.isReportableEvent())
            sum.intervals.add(sum.right.getLoc());

        if (FilenameUtils.getExtension(out.getName()).equals("interval_list")) {
            final SAMFileHeader masterSequenceDictionaryHeader = new SAMFileHeader();
            masterSequenceDictionaryHeader.setSequenceDictionary(getBestAvailableSequenceDictionary());
            final IntervalList intervalList = new IntervalList(masterSequenceDictionaryHeader);
            for (GenomeLoc loc : sum.intervals) {
                intervalList.add(new Interval(loc.getContig(), loc.getStart(), loc.getStop()));
            }
            intervalList.write(out);
        } else {
            try (BufferedWriter bufferedWriter = IOUtil.openFileForBufferedWriting(out)) {
                for (GenomeLoc loc : sum.intervals) {
                    bufferedWriter.write(loc.toString());
                    bufferedWriter.newLine();
                }
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(out, "error writing out intervals to file", e);
            }
        }
        return null;
    }

    public EventPair reduce(Event value, EventPair sum) {
        if (value == null) {
            ; // do nothing
        } else if (sum.left == null) {
            sum.left = value;
        } else if (sum.right == null) {
            if (canBeMerged(sum.left, value))
                sum.left = mergeEvents(sum.left, value);
            else
                sum.right = value;
        } else {
            if (canBeMerged(sum.right, value))
                sum.right = mergeEvents(sum.right, value);
            else {
                if (sum.right.isReportableEvent())
                    sum.intervals.add(sum.right.getLoc());
                sum.right = value;
            }
        }

        return sum;
    }

    static private boolean canBeMerged(Event left, Event right) {
        return left.loc.getContig().equals(right.loc.getContig()) && left.furthestStopPos >= right.loc.getStart();
    }

    static private Event mergeEvents(Event left, Event right) {
        Utils.nonNull(left);
        Utils.nonNull(right);
        left.merge(right);
        return left;
    }

    private enum EVENT_TYPE {POINT_EVENT, INDEL_EVENT, BOTH}

    static class EventPair {
        public Event left, right;
        public TreeSet<GenomeLoc> intervals = new TreeSet<GenomeLoc>();

        public EventPair(Event left, Event right) {
            this.left = left;
            this.right = right;
        }

        public EventPair(Event left, Event right, TreeSet<GenomeLoc> set1, TreeSet<GenomeLoc> set2) {
            this.left = left;
            this.right = right;
            intervals.addAll(set1);
            intervals.addAll(set2);
        }
    }

    class Event {
        public int furthestStopPos;

        private GenomeLoc loc;
        private int eventStartPos;
        private int eventStopPos;
        private EVENT_TYPE type;
        private ArrayList<Integer> pointEvents = new ArrayList<Integer>();

        public Event(GenomeLoc loc, int furthestStopPos, EVENT_TYPE type) {
            this.loc = loc;
            this.furthestStopPos = furthestStopPos;
            this.type = type;

            if (type == EVENT_TYPE.INDEL_EVENT || type == EVENT_TYPE.BOTH) {
                eventStartPos = loc.getStart();
                eventStopPos = loc.getStop();
            } else {
                eventStartPos = -1;
                eventStopPos = -1;
            }

            if (type == EVENT_TYPE.POINT_EVENT || type == EVENT_TYPE.BOTH) {
                pointEvents.add(loc.getStart());
            }
        }

        public void merge(Event e) {

            // merges only get called for events with certain types
            if (e.type == EVENT_TYPE.INDEL_EVENT || e.type == EVENT_TYPE.BOTH) {
                if (eventStartPos == -1)
                    eventStartPos = e.eventStartPos;
                eventStopPos = e.eventStopPos;
                furthestStopPos = e.furthestStopPos;
            }

            if (e.type == EVENT_TYPE.POINT_EVENT || e.type == EVENT_TYPE.BOTH) {
                int newPosition = e.pointEvents.get(0);
                if (pointEvents.size() > 0) {
                    int lastPosition = pointEvents.get(pointEvents.size() - 1);
                    if (newPosition - lastPosition < windowSize) {
                        eventStopPos = Math.max(eventStopPos, newPosition);
                        furthestStopPos = e.furthestStopPos;

                        if (eventStartPos == -1)
                            eventStartPos = lastPosition;
                        else
                            eventStartPos = Math.min(eventStartPos, lastPosition);
                    } else if (eventStartPos == -1 && e.eventStartPos != -1) {
                        eventStartPos = e.eventStartPos;
                        eventStopPos = e.eventStopPos;
                        furthestStopPos = e.furthestStopPos;
                    }
                }
                pointEvents.add(newPosition);
            }
        }

        public boolean isReportableEvent() {
            return genomeLocParser.isValidGenomeLoc(loc.getContig(), eventStartPos, eventStopPos, true) && eventStopPos >= 0 && eventStopPos - eventStartPos < maxIntervalSize;
        }

        public GenomeLoc getLoc() {
            return genomeLocParser.createGenomeLoc(loc.getContig(), eventStartPos, eventStopPos);
        }
    }
}