package org.broadinstitute.hellbender.tools.walkers.chimericreads;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadPairWalker;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Examines Readpairs that have been annotated with AnnotateChimericReads, and alignes the embedded references to each other. " +
                "Can be used to figure out if references located near read pairs have similar sub-sequences. Annotates reads with alignment," +
                "and emits a metric with a summary of results. Expects queryname-sorted input, but emits coordinte-sorted file.",
        oneLineSummary = "Examines pairs of chimeric reads and aligns embedded references to each other.",
        programGroup = OtherProgramGroup.class
)
@ExperimentalFeature
public class ExamineChimericReads extends ReadPairWalker {
    private final SmithWatermanAligner smithWatermanAligner = SmithWatermanJavaAligner.getInstance();
    private final SWOverhangStrategy swOverhangStrategy = SWOverhangStrategy.SOFTCLIP;

    private SWParameters swParameters = null;

    @Argument(fullName = "match_value")
    int MATCH_VALUE = 5;
    @Argument(fullName = "mismatch_value")
    int MISMATCH_PENALTY_VALUE = -4;
    @Argument(fullName = "gap_open_value")
    int GAP_OPEN_VALUE = -12;
    @Argument(fullName = "gap_extend_value")
    int GAP_EXTEND_VALUE = -2;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Write output Reads to this file")
    public String outputBam;

    @Argument(fullName = StandardArgumentDefinitions.METRICS_FILE_LONG_NAME,
            shortName = StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME,
            doc = "Write output metrics to this file")
    public String outputMetrics;

    @Argument(fullName = "reference-bases-tag-name", shortName = "tn")
    String tagName = "rb";

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        swParameters = new SWParameters(MATCH_VALUE, MISMATCH_PENALTY_VALUE, GAP_OPEN_VALUE, GAP_EXTEND_VALUE);
        outputWriter = createSAMWriter(IOUtils.getPath(outputBam), false);
    }

    final private Histogram<Integer> cisMatchLength = new Histogram<>();
    final private Histogram<Integer> transMatchLength = new Histogram<>();
    private SAMFileGATKReadWriter outputWriter;

    @Override
    protected SAMFileHeader getHeaderForSAMWriter() {
        final SAMFileHeader headerForSAMWriter = super.getHeaderForSAMWriter();

        headerForSAMWriter.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return headerForSAMWriter;
    }

    @Override
    public void apply(final Set<GATKRead> reads) {

        final GATKRead readOne = reads.stream().filter(GATKRead::isFirstOfPair).findFirst().orElse(null);
        final GATKRead readTwo = reads.stream().filter(gatkRead -> !gatkRead.isFirstOfPair()).findFirst().orElse(null);

        if (readOne == null || readTwo == null ||
                !readOne.hasAttribute(tagName) ||
                !readTwo.hasAttribute(tagName)) {
            return;
        }

        final String refOne = readOne.getAttributeAsString(tagName);
        final String refTwo = readTwo.getAttributeAsString(tagName);

        final byte[] seq1 = refOne.getBytes();
        {
            final byte[] seq2 = refTwo.getBytes();
            final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);

            final int matches = countEquals(alignmentMis.getCigar());
            cisMatchLength.increment(matches);
            readOne.setAttribute("c1", alignmentMis.getCigar().toString());
            readTwo.setAttribute("c1", alignmentMis.getCigar().toString());
            readOne.setAttribute("m1", matches);
            readTwo.setAttribute("m1", matches);
        }
        {
            final byte[] seq2 = SequenceUtil.reverseComplement(refTwo).getBytes();
            final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);
            readOne.setAttribute("c2", alignmentMis.getCigar().toString());
            readTwo.setAttribute("c2", alignmentMis.getCigar().toString());

            final int revCompMatches = countEquals(alignmentMis.getCigar());
            readOne.setAttribute("m2", revCompMatches);
            readTwo.setAttribute("m2", revCompMatches);
            transMatchLength.increment(revCompMatches);
        }

        outputWriter.addRead(readOne);
        outputWriter.addRead(readTwo);
    }

    private int countEquals(final Cigar cigar) {
        return countMaxOps(cigar, CigarOperator.EQ);
    }

    private int countMaxOps(final Cigar cigar, final CigarOperator op) {
        final List<CigarElement> elements = cigar.getCigarElements();
        int maxOps = 0;
        for (int i = 0; i < elements.size(); i++) {
            int current;
            if (elements.get(i).getOperator() != op) {
                continue;
            }
            current = elements.get(i).getLength();
            final int next = i + 1;
            final int nextNext = next + 1;

            if (next < elements.size() &&
                    elements.get(next).getLength() == 1 &&
                    nextNext < elements.size() &&
                    elements.get(nextNext).getOperator() == op) {
                current += elements.get(nextNext).getLength();
            }
            maxOps = Math.max(maxOps, current);

        }
        return maxOps;
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();

        MetricsFile<?, Integer> metricsFile = getMetricsFile();
        cisMatchLength.setBinLabel("match_length");
        cisMatchLength.setValueLabel("cis_count");
        metricsFile.addHistogram(cisMatchLength);
        transMatchLength.setBinLabel("match_length");
        transMatchLength.setValueLabel("trans_count");
        metricsFile.addHistogram(transMatchLength);

        metricsFile.write(new File(outputMetrics));
        outputWriter.close();

        return null;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.addAll(super.getDefaultReadFilters());
        filters.add(new ReadFilterLibrary.NotDuplicateReadFilter());
        filters.add(new MappingQualityReadFilter(60));
        return filters;
    }
}
