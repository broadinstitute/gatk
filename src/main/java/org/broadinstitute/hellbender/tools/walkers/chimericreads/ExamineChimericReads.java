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
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = OtherProgramGroup.class,
        omitFromCommandLine = true
)
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
            doc="Write output to this file")
    public String outputBam;

    @Argument(fullName = StandardArgumentDefinitions.METRICS_FILE_LONG_NAME,
            shortName = StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME,
            doc="Write output to this file")
    public String outputMetrics;

    @Argument(fullName = "reference-bases-tag-name", shortName = "tn")
    String tagName = "RB";

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        swParameters = new SWParameters(MATCH_VALUE, MISMATCH_PENALTY_VALUE, GAP_OPEN_VALUE, GAP_EXTEND_VALUE);
        outputWriter = createSAMWriter(IOUtils.getPath(outputBam), false);
    }

    final private Histogram<Integer> minMatchLength = new Histogram<>();
    final private Histogram<Integer> maxMatchLength = new Histogram<>();
    private SAMFileGATKReadWriter outputWriter;

    @Override
    protected SAMFileHeader getHeaderForSAMWriter() {
        final SAMFileHeader headerForSAMWriter = super.getHeaderForSAMWriter();

        headerForSAMWriter.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return headerForSAMWriter;
    }

    @Override
    public void apply(final Set<GATKRead> reads) {

        final Map<Boolean, List<GATKRead>> firstOrSecond = reads.stream()
                .limit(2)
                .collect(Collectors.partitioningBy(GATKRead::isFirstOfPair));

        final GATKRead readOne = firstOrSecond.get(true).get(0);
        final GATKRead readTwo = firstOrSecond.get(false).get(0);

        if (readOne == null || readTwo == null ||
                !readOne.hasAttribute(tagName) ||
                !readTwo.hasAttribute(tagName)) {
            return;
        }

        final String refOne = readOne.getAttributeAsString(tagName);
        final String refTwo = readTwo.getAttributeAsString(tagName);

        final int matches, revCompMatches;
        {
            final byte[] seq1 = refOne.getBytes();
            final byte[] seq2 = refTwo.getBytes();
            final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);

            matches = countEquals(alignmentMis.getCigar());
        }
        {
            final byte[] seq1 = refOne.getBytes();
            final byte[] seq2 = SequenceUtil.reverseComplement(refTwo).getBytes();
            final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);

            revCompMatches = countEquals(alignmentMis.getCigar());
        }
        minMatchLength.increment(Math.min(matches, revCompMatches));
        maxMatchLength.increment(Math.max(matches, revCompMatches));

        readOne.setAttribute("mi", Math.min(matches, revCompMatches));
        readTwo.setAttribute("mi", Math.min(matches, revCompMatches));
        readOne.setAttribute("ma", Math.max(matches, revCompMatches));
        readTwo.setAttribute("ma", Math.max(matches, revCompMatches));

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

        MetricsFile<?,Integer> metricsFile = getMetricsFile();
        metricsFile.addHistogram(minMatchLength);
        metricsFile.addHistogram(maxMatchLength);

        metricsFile.write(new File(outputMetrics));

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
