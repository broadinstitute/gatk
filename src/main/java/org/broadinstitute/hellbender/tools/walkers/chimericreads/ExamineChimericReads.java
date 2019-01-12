package org.broadinstitute.hellbender.tools.walkers.chimericreads;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.ReadPairWalker;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadWriter;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = ExampleProgramGroup.class,
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
    public String output;

    private GATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        swParameters = new SWParameters(MATCH_VALUE, MISMATCH_PENALTY_VALUE, GAP_OPEN_VALUE, GAP_EXTEND_VALUE);
            outputWriter = createSAMWriter(IOUtils.getPath(output), true);
    }

    final private Histogram<Integer> minMatchLength=new Histogram<>();
    final private Histogram<Integer> maxMatchLength=new Histogram<>();

    @Override
    public void apply(final Set<GATKRead> reads) {
        final GATKRead readOne = reads.stream().filter(GATKRead::isFirstOfPair).findFirst().orElse(null);
        final GATKRead readTwo = reads.stream().filter(gatkRead -> !gatkRead.isFirstOfPair()).findFirst().orElse(null);

        if (readOne == null || readTwo == null ||
                !readOne.hasAttribute("RB") ||
                !readTwo.hasAttribute("RB")) {
            return;
        }

        final String refOne = readOne.getAttributeAsString("RB");
        final String refTwo = readTwo.getAttributeAsString("RB");

        final int matches, revCompMatches;
        {
            final byte[] seq1 = refOne.getBytes();
            final byte[] seq2 = refTwo.getBytes();

            final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);

            matches=countEquals(alignmentMis.getCigar());
        }
        {
            final byte[] seq1 = refOne.getBytes();
            final byte[] seq2 = SequenceUtil.reverseComplement(refTwo).getBytes();

            final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);

            revCompMatches=countEquals(alignmentMis.getCigar());
        }
        minMatchLength.increment(Math.min(matches, revCompMatches));
        maxMatchLength.increment(Math.max(matches, revCompMatches));
    }

    private int countEquals(final Cigar cigar) {
        return countMaxOps(cigar, CigarOperator.EQ);
    }

    private int countMaxOps(final Cigar cigar, final CigarOperator op) {
        final List<CigarElement> elements = cigar.getCigarElements();
        int maxOps = 0;
        for (int i = 0; i < elements.size(); i++) {
            int current = 0;
            if (elements.get(i).getOperator() != op) {
                continue;
            }
            current = elements.get(i).getLength();

            if (i + 1 < elements.size() && elements.get(i + 1).getLength() == 1 &&
                    i + 2 < elements.size() && elements.get(i + 2).getOperator() == op) {
                current += elements.get(i + 2).getLength();
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

        metricsFile.write(new File(output));

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
