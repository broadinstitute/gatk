package org.broadinstitute.hellbender.tools.walkers.chimericreads;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.ReadPairWalker;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;

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
    private final SWOverhangStrategy swOverhangStrategy = SWOverhangStrategy.LEADING_INDEL;

    private SWParameters swParameters = null;

    @Argument(fullName = "match_value")
    int MATCH_VALUE = SmithWatermanAligner.STANDARD_NGS.getMatchValue();
    @Argument(fullName = "mismatch_value")
    int MISMATCH_PENALTY_VALUE = SmithWatermanAligner.STANDARD_NGS.getMismatchPenalty();
    @Argument(fullName = "gap_open_value")
    int GAP_OPEN_VALUE = SmithWatermanAligner.STANDARD_NGS.getGapOpenPenalty();
    @Argument(fullName = "gap_extend_value")
    int GAP_EXTEND_VALUE = SmithWatermanAligner.STANDARD_NGS.getGapExtendPenalty();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        swParameters = new SWParameters(MATCH_VALUE, MISMATCH_PENALTY_VALUE, GAP_OPEN_VALUE, GAP_EXTEND_VALUE);
    }

    @Override
    public void apply(final Set<GATKRead> reads) {
        final GATKRead readOne = reads.stream().filter(GATKRead::isFirstOfPair).findFirst().orElse(null);
        final GATKRead readTwo = reads.stream().filter(gatkRead -> !gatkRead.isFirstOfPair()).findFirst().orElse(null);

        if (readOne == null || readTwo == null ||
                !readOne.hasAttribute("RB") ||
                !readTwo.hasAttribute("RB")) {
            return;
        }

        logger.info(readOne.getName());
        final String refOne = readOne.getAttributeAsString("RB");
        final String refTwo = readTwo.getAttributeAsString("RB");

        {
            final SmithWatermanAlignment alignment = smithWatermanAligner.alignWithMismatches(refOne.getBytes(), refTwo.getBytes(),
                    swParameters, swOverhangStrategy);

            logger.info("ref1:ref2: "+ alignment.getCigar().toString());
        }
        {
            final SmithWatermanAlignment alignment = smithWatermanAligner.alignWithMismatches(refOne.getBytes(), SequenceUtil.reverseComplement(refTwo).getBytes(),
                    swParameters, swOverhangStrategy);

            logger.info("ref1:ref2': "+alignment.getCigar().toString());
        }
        {
            final SmithWatermanAlignment alignment = smithWatermanAligner.alignWithMismatches(refTwo.getBytes(), readOne.getBases(),
                    swParameters, swOverhangStrategy);

            logger.info("ref2:read1: " + alignment.getCigar().toString());
        }
        {
            final SmithWatermanAlignment alignment = smithWatermanAligner.alignWithMismatches(refTwo.getBytes(), SequenceUtil.reverseComplement(readOne.getBasesString()).getBytes(),
                    swParameters, swOverhangStrategy);

            logger.info("ref2:read1': " + alignment.getCigar().toString());
        }
        {
            final SmithWatermanAlignment alignment = smithWatermanAligner.alignWithMismatches(refOne.getBytes(), readTwo.getBases(),
                    swParameters, swOverhangStrategy);

            logger.info("ref1:read2: " + alignment.getCigar().toString());
        }
        {
            final SmithWatermanAlignment alignment = smithWatermanAligner.align(refTwo.getBytes(), SequenceUtil.reverseComplement(readOne.getBasesString()).getBytes(),
                    swParameters, swOverhangStrategy);

            logger.info("ref1:read2': " + alignment.getCigar().toString());
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.addAll(super.getDefaultReadFilters());
        filters.add(new MappingQualityReadFilter(60));
        return filters;
    }
}
