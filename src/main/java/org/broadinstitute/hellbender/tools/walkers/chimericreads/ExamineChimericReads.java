package org.broadinstitute.hellbender.tools.walkers.chimericreads;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignerUtils;
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
            final byte[] seq1 = refOne.getBytes();
            final byte[] seq2 = refTwo.getBytes();

            alignAndInfo(seq1, seq2, "ref1:ref2: ");
        }
        {
            final byte[] seq1 = refOne.getBytes();
            final byte[] seq2 = SequenceUtil.reverseComplement(refTwo).getBytes();

            alignAndInfo(seq1, seq2, "ref1:ref2': ");
        }

//        {
//            final byte[] seq1 = refTwo.getBytes();
//            final byte[] seq2 = readOne.getBases();
//
//            alignAndInfo(seq1, seq2, "ref2:read1: ");
//        }
//        {
//            final byte[] seq1 = refTwo.getBytes();
//            final byte[] seq2 = SequenceUtil.reverseComplement(readOne.getBasesString()).getBytes();
//
//            alignAndInfo(seq1, seq2, "ref2:read1': ");
//        }
//        {
//            final byte[] seq1 = refOne.getBytes();
//            final byte[] seq2 = SequenceUtil.reverseComplement(readTwo.getBasesString()).getBytes();
//
//            alignAndInfo(seq1, seq2, "ref1:read2': ");
//        }
//        {
//            final byte[] seq1 = refTwo.getBytes();
//            final byte[] seq2 = SequenceUtil.reverseComplement(readOne.getBasesString()).getBytes();
//
//            alignAndInfo(seq1, seq2, "ref2:read1': ");
//        }
    }

    private void alignAndInfo(final byte[] seq1, final byte[] seq2, final String what) {
        final SmithWatermanAlignment alignmentMis = smithWatermanAligner.alignWithMismatches(seq1, seq2, swParameters, swOverhangStrategy);
        final SmithWatermanAlignment alignment = smithWatermanAligner.align(seq1, seq2, swParameters, swOverhangStrategy);

        logger.info(what + alignmentMis.getCigar().toString());
        logger.info(alignment.getCigar().toString());
        logger.info(String.format("%d matches, and %d equals, rate= %g", countMatches(alignment.getCigar()),
                countEquals(alignmentMis.getCigar()),
                countEquals(alignmentMis.getCigar()) / (double) countMatches(alignment.getCigar())));
        SmithWatermanAlignerUtils.printAlignment(seq1, seq2, alignmentMis, swOverhangStrategy);
    }

    private int countMatches(final Cigar cigar) {
        return countOps(cigar,CigarOperator.M);
    }

    private int countEquals(final Cigar cigar){
        return countOps(cigar,CigarOperator.EQ);
    }

    private int countOps(final Cigar cigar, final CigarOperator op){
        return cigar.getCigarElements().stream()
                .filter(e->e.getOperator()== op)
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.addAll(super.getDefaultReadFilters());
        filters.add(new MappingQualityReadFilter(60));
        return filters;
    }
}
