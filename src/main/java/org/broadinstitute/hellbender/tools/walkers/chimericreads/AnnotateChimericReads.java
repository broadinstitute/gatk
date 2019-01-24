package org.broadinstitute.hellbender.tools.walkers.chimericreads;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Prints chimeric reads from the provided file or files with corresponding reference bases
 * in an auxililary tag.
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = OtherProgramGroup.class
)
@ExperimentalFeature
public final class AnnotateChimericReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String OUTPUT_FILE = null;

    @SuppressWarnings("FieldCanBeLocal")
    @Argument(fullName = "reference_length", doc = "The desired length of the reference sequence to annotate chimeric reads with. " +
            "Reference sequence is taken from the 3' end of the read in the direction of the 5' end. " +
            "Unmapped reads are ignored.")
    private Integer REFERENCE_LENGTH = 400;

    @Argument(fullName = "filter_nonchimeric_reads", doc = "Should non-chimeric reads be filtered or not.")
    private Boolean FILTER_NONCHIMERIC_READS = Boolean.TRUE;

    @Argument(fullName = "max_read_pairs", doc = "Stop after this many chimeric read-pairs (mates of chimeras will be looked for and emitted). Set to 0 to process all the reads.")
    private Integer MAX_READ_PAIRS = 1_000_000;

    @Argument(fullName = "reference-bases-tag-name", shortName = "tn")
    String tagName = "rb";

    private SAMFileGATKReadWriter outputWriter;
    private final ReadFilter chimericReadFilter = new ChimericReadFilter();
    private final Set<String> readPairs = new HashSet<>(MAX_READ_PAIRS / 10);
    private int processedReads = 0;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultReadFilters = super.getDefaultReadFilters();
        if (FILTER_NONCHIMERIC_READS) {
            final List<ReadFilter> filters = new ArrayList<>(defaultReadFilters);
            filters.add(chimericReadFilter);
            return filters;
        }
        return defaultReadFilters;
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(OUTPUT_FILE), true);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        final boolean earlyBreak = MAX_READ_PAIRS > 0 && processedReads >= MAX_READ_PAIRS;
        final boolean addRead = MAX_READ_PAIRS > 0 && chimericReadFilter.test(read) && !readPairs.contains(read.getName());

        if (earlyBreak) {
            readPairs.remove(read.getName());
            if (addRead) {
                return;
            }
        } else {
            if (addRead) {
                readPairs.add(read.getName());
            }
        }

        if (chimericReadFilter.test(read)) {
            if (processedReads++ == MAX_READ_PAIRS) {
                logger.info("Processed requested chimeric reads. Continuing to read through input file to find all the mates.");
            }

            final int windowLength = referenceContext.getWindow().getLengthOnReference();
            referenceContext.setWindow(read.isReverseStrand() ? REFERENCE_LENGTH - windowLength : 0,
                    read.isReverseStrand() ? 0 : REFERENCE_LENGTH - windowLength);

            // useful to reverse-complement the reference to a standard form.
            final boolean orientedCorrectly = read.isFirstOfPair() ^ read.isReverseStrand();

            read.setAttribute(tagName, maybeRevComp(new String(referenceContext.getBases()), !orientedCorrectly));
        }
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if (outputWriter != null) {
            outputWriter.close();
        }
    }

    private String maybeRevComp(final String input, boolean revComp) {
        return revComp ? SequenceUtil.reverseComplement(input) : input;
    }

    private static class ChimericReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test(final GATKRead read) {
            return read != null &&
                    !read.isSecondaryAlignment() &&
                    !read.isSupplementaryAlignment() &&
                    read.isPaired() &&
                    !read.isProperlyPaired() &&
                    !read.isUnmapped() &&
                    !read.mateIsUnmapped() &&
                    !read.getContig().equals(read.getMateContig());
        }
    }
}
