package org.broadinstitute.hellbender.tools.walkers.chimericreads;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Prints chimeric reads from the provided file or files with corresponding reference bases
 * in an auxililary tag.
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class AnnotateChimericReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private Path OUTPUT_FILE = null;

    @SuppressWarnings("FieldCanBeLocal")
    @Argument(fullName = "reference_length", doc = "The desired length of the reference sequence to annotate chimeric reads with. " +
            "Reference sequence is taken from the 3' end of the read in the direction of the 5' end. " +
            "Unmapped reads are ignored.")
    private Integer REFERENCE_LENGTH = 400;

    @Argument(fullName = "filter_nonchimeric_reads", doc = "Should non-chimeric reads be filtered or not.")
    private Boolean FILTER_NONCHIMERIC_READS = Boolean.TRUE;

    @Argument(fullName = "max_read_pairs", doc = "Stop after this many chimeric read-pairs (mates of chimeras will be looked for and emitted). Set to 0 to process all the reads.")
    private Integer MAX_READ_PAIRS = 1_000_000;

    private SAMFileGATKReadWriter outputWriter;
    private final ReadFilter chimericReadFilter = new ChimericReadFilter();
    private final Set<String> readPairs = new HashSet<>(MAX_READ_PAIRS/10);
    private int processedReads=0;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultReadFilters = super.getDefaultReadFilters();
        if (FILTER_NONCHIMERIC_READS) {
            defaultReadFilters.add(chimericReadFilter);
        }
        return defaultReadFilters;
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT_FILE, true);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        final boolean earlyBreak = MAX_READ_PAIRS > 0 && processedReads > MAX_READ_PAIRS;
        final boolean addRead = MAX_READ_PAIRS > 0 && !readPairs.contains(read.getName());

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
            processedReads++;
            referenceContext.setWindow(read.isReverseStrand() ? REFERENCE_LENGTH - referenceContext.getBases().length : 0,
                    read.isReverseStrand() ? 0 : REFERENCE_LENGTH - referenceContext.getBases().length);
            read.setAttribute("RB", referenceContext.getBases());

        }
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if (outputWriter != null) {
            outputWriter.close();
        }
    }

    private static class ChimericReadFilter extends ReadFilter {
        @Override
        public boolean test(final GATKRead read) {
            return !read.isProperlyPaired() && !read.getContig().equals(read.getMateContig());
        }
    }
}
