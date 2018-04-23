package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * Clears the 0x400 duplicate SAM flag from reads.
 * <p>
 * Most GATK tools employ the NotDuplicateReadFilter that removes duplicate reads from analysis.
 * For these GATK tools, it is possible to disable the engine-level NotDuplicateReadFilter with the --disable-read-filter argument.
 * Disabling the filter allows a tool to then include duplicate reads in its analysis.
 * Certain data types, e.g. amplicon data, need to include reads flagged as duplicate in downstream analyses.
 * </p>
 *
 * <p>
 * Removing the duplicate flag in its entirety may be desirable for convenience
 * or for analysis with programs that do not allow disabling their duplicate read filter.
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 * <li>
 * A SAM/BAM/CRAM file marked for duplicates
 * </li>
 * </ul>
 * <h3>Output</h3>
 * <ul>
 * <li>
 * A SAM/BAM/CRAM file lacking 0x400 flags on reads
 * </li>
 * </ul>
 *
 * <h3>Usage example: </h3>
 * <pre>
 * gatk UnmarkDuplicates \
 *   -I clean.bam \
 *   -O unmarked.bam
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = UnmarkDuplicates.USAGE_DETAILS,
        oneLineSummary = UnmarkDuplicates.USAGE_SUMMARY,
        usageExample = "gatk UnmarkDuplicates -I marked_duplicates.bam -O unmarked_duplicates.bam",
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class UnmarkDuplicates extends ReadWalker {

    static final String USAGE_SUMMARY = "Clears the 0x400 duplicate SAM flag";
    static final String USAGE_DETAILS = "Simple tool to unmark duplicates in a SAM/BAM/CRAM file. "+
            "Clears the 0x400 SAM flag bit on all reads.";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public String OUTPUT;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(OUTPUT), true);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        read.setIsDuplicate(false);
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
