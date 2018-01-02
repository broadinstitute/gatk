package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
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
 * Removes duplicate read tags.
 *
 * <p> This tool "unmarks" duplicate reads in a SAM/BAM/CRAM file. Duplicate reads are defined as originated from a
 * single fragment of DNA and can arise during sample preparation. </p>
 *
 *
 * <h3>Input</h3>
 * <p>
 * A SAM/BAM/CRAM file that has been marked for duplicates.
 * </p>
 * <h3>Output</h3>
 * <p>
 * A SAM/BAM/CRAM file without isDuplicate on reads.
 * </p>
 *
 * <h3>Usage example: </h3>
 * <pre>
 *     gatk UnmarkDuplicates \
 *     -I marked_duplicates.bam \
 *     -O unmarked_duplicates.bam
 * </pre>
 *
 */


@CommandLineProgramProperties(
        summary = UnmarkDuplicates.USAGE_DETAILS,
        oneLineSummary = UnmarkDuplicates.USAGE_SUMMARY,
        usageExample = "gatk UnmarkDuplicates -I marked_duplicates.bam -O unmarked_duplicates.bam",
        programGroup = ReadProgramGroup.class)
@DocumentedFeature
public class UnmarkDuplicates extends ReadWalker {

    static final String USAGE_SUMMARY = "Unmark duplicates in a SAM/BAM/CRAM file";
    static final String USAGE_DETAILS = "Simple tool to \"unmark\" duplicates in a SAM/BAM/CRAM file."+
            "Clears the isDuplicate bit on all reads.";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT, true);
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
