package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Simple tool to \"unmark\" duplicates in a SAM/BAM file. Clears the isDuplicate bit on all reads.",
        oneLineSummary = "Unmark Duplicates",
        usageExample = "hellbender UnmarkDuplicates -I input.bam -O output.bam",
        programGroup = ReadProgramGroup.class
)
public class UnmarkDuplicates extends ReadWalker {

    @Argument(fullName = "output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        final SAMFileHeader outputHeader = ReadUtils.cloneSAMFileHeader(getHeaderForReads());
        outputWriter = new SAMFileGATKReadWriter(new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, referenceArguments.getReferenceFile()));
    }

    @Override
    public CountingReadFilter makeReadFilter() {
        return new CountingReadFilter("Allow all reads", ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        read.setIsDuplicate(false);
        outputWriter.addRead(read);
    }

    @Override
    public Object onTraversalDone() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
        return null;
    }
}
