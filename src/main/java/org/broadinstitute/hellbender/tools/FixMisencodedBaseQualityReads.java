package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.transformers.MisencodedBaseQualityReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;

import java.io.File;

@CommandLineProgramProperties(
	usage = "Fixes Illumina base quality scores.",
    usageShort = "Fix Illumina base quality scores.",
    programGroup = ReadProgramGroup.class
)
public class FixMisencodedBaseQualityReads extends ReadWalker {

    @Argument(fullName = "output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private SAMFileWriter outputWriter;

    private ReadTransformer transform;

    @Override
    public void onTraversalStart() {
        final SAMFileHeader outputHeader = getHeaderForReads().clone();
        outputWriter = new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, referenceArguments.getReferenceFile());
        transform = new MisencodedBaseQualityReadTransformer();
    }

    @Override
    public void apply( SAMRecord read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addAlignment(transform.apply(read));
    }

    @Override
    public Object onTraversalDone() {
        CloserUtil.close(outputWriter);
        return null;
    }
}
