package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Applies the BQSR table to the input BAM.",
        oneLineSummary = "Applies the BQSR table to the input BAM.",
        programGroup = ReadProgramGroup.class
)
public final class ApplyBQSR extends ReadWalker{

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    /**
     * Enables recalibration of base qualities.
     * The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Argument(fullName="bqsr_recal_file", shortName="bqsr", doc="Input covariates table file for base quality score recalibration")
    public File BQSR_RECAL_FILE;

    /**
     * command-line arguments to fine tune the recalibration.
     */
    @ArgumentCollection
    public ApplyBQSRArgumentCollection bqsrArgs = new ApplyBQSRArgumentCollection();
    
    private SAMFileGATKReadWriter outputWriter;

    private ReadTransformer transform;

    @Override
    public void onTraversalStart() {
        final SAMFileHeader outputHeader = ReadUtils.cloneSAMFileHeader(getHeaderForReads());
        outputWriter = new SAMFileGATKReadWriter(new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, referenceArguments.getReferenceFile()));
        transform = new BQSRReadTransformer(outputHeader, BQSR_RECAL_FILE, bqsrArgs.quantizationLevels, bqsrArgs.disableIndelQuals, bqsrArgs.PRESERVE_QSCORES_LESS_THAN, bqsrArgs.emitOriginalQuals, bqsrArgs.globalQScorePrior);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(transform.apply(read));
    }

    @Override
    public Object onTraversalDone() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
        return null;
    }

}
