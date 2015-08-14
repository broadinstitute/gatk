package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.ApplyBQSRTransform;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalOutputSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Applies the BQSR table to the input BAM.",
        oneLineSummary = "Applies the BQSR table to the input BAM.",
        programGroup = ReadProgramGroup.class
)
public final class ApplyBQSRDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    private static final Logger logger = LogManager.getLogger(ApplyBQSRDataflow.class);

    @ArgumentCollection
    public final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    /**
     * Output path. Can be local or gs://
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write recalibrated reads to this new BAM file")
    public String OUTPUT;

    /**
     * Enables recalibration of base qualities.
     * The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Argument(fullName="bqsr_recal_file", shortName="bqsr", doc="Input covariates table file for base quality score recalibration")
    public String BQSR_RECAL_FILE_NAME;

    /**
     * command-line options to fine tune the recalibration.
     */
    @ArgumentCollection
    public ApplyBQSRArgumentCollection bqsrOpts = new ApplyBQSRArgumentCollection();

    private String intermediateRemoteBam;

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        if (readArguments.getReadFilesNames().size()>1) {
            throw new UserException("Sorry, we only support a single input file for now.");
        }
        String filename = readArguments.getReadFilesNames().get(0);
        ReadsDataflowSource readsSource = new ReadsDataflowSource(filename, pipeline);
        SAMFileHeader header = readsSource.getHeader();
        PCollectionView<SAMFileHeader> headerView = pipeline.apply(Create.of(header)).apply(View.asSingleton());

        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary) :
                IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        PCollectionView<BaseRecalOutput> recalInfoSingletonView = BaseRecalOutputSource.loadFileOrRemote(pipeline, BQSR_RECAL_FILE_NAME).apply(View.asSingleton());
        PCollection<GATKRead> output = readsSource.getReadPCollection(intervals, ValidationStringency.SILENT, false)
                .apply(new ApplyBQSRTransform(headerView, recalInfoSingletonView, bqsrOpts));
        intermediateRemoteBam = OUTPUT;
        if (needsIntermediateCopy()) {
            // The user specified remote execution and provided a local file name. So we're going to have to save to remote storage as a go-between.
            // Note that this may require more permissions
            intermediateRemoteBam = BucketUtils.randomRemotePath(stagingLocation, "temp-applyBqsr-output-", ".bam");
            logger.info("Staging results at " + intermediateRemoteBam);
        }
        SmallBamWriter.writeToFile(pipeline, output, header, intermediateRemoteBam);
    }

    @Override
    protected void afterPipeline(Pipeline pipeline) {
        if (!needsIntermediateCopy()) return;
        try {
            logger.info("Copying results from " + intermediateRemoteBam + " to " + OUTPUT + ".");
            BucketUtils.copyFile(intermediateRemoteBam, pipeline.getOptions(), OUTPUT);
        } catch (IOException x) {
            // keep the intermediate file if anything goes wrong, so we can investigate
            throw new UserException.CouldNotCreateOutputFile("Error writing to '" + OUTPUT + "'.",x);
        }
        try {
            BucketUtils.deleteFile(intermediateRemoteBam, pipeline.getOptions());
        } catch (Exception x) {
            // log error but continue since this error isn't fatal.
            logger.warn("Unable to delete temporary file '" + intermediateRemoteBam + "'.", x);
        }

    }

    // Specified remote execution and a local output.
    // The user probably didn't mean for the output to end up on the worker's local disk.
    private boolean needsIntermediateCopy() {
        return isRemote() && !BucketUtils.isRemoteStorageUrl(OUTPUT);
    }

}
