package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.sun.security.auth.UserPrincipal;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
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
import org.broadinstitute.hellbender.dev.pipelines.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.BaseRecalOutputSource;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.SmallBamWriter;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.channels.Channels;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@CommandLineProgramProperties(
        usage = "Applies the BQSR table to the input BAM.",
        usageShort = "Applies the BQSR table to the input BAM.",
        programGroup = ReadProgramGroup.class
)
public final class ApplyBQSRDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    private final static Logger logger = LogManager.getLogger(ApplyBQSRDataflow.class);

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

    private String intermediateGCSBam;

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        if (readArguments.getReadFilesNames().size()>1) {
            throw new UserException("Sorry, we only support a single input file for now.");
        }
        String filename = readArguments.getReadFilesNames().get(0);
        ReadsSource readsSource = new ReadsSource(filename, pipeline);
        SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary) :
                IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        PCollection<BaseRecalOutput> recalInfoSingletonCollection = BaseRecalOutputSource.loadFileOrGcs(pipeline, BQSR_RECAL_FILE_NAME);
        PCollection<Read> output = readsSource.getReadPCollection(intervals, ValidationStringency.SILENT)
                .apply(new ApplyBQSRTransform(header, recalInfoSingletonCollection, bqsrOpts));
        intermediateGCSBam = OUTPUT;
        if (needsIntermediateCopy()) {
            // The user specified remote execution and provided a local file name. So we're going to have to save to GCS as a go-between.
            // Note that this may require more permissions
            intermediateGCSBam = BucketUtils.randomGcsPath(stagingLocation, "temp-applyBqsr-output-", ".bam");
            logger.info("Staging results at " + intermediateGCSBam);
        }
        SmallBamWriter.writeToFile(pipeline, output, header, intermediateGCSBam);
    }

    @Override
    protected void afterPipeline(Pipeline pipeline) {
        if (!needsIntermediateCopy()) return;
        try {
            logger.info("Copying results from " + intermediateGCSBam + " to " + OUTPUT + ".");
            BucketUtils.copyFile(intermediateGCSBam, pipeline.getOptions(), OUTPUT);
        } catch (IOException x) {
            // keep the intermediate file if anything goes wrong, so we can investigate
            throw new UserException.CouldNotCreateOutputFile("Error writing to '" + OUTPUT + "'.",x);
        }
        try {
            BucketUtils.deleteFile(intermediateGCSBam, pipeline.getOptions());
        } catch (Exception x) {
            // log error but continue since this error isn't fatal.
            logger.warn("Unable to delete temporary file '" + intermediateGCSBam + "'.", x);
        }

    }

    // Specified remote execution and a local output.
    // The user probably didn't mean for the output to end up on the worker's local disk.
    private boolean needsIntermediateCopy() {
        return isRemote() && !BucketUtils.isCloudStorageUrl(OUTPUT);
    }

}
