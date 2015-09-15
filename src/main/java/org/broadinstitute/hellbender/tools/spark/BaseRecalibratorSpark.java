package org.broadinstitute.hellbender.tools.spark;

import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.dataflow.pipelines.BaseRecalibratorDataflow;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "Generates recalibration table",
        programGroup = ReadProgramGroup.class
)
public class BaseRecalibratorSpark extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final BaseRecalibrationArgumentCollection bqsrArgs = new BaseRecalibrationArgumentCollection();

    @Argument(doc = "the known variants", shortName = "BQSRKnownVariants", fullName = "baseRecalibrationKnownVariants", optional = false)
    protected List<String> baseRecalibrationKnownVariants;

    /**
     * Path to save the serialized tables to. Local or GCS.
     */
    @Argument(doc = "Path to save the serialized recalibrationTables to. If running on the cloud, either leave outputTablesPath unset or point it to a GCS location.",
            shortName = "O", fullName = "outputTablesPath", optional = true)
    protected String outputTablesPath = null;

    @Override
    protected void runPipeline( JavaSparkContext ctx ) {
        if ( bqsrArgs.readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for now.");
        }
        final String bam = bqsrArgs.readArguments.getReadFilesNames().get(0);

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        final List<SimpleInterval> intervals = bqsrArgs.intervalArgumentCollection.intervalsSpecified() ? bqsrArgs.intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        JavaRDD<GATKRead> initialReads = readSource.getParallelReads(bam, intervals);

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<Variant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(baseRecalibrationKnownVariants);

        GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
        options.setApiKey(apiKey);
        final String referenceURL = bqsrArgs.referenceArguments.getReferenceFileName();
        final ReferenceDataflowSource referenceDataflowSource = new ReferenceDataflowSource(options, referenceURL, BaseRecalibratorDataflow.BQSR_REFERENCE_WINDOW_FUNCTION);
        //checkSequenceDictionaries(referenceDataflowSource.getReferenceSequenceDictionary(null), readsHeader.getSequenceDictionary());

        // TODO: Look into broadcasting the reference to all of the workers. This would make AddContextDataToReadSpark
        // TODO: and ApplyBQSRStub simpler (#855).
        JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(initialReads, referenceDataflowSource, bqsrKnownVariants);
        // TODO: broadcast the reads header?
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, readsHeader, referenceDataflowSource.getReferenceSequenceDictionary(null), bqsrArgs);

        try ( final PrintStream reportStream = new PrintStream(BucketUtils.createFile(outputTablesPath, options)) ) {
            bqsrArgs.RAC.RECAL_TABLE = reportStream;
            bqsrArgs.RAC.RECAL_TABLE_FILE = new File(outputTablesPath);
            RecalUtils.outputRecalibrationReport(bqsrArgs.RAC, bqsrReport.getQuantizationInfo(), bqsrReport.getRecalibrationTables(), new StandardCovariateList(bqsrArgs.RAC, readsHeader), true);
        }
    }


    @Override
    protected String getProgramName() {
        return "BaseRecalibrator";
    }

    private void checkSequenceDictionaries(final SAMSequenceDictionary refDictionary, SAMSequenceDictionary readsDictionary) {
        Utils.nonNull(refDictionary);
        Utils.nonNull(readsDictionary);
        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);
    }
}
