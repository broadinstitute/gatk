package org.broadinstitute.hellbender.tools.spark;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.PrintStream;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "Generates recalibration table",
        programGroup = SparkProgramGroup.class
)
public class BaseRecalibratorSpark extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    @ArgumentCollection
    private final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    private final IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    private final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();

    @Argument(doc = "the known variants", shortName = "knownSites", fullName = "knownSites", optional = false)
    private List<String> knownVariants;

    @Argument(doc = "the join strategy for reference bases", shortName = "joinStrategy", fullName = "joinStrategy", optional = true)
    private JoinStrategy joinStrategy = JoinStrategy.SHUFFLE;

    @Argument(doc = "Path to save the final recalibration tables to.",
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputTablesPath = null;

    @Override
    protected void runPipeline( JavaSparkContext ctx ) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for now.");
        }
        final String bam = readArguments.getReadFilesNames().get(0);

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        JavaRDD<GATKRead> initialReads = readSource.getParallelReads(bam, intervals);

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        if ( knownVariants.size() > 1 ) {
            throw new GATKException("Cannot currently handle more than one known sites file, " +
                    "as getParallelVariants(List) is broken");
        }
        JavaRDD<Variant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(knownVariants.get(0));

        final GCSOptions gcsOptions = getAuthenticatedGCSOptions(); // null if we have no api key
        final String referenceURL = referenceArguments.getReferenceFileName();
        final ReferenceDataflowSource referenceDataflowSource = new ReferenceDataflowSource(gcsOptions, referenceURL, BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION);
        final SAMSequenceDictionary referenceDictionary = referenceDataflowSource.getReferenceSequenceDictionary(readsHeader.getSequenceDictionary());
        checkSequenceDictionaries(referenceDictionary, readsHeader.getSequenceDictionary());

        // TODO: Look into broadcasting the reference to all of the workers. This would make AddContextDataToReadSpark
        // TODO: and ApplyBQSRStub simpler (#855).
        JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(
                initialReads, referenceDataflowSource, bqsrKnownVariants, joinStrategy);
        // TODO: broadcast the reads header?
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, readsHeader, referenceDictionary, bqsrArgs);

        try ( final PrintStream reportStream = new PrintStream(BucketUtils.createFile(outputTablesPath, gcsOptions)) ) {
            RecalUtils.outputRecalibrationReport(reportStream, bqsrArgs, bqsrReport.getQuantizationInfo(), bqsrReport.getRecalibrationTables(), new StandardCovariateList(bqsrArgs, readsHeader), true);
        }
    }

    private void checkSequenceDictionaries(final SAMSequenceDictionary refDictionary, SAMSequenceDictionary readsDictionary) {
        Utils.nonNull(refDictionary);
        Utils.nonNull(readsDictionary);
        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);
    }
}
