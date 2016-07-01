package org.broadinstitute.hellbender.tools.spark;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.base.Stopwatch;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.ContextShard;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.VariantsSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSparkOptimized;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.transforms.bqsr.BaseRecalibratorEngineSparkWrapper;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

@CommandLineProgramProperties(
        summary = "Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "BaseRecalibrator on Spark (experimental sharded implementation)",
        programGroup = SparkProgramGroup.class
)
public class BaseRecalibratorSparkSharded extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    // inputs can be local, GCS, or HDFS.
    @ArgumentCollection
    private final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    private final IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    private final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();

    @Argument(doc = "the known variants. Must be local.", shortName = "knownSites", fullName = "knownSites", optional = false)
    private List<String> knownVariants;

    // output can be local or GCS.
    @Argument(doc = "Path to save the final recalibration tables to.",
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputTablesPath = null;


    private AuthHolder auth;

    @Override
    protected void runPipeline( JavaSparkContext ctx ) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for now.");
        }
        final String bam = readArguments.getReadFilesNames().get(0);
        final String referenceURL = referenceArguments.getReferenceFileName();

        auth = getAuthHolder();

        final ReferenceMultiSource rds = new ReferenceMultiSource(auth, referenceURL, BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION);

        SAMFileHeader readsHeader = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency()).getHeader(bam, referenceURL, auth);
        final SAMSequenceDictionary readsDictionary = readsHeader.getSequenceDictionary();
        final SAMSequenceDictionary refDictionary = rds.getReferenceSequenceDictionary(readsDictionary);
        final ReadFilter readFilterToApply = BaseRecalibrator.getStandardBQSRReadFilter(readsHeader);

        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);

        Broadcast<SAMFileHeader> readsHeaderBcast = ctx.broadcast(readsHeader);
        Broadcast<SAMSequenceDictionary> refDictionaryBcast = ctx.broadcast(refDictionary);

        List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        List<String> localVariants = knownVariants;
        localVariants = hackilyCopyFromGCSIfNecessary(localVariants);
        List<GATKVariant> variants = VariantsSource.getVariantsList(localVariants);

        // get reads, reference, variants
        JavaRDD<ContextShard> readsWithContext = AddContextDataToReadSparkOptimized.add(ctx, intervals, bam, variants, auth, readFilterToApply, rds);

        // run BaseRecalibratorEngine.
        BaseRecalibratorEngineSparkWrapper recal = new BaseRecalibratorEngineSparkWrapper(readsHeaderBcast, refDictionaryBcast, bqsrArgs);
        JavaRDD<RecalibrationTables> tables = readsWithContext.mapPartitions(s->recal.apply(s));

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(bqsrArgs, readsHeader));
        final RecalibrationTables table = tables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(tables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(table);

        try {
            BaseRecalibratorEngineSparkWrapper.saveTextualReport(outputTablesPath, readsHeader, table, bqsrArgs, auth);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(new File(outputTablesPath), e);
        }
    }


    // please add support for reading variant files from GCS.
    private ArrayList<String> hackilyCopyFromGCSIfNecessary(List<String> localVariants) {
        int i=0;
        Stopwatch hacking = Stopwatch.createStarted();
        boolean copied = false;
        ArrayList<String> ret = new ArrayList<>();
        for (String v : localVariants) {
            if (BucketUtils.isCloudStorageUrl(v)) {
                if (!copied) {
                    logger.info("(HACK): copying the GCS variant file to local just so we can read it back.");
                    copied=true;
                }
                // this only works with the API_KEY, but then again it's a hack so there's no point in polishing it. Please don't make me.
                PipelineOptions popts = auth.asPipelineOptionsDeprecated();
                String d = IOUtils.createTempFile("knownVariants-"+i,".vcf").getAbsolutePath();
                try {
                    BucketUtils.copyFile(v, popts, d);
                } catch (IOException x) {
                    throw new UserException.CouldNotReadInputFile(v,x);
                }
                ret.add(d);
            } else {
                ret.add(v);
            }
        }
        hacking.stop();
        if (copied) {
            logger.info("Copying the vcf took "+hacking.elapsed(TimeUnit.MILLISECONDS)+" ms.");
        }
        return ret;
    }

}
