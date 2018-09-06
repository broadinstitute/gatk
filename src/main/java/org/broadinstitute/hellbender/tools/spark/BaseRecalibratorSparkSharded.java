package org.broadinstitute.hellbender.tools.spark;

import com.google.common.base.Stopwatch;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
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
import org.broadinstitute.hellbender.utils.recalibration.covariates.BQSRCovariateList;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Experimental sharded spark implementation of the first pass of the base quality score recalibration.
 * Generates a recalibration table based on various covariates.
 * The default covariates are read group, reported quality score, machine cycle, and nucleotide context.
 *
 *
 * <h3>Input</h3>
 * <p>
 * The input read data whose base quality scores need to be assessed.
 * <p>
 * A database of known polymorphic sites to skip over.
 * </p>
 * <p>
 * <h3>Output</h3>
 * <p>
 * A GATK Report file with many tables:
 * <ol>
 * <li>The list of arguments</li>
 * <li>The quantized qualities table</li>
 * <li>The recalibration table by read group</li>
 * <li>The recalibration table by quality score</li>
 * <li>The recalibration table for all the optional covariates</li>
 * </ol>
 * <p>
 * The GATK Report is intended to be easy to read by humans or computers. Check out the documentation of the GATKReport to learn how to manipulate this table.
 * </p>
 * <p>
 * <h3>Examples</h3>
 * <pre>
 * gatk BaseRecalibratorSparkSharded \
 *   -I gs://my-gcs-bucket/my_reads.bam \
 *   -R gs://my-gcs-bucket/reference.fasta \
 *   --known-sites gs://my-gcs-bucket/sites_of_variation.vcf \
 *   --known-sites gs://my-gcs-bucket/another/optional/setOfSitesToMask.vcf \
 *   -O gs://my-gcs-bucket/recal_data.table \
 *   -- \
 *   --sparkRunner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 */

@CommandLineProgramProperties(
        summary = BaseRecalibratorSparkSharded.USAGE_SUMMARY,
        oneLineSummary = BaseRecalibratorSparkSharded.USAGE_ONE_LINE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class BaseRecalibratorSparkSharded extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;
    static final String USAGE_ONE_LINE_SUMMARY = "BaseRecalibrator on Spark (experimental sharded implementation)";
    static final String USAGE_SUMMARY = "Experimental sharded implementation of the first pass of the Base Quality Score Recalibration (BQSR)" +
            " -- Generates recalibration table based on various user-specified covariates " +
            "(such as read group, reported quality score, machine cycle, and nucleotide context).";
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

    @Argument(doc = "the known variants. Must be local.", fullName = BaseRecalibrator.KNOWN_SITES_ARG_FULL_NAME, optional = false)
    private List<String> knownVariants;

    // output can be local or GCS.
    @Argument(doc = "Path to save the final recalibration tables to.",
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputTablesPath = null;

    @Override
    protected void runPipeline( JavaSparkContext ctx ) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for now.");
        }
        final String bam = readArguments.getReadFilesNames().get(0);
        final String referenceURL = referenceArguments.getReferenceFileName();

        final ReferenceMultiSource rds = new ReferenceMultiSource(referenceURL, BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION);

        SAMFileHeader readsHeader = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency()).getHeader(bam, referenceURL);
        final SAMSequenceDictionary readsDictionary = readsHeader.getSequenceDictionary();
        final SAMSequenceDictionary refDictionary = rds.getReferenceSequenceDictionary(readsDictionary);
        final ReadFilter readFilterToApply = ReadFilter.fromList(BaseRecalibrator.getStandardBQSRReadFilterList(), readsHeader);

        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);

        Broadcast<SAMFileHeader> readsHeaderBcast = ctx.broadcast(readsHeader);
        Broadcast<SAMSequenceDictionary> refDictionaryBcast = ctx.broadcast(refDictionary);

        List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        List<String> localVariants = knownVariants;
        localVariants = hackilyCopyFromGCSIfNecessary(localVariants);
        List<GATKVariant> variants = VariantsSource.getVariantsList(localVariants);

        // get reads, reference, variants
        JavaRDD<ContextShard> readsWithContext = AddContextDataToReadSparkOptimized.add(ctx, intervals, bam, variants,
            readFilterToApply, rds);

        // run BaseRecalibratorEngine.
        BaseRecalibratorEngineSparkWrapper recal = new BaseRecalibratorEngineSparkWrapper(readsHeaderBcast, refDictionaryBcast, bqsrArgs);
        JavaRDD<RecalibrationTables> tables = readsWithContext.mapPartitions(s->recal.apply(s));

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new BQSRCovariateList(bqsrArgs, readsHeader));
        final RecalibrationTables table = tables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(tables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(table);

        try {
            BaseRecalibratorEngineSparkWrapper.saveTextualReport(outputTablesPath, readsHeader, table, bqsrArgs);
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
                String d = IOUtils.createTempFile("knownVariants-"+i,".vcf").getAbsolutePath();
                try {
                    BucketUtils.copyFile(v, d);
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
