package org.broadinstitute.hellbender.tools.spark.transforms;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;

public final class BaseRecalibratorSparkFn {

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary, final RecalibrationArgumentCollection recalArgs ) {
        JavaRDD<RecalibrationTables> unmergedTables = readsWithContext.mapPartitions(readWithContextIterator -> {
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);
            bqsr.logCovariatesUsed();

            while ( readWithContextIterator.hasNext() ) {
                final Tuple2<GATKRead, ReadContextData> readWithData = readWithContextIterator.next();
                Iterable<GATKVariant> variants = readWithData._2().getOverlappingVariants();
                final ReferenceBases refBases = readWithData._2().getOverlappingReferenceBases();
                ReferenceDataSource refDS = new ReferenceMemorySource(refBases, referenceDictionary);

                bqsr.processRead(readWithData._1(), refDS, variants);
            }
            return Arrays.asList(bqsr.getRecalibrationTables()).iterator();
        });

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(recalArgs, header));
        final RecalibrationTables combinedTables = unmergedTables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(unmergedTables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        final QuantizationInfo quantizationInfo = new QuantizationInfo(combinedTables, recalArgs.QUANTIZING_LEVELS);

        final StandardCovariateList covariates = new StandardCovariateList(recalArgs, header);
        return RecalUtils.createRecalibrationReport(recalArgs.generateReportTable(covariates.covariateNames()), quantizationInfo.generateReportTable(), RecalUtils.generateReportTables(combinedTables, covariates));
    }

    /**
     * Run the {@link BaseRecalibrationEngine} on reads and overlapping variants.
     * @param readsWithVariants the RDD of reads with overlapping variants
     * @param header the reads header
     * @param referenceFileName the name of the reference file added via {@code SparkContext#addFile()}
     * @param recalArgs arguments to use during recalibration
     * @return the recalibration report object
     */
    public static RecalibrationReport apply(final JavaPairRDD<GATKRead, Iterable<GATKVariant>> readsWithVariants, final SAMFileHeader header, final String referenceFileName, final RecalibrationArgumentCollection recalArgs) {
        JavaRDD<RecalibrationTables> unmergedTables = readsWithVariants.mapPartitions(readsWithVariantsIterator -> {
            String pathOnExecutor = SparkFiles.get(referenceFileName);
            ReferenceDataSource referenceDataSource = new ReferenceFileSource(IOUtils.getPath(pathOnExecutor));
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);
            bqsr.logCovariatesUsed();
            Utils.stream(readsWithVariantsIterator).forEach(t -> bqsr.processRead(t._1, referenceDataSource, t._2));
            return Iterators.singletonIterator(bqsr.getRecalibrationTables());
        });

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(recalArgs, header));
        final RecalibrationTables combinedTables = unmergedTables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(unmergedTables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        final QuantizationInfo quantizationInfo = new QuantizationInfo(combinedTables, recalArgs.QUANTIZING_LEVELS);

        final StandardCovariateList covariates = new StandardCovariateList(recalArgs, header);
        return RecalUtils.createRecalibrationReport(recalArgs.generateReportTable(covariates.covariateNames()), quantizationInfo.generateReportTable(), RecalUtils.generateReportTables(combinedTables, covariates));
    }
}
