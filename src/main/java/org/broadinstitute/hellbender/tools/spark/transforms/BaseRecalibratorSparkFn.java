package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
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
}
