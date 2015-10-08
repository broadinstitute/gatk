package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

public class BaseRecalibratorSparkFn {

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary, final RecalibrationArgumentCollection recalArgs ) {
        final ReadFilter filter = BaseRecalibrator.getStandardBQSRReadFilter(header);
        final JavaPairRDD<GATKRead, ReadContextData> filtered = readsWithContext.filter(readWithContext -> filter.apply(readWithContext._1()));

        JavaRDD<RecalibrationTables> unmergedTables = filtered.mapPartitions(readWithContextIterator -> {
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);

            while ( readWithContextIterator.hasNext() ) {
                final Tuple2<GATKRead, ReadContextData> readWithData = readWithContextIterator.next();
                Iterable<Variant> variants = readWithData._2().getOverlappingVariants();
                final ReferenceBases refBases = readWithData._2().getOverlappingReferenceBases();
                ReferenceDataSource refDS = new ReferenceMemorySource(refBases, referenceDictionary);

                bqsr.processRead(readWithData._1(), refDS, variants);
            }
            // Need to wrap in ArrayList due to our current inability to serialize the return value of Arrays.asList() directly
            return new ArrayList<>(Arrays.asList(bqsr.getRecalibrationTables()));
        });

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(recalArgs, header));
        final RecalibrationTables combinedTables = unmergedTables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(unmergedTables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        try {
            final QuantizationInfo quantizationInfo = new QuantizationInfo(combinedTables, recalArgs.QUANTIZING_LEVELS);

            // TODO: we write the report out and read it back in only because the report changes when written
            // TODO: out and read back in -- remove this step once we fix this issue.
            File tempReport = IOUtils.createTempFile("temp-recalibrationtable-", ".tmp");
            try ( PrintStream reportStream = new PrintStream(tempReport) ) {
                RecalUtils.outputRecalibrationReport(reportStream, recalArgs, quantizationInfo, combinedTables, new StandardCovariateList(recalArgs, header), recalArgs.SORT_BY_ALL_COLUMNS);
            }
            return new RecalibrationReport(tempReport);
        } catch (FileNotFoundException e) {
            throw new GATKException("can't find my own temporary file", e);
        }
    }
}
