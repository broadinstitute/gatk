package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorFn;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.walkers.bqsr.RecalibrationEngine;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class BaseRecalibratorSparkFn {

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary, final BaseRecalibrationArgumentCollection args ) {
        final BQSRSparkWorker bqsr = new BQSRSparkWorker(header, referenceDictionary, args);
        bqsr.onTraversalStart();

        final ReadFilter filter = readFilter(header);
        readsWithContext.filter(readWithContext -> filter.apply(readWithContext._1()));

        JavaRDD<RecalibrationTables> unmergedTables = readsWithContext.mapPartitions(readWithContextIterator -> {
            while ( readWithContextIterator.hasNext() ) {
                final Tuple2<GATKRead, ReadContextData> readWithData = readWithContextIterator.next();
                Iterable<Variant> variants = readWithData._2().getOverlappingVariants();
                final ReferenceBases refBases = readWithData._2().getOverlappingReferenceBases();
                ReferenceDataSource refDS = new ReferenceMemorySource(refBases, referenceDictionary);

                bqsr.apply(readWithData._1(), refDS, variants);
            }
            return new ArrayList<>(Arrays.asList(bqsr.getRecalibrationTable()));
        });

        final RecalibrationTables combinedTables = unmergedTables.reduce(RecalibrationTables::safeCombine);
        RecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        File temp = IOUtils.createTempFile("temp-recalibrationtable-", ".tmp");
        try {
            BaseRecalibratorFn.SaveTextualReport(temp, header, combinedTables, args);
            return new RecalibrationReport(temp);
        } catch (FileNotFoundException e) {
            throw new GATKException("can't find my own temporary file", e);
        } catch (IOException e) {
            throw new GATKException("unable to save temporary report to " + temp.getPath(), e);
        }
    }


    private static CountingReadFilter readFilter( final SAMFileHeader header ) {
        return new CountingReadFilter("Wellformed", new WellformedReadFilter(header))
                .and(new CountingReadFilter("Mapping_Quality_Not_Zero", ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO))
                .and(new CountingReadFilter("Mapping_Quality_Available", ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE))
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED))
                .and(new CountingReadFilter("Primary_Alignment", ReadFilterLibrary.PRIMARY_ALIGNMENT))
                .and(new CountingReadFilter("Not_Duplicate", ReadFilterLibrary.NOT_DUPLICATE))
                .and(new CountingReadFilter("Passes_Vendor_Quality_Check", ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK));
    }
}
