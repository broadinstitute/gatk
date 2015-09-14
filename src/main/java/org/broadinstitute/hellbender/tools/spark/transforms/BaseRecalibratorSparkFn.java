package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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
            return new ArrayList<>(Arrays.asList(bqsr.getRecalibrationTables()));
        });

        final RecalibrationTables combinedTables = unmergedTables.reduce(RecalibrationTables::safeCombine);
        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        File temp = IOUtils.createTempFile("temp-recalibrationtable-", ".tmp");
        try {
            BaseRecalibratorFn.saveTextualReport(temp, header, combinedTables, recalArgs);
            return new RecalibrationReport(temp);
        } catch (FileNotFoundException e) {
            throw new GATKException("can't find my own temporary file", e);
        } catch (IOException e) {
            throw new GATKException("unable to save temporary report to " + temp.getPath(), e);
        }
    }
}
