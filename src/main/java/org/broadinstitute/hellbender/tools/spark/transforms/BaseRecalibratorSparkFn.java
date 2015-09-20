package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
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
import java.util.Map;

public class BaseRecalibratorSparkFn {

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, Iterable<Variant>> readsWithVariants, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary, final RecalibrationArgumentCollection recalArgs, final Map<String, ReferenceBases> allReferenceBases, final JavaSparkContext ctx) {

        final Broadcast<Map<String, ReferenceBases>> refBroadcast = ctx.broadcast(allReferenceBases);

        final ReadFilter filter = BaseRecalibrator.getStandardBQSRReadFilter(header);
        final JavaPairRDD<GATKRead, Iterable<Variant>> filtered = readsWithVariants.filter(readWithVariants -> filter.apply(readWithVariants._1()));

        JavaRDD<RecalibrationTables> unmergedTables = filtered.mapPartitions(readWithVariants -> {
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);

            while ( readWithVariants.hasNext() ) {
                final Tuple2<GATKRead, Iterable<Variant>> readWithData = readWithVariants.next();
                GATKRead read = readWithData._1();
                Iterable<Variant> variants = readWithData._2();
                final ReferenceBases referenceBases = refBroadcast.getValue().get(read.getContig());
                ReferenceDataSource refDS = new ReferenceMemorySource(referenceBases, referenceDictionary);
                bqsr.processRead(read, refDS, variants);
            }
            return new ArrayList<>(Arrays.asList(bqsr.getRecalibrationTables()));
        });

        final RecalibrationTables combinedTables = unmergedTables.reduce(RecalibrationTables::safeCombine);
        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        File temp = IOUtils.createTempFile("temp-recalibrationtable-", ".tmp");
        try {
            BaseRecalibratorFn.SaveTextualReport(temp, header, combinedTables, recalArgs);
            return new RecalibrationReport(temp);
        } catch (FileNotFoundException e) {
            throw new GATKException("can't find my own temporary file", e);
        } catch (IOException e) {
            throw new GATKException("unable to save temporary report to " + temp.getPath(), e);
        }
    }
}
