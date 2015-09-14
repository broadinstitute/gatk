package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
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

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary ) {
        final BaseRecalibrationArgumentCollection brac = new BaseRecalibrationArgumentCollection();
        final BQSRSpark bqsr = new BQSRSpark(header, referenceDictionary, brac);
        bqsr.onTraversalStart();
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
            BaseRecalibratorFn.SaveTextualReport(temp, header, combinedTables, brac);
            return new RecalibrationReport(temp);
        } catch (FileNotFoundException e) {
            throw new GATKException("can't find my own temporary file", e);
        } catch (IOException e) {
            throw new GATKException("unable to save temporary report to " + temp.getPath(), e);
        }
    }

}
