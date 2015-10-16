package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
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
import java.util.List;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;


public class BaseRecalibratorSparkFn {

    public static Logger log = LogManager.getLogger(BaseRecalibratorSparkFn.class);

    static class ComparisonResult {
        long matching = 0;
        long deviate = 0;
        long bigDeviate = 0;
        // default: 1% errors are accounted separately.
        double bigThreshold = 0.01;

        public void addData(double expected, double got) {
            if (expected==got) {
                matching++;
            } else {
                deviate++;
                double diffRatio = Math.abs( (expected-got) / ((expected+got)/2.0));
                if (diffRatio > bigThreshold) bigDeviate++;
            }
        }

        public String formatResult() {
            long total = matching+deviate;
            double matchingPct = 100.0 * (matching / (double)total);
            return "Matching: "+matching+"/"+total+" ("+matchingPct+"%). Errors: "+deviate+" Big errors (>"+(100*bigThreshold)+"%): "+bigDeviate;
        }
    }

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary, final RecalibrationArgumentCollection recalArgs ) {
        final ReadFilter filter = BaseRecalibrator.getStandardBQSRReadFilter(header);
        final JavaPairRDD<GATKRead, ReadContextData> filtered = readsWithContext.filter(readWithContext -> filter.apply(readWithContext._1()));

        JavaRDD<RecalibrationTables> unmergedTables = filtered.mapPartitions(readWithContextIterator -> {
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);
            bqsr.logCovariatesUsed();

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
            RecalibrationReport ret = new RecalibrationReport(tempReport);

            checkForRoundtripDifference(combinedTables, ret.getRecalibrationTables(), "");

            RecalibrationTables victim = ret.getRecalibrationTables();
            roundToTwoDecimalPlaces(victim);
            checkForRoundtripDifference(combinedTables, victim, "(after rounding) ");

            return ret;
        } catch (FileNotFoundException e) {
            throw new GATKException("can't find my own temporary file", e);
        }
    }

    // if report save/reload is a no-op then this should output "100% matching".
    private static void checkForRoundtripDifference(RecalibrationTables original, RecalibrationTables ret, String msg) {

        ComparisonResult[] comp = new ComparisonResult[3];
        for (int i=0; i<comp.length; i++) comp[i] = new ComparisonResult();

        for (int i=0; i<original.numTables(); i++) {
            List<NestedIntegerArray.Leaf<RecalDatum>> l1 = original.getTable(i).getAllLeaves();
            List<NestedIntegerArray.Leaf<RecalDatum>> l2 = ret.getTable(i).getAllLeaves();
            for (int j = 0; j<l1.size(); j++) {
                RecalDatum a = l1.get(j).value;
                RecalDatum b = l2.get(j).value;
                comp[0].addData(a.getNumMismatches(), b.getNumMismatches());
                comp[1].addData(a.getNumObservations(), b.getNumObservations());
                comp[2].addData(a.getEmpiricalQuality(), b.getEmpiricalQuality());
            }
        }

        log.info("RecalibrationTables result. "+ msg + "NumMismatches: " + comp[0].formatResult());
        log.info("RecalibrationTables result. "+ msg + "NumObservations: " + comp[1].formatResult());
        log.info("RecalibrationTables result. "+ msg + "Empirical quality: " + comp[2].formatResult());

    }

    private static void roundToTwoDecimalPlaces(RecalibrationTables rt) {
        for (int i=0; i<rt.numTables(); i++) {
            for (NestedIntegerArray.Leaf<RecalDatum> l : rt.getTable(i).getAllLeaves()) {
                l.value.setNumMismatches( Math.round(l.value.getNumMismatches()*100)/100.0 );
            }
        }
    }
}
