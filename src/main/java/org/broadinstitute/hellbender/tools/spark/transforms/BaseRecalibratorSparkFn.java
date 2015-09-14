package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorFn;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.walkers.bqsr.ReadRecalibrationInfo;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Iterator;

public class BaseRecalibratorSparkFn {

    public static RecalibrationTables apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary ) {
        final BaseRecalibrationArgumentCollection brac = new BaseRecalibrationArgumentCollection();
        final BQSRSpark bqsr = new BQSRSpark(header, referenceDictionary, brac);
        bqsr.onTraversalStart();
        JavaRDD<RecalibrationTables> unmergedTables = readsWithContext.mapPartitions(new FlatMapFunction<Iterator<Tuple2<GATKRead, ReadContextData>>, RecalibrationTables>() {
            @Override
            public Iterable<RecalibrationTables> call( Iterator<Tuple2<GATKRead, ReadContextData>> readWithContextIterator ) throws Exception {
                while ( readWithContextIterator.hasNext() ) {
                    final Tuple2<GATKRead, ReadContextData> readWithData = readWithContextIterator.next();
                    Iterable<Variant> variants = readWithData._2().getOverlappingVariants();
                    final ReferenceBases refBases = readWithData._2().getOverlappingReferenceBases();
                    ReferenceDataSource refDS = new ReferenceMemorySource(refBases, header.getSequenceDictionary());

                    bqsr.apply(readWithData._1(), refDS, variants);
                }
                return Arrays.asList(bqsr.getRecalibrationTable());
            }
        });


    }

}
