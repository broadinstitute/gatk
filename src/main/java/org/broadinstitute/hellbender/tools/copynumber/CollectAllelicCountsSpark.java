package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerContext;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerSpark;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollector;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleNameUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.util.Collections;
import java.util.Iterator;

public class CollectAllelicCountsSpark extends LocusWalkerSpark {


    @Override
    protected void processAlignments(JavaRDD<LocusWalkerContext> rdd, JavaSparkContext ctx) {
        final String sampleName = SampleNameUtils.readSampleName(getHeaderForReads());
        final SampleMetadata sampleMetadata = new SimpleSampleMetadata(sampleName);
        final Broadcast<SampleMetadata> sampleMetadataBroadcast = ctx.broadcast(sampleMetadata);
        // rdd.map(pileupFunction(metadata, outputInsertLength, showVerbose)).saveAsTextFile(outputFile);
        rdd.mapPartitions(distributedCount(sampleMetadataBroadcast.getValue(), 20)).reduce();
    }

    private static FlatMapFunction<Iterator<LocusWalkerContext>, AllelicCountCollector> distributedCount(final SampleMetadata sampleMetadata,
                                                                                                         final int minimumBaseQuality) {
        return (FlatMapFunction<Iterator<LocusWalkerContext>, AllelicCountCollector>) contextIterator -> {
            final AllelicCountCollector result = new AllelicCountCollector(sampleMetadata);

            contextIterator.forEachRemaining( ctx -> {
                final ReferenceContext referenceContext = ctx.getReferenceContext();
                final byte refAsByte = referenceContext.getBase();
                result.collectAtLocus(Nucleotide.valueOf(refAsByte), ctx.getAlignmentContext().getBasePileup(),
                        ctx.getAlignmentContext().getLocation(), minimumBaseQuality);
                }
            );
            return Collections.singletonList(result).iterator();
        };
    }

    private static AllelicCountCollector combineAllelicCountCollectors(final AllelicCountCollector allelicCountCollector1, final AllelicCountCollector allelicCountCollector2) {

    }
}
