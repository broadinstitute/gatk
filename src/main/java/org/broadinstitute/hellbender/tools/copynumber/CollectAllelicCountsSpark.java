package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerContext;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerSpark;
import org.broadinstitute.hellbender.tools.copynumber.datacollection.AllelicCountCollector;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.*;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * See {@link CollectAllelicCounts}.  This behaves the same, except that it supports spark.
 */
@CommandLineProgramProperties(
        summary = "Collects ref/alt counts at sites.",
        oneLineSummary = "Collects ref/alt counts at sites.",
        programGroup = CopyNumberProgramGroup.class
)
public class CollectAllelicCountsSpark extends LocusWalkerSpark {

    private static final long serialVersionUID = 1L;

    private static final Logger logger = LogManager.getLogger(CollectAllelicCounts.class);

    @Argument(
            doc = "Output allelic-counts file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputAllelicCountsFile;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            fullName = "minimumBaseQuality",
            shortName = "minBQ",
            minValue = 0,
            optional = true
    )
    private int minimumBaseQuality = 20;

    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

    @Override
    public boolean emitEmptyLoci() {return true;}

    @Override
    public boolean requiresReference() {return true;}

    @Override
    public boolean requiresIntervals() {return true;}

    @Override
    protected void processAlignments(JavaRDD<LocusWalkerContext> rdd, JavaSparkContext ctx) {
        final SampleLocatableMetadata metadata = MetadataUtils.fromHeader(getHeaderForReads(), Metadata.Type.SAMPLE_LOCATABLE);
        final Broadcast<SampleLocatableMetadata> sampleMetadataBroadcast = ctx.broadcast(metadata);

        final AllelicCountCollector finalAllelicCountCollector =
                rdd.mapPartitions(distributedCount(sampleMetadataBroadcast, minimumBaseQuality))
                .reduce((a1, a2) -> combineAllelicCountCollectors(a1, a2, sampleMetadataBroadcast.getValue()));
        finalAllelicCountCollector.getAllelicCounts().write(outputAllelicCountsFile);
    }

    private static FlatMapFunction<Iterator<LocusWalkerContext>, AllelicCountCollector> distributedCount(final Broadcast<SampleLocatableMetadata> sampleMetadataBroadcast,
                                                                                                         final int minimumBaseQuality) {
        return (FlatMapFunction<Iterator<LocusWalkerContext>, AllelicCountCollector>) contextIterator -> {
            final AllelicCountCollector result = new AllelicCountCollector(sampleMetadataBroadcast.getValue());

            contextIterator.forEachRemaining( ctx -> {
                final byte refAsByte = ctx.getReferenceContext().getBase();
                result.collectAtLocus(Nucleotide.decode(refAsByte), ctx.getAlignmentContext().getBasePileup(),
                        ctx.getAlignmentContext().getLocation(), minimumBaseQuality);
                }
            );
            return Collections.singletonList(result).iterator();
        };
    }

    private static AllelicCountCollector combineAllelicCountCollectors(final AllelicCountCollector allelicCountCollector1,
                                                                       final AllelicCountCollector allelicCountCollector2,
                                                                       final SampleLocatableMetadata sampleMetadata) {
        return AllelicCountCollector.combine(allelicCountCollector1, allelicCountCollector2, sampleMetadata);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> initialReadFilters = new ArrayList<>(super.getDefaultReadFilters());
        initialReadFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));

        return initialReadFilters;
    }
}
