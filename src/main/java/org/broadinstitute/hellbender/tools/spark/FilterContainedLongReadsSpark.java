package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;

@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Retain only alignment records that are not \"contained\" in the tart contig they are mapped to",
        summary =
                "Retain only alignment records that are not \"contained\" in the tart contig they are mapped to",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class FilterContainedLongReadsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    private static final int REF_CONTIG_ROOM = 10;

    @Argument(doc = "output BAM name",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outName;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.MAPPED);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final SAMSequenceDictionary referenceSequenceDictionary = getReferenceSequenceDictionary();
        JavaRDD<GATKRead> notContainedReads = reads.filter(read -> {
            final String contig = read.getContig();
            if (null != contig) {
                int end = read.getEnd();
                int sequenceLength = referenceSequenceDictionary.getSequence(contig).getSequenceLength();

                int start = read.getStart();
                return start < REF_CONTIG_ROOM ||
                        sequenceLength - end < REF_CONTIG_ROOM ; // leave some room?
            } else {
                return false;
            }
        });

        writeReads(ctx, outName, notContainedReads);
    }
}
