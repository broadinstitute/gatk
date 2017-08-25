package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.barclay.help.DocumentedFeature;

@CommandLineProgramProperties(summary="Find count of each contig in a BAM file using Spark",
                oneLineSummary="Find count of each contig in a BAM file using Spark",
                programGroup = SparkProgramGroup.class)
@DocumentedFeature
public class ViewSequenceDictionary extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        System.out.println("Ref");
        SAMSequenceDictionary referenceSequenceDictionary = getReferenceSequenceDictionary();
        referenceSequenceDictionary.getSequences().forEach(System.out::println);

        System.out.println("Reads");
        SAMSequenceDictionary readsSequenceDictionary = getHeaderForReads().getSequenceDictionary();
        readsSequenceDictionary.getSequences().forEach(System.out::println);

        referenceSequenceDictionary.assertSameDictionary(readsSequenceDictionary);
    }
}
