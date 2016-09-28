package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMTag;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.Map;

@CommandLineProgramProperties(summary="Compute pct chimeras ala picard",
        oneLineSummary="Compute pct chimeras ala picard",
        programGroup = SparkProgramGroup.class)
public class CalculatePctChimerasSpark extends GATKSparkTool {
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
    protected void runTool(final JavaSparkContext ctx) {
        final Map<Integer, Long> chimeraCountsByType =
                getUnfilteredReads().
                        filter(r -> !r.isSecondaryAlignment() &&  !r.isSupplementaryAlignment() && ! r.failsVendorQualityCheck() && ! r.isUnmapped()).
                        map(r -> isPartOfChimera(r)).
                        countByValue();
        for (Integer key : chimeraCountsByType.keySet()) {
            System.out.println(key + "\t" + chimeraCountsByType.get(key));
        }

        long denom = chimeraCountsByType.get(0);
        long totalChimeras = 0;
        for (int i = 1; i <= 5; i++) {
            totalChimeras = totalChimeras + chimeraCountsByType.get(i);
            denom = denom + chimeraCountsByType.get(i);
        }

        System.out.println("denom: " + denom);
        System.out.println("chimeras: " + totalChimeras);
        System.out.println("pct: " + ((double) totalChimeras) / denom);
    }

    public static int isPartOfChimera(final GATKRead read) {
        if (read.isUnmapped()) {
            return -1;
        }
        if (read.isPaired() && !read.mateIsUnmapped()) {
            if (Math.abs(read.getFragmentLength()) > 100000) {
                return 4;
            }
            if (!read.getContig().equals(read.getMateContig())) {
                return 1;
            }
            if (read.hasAttribute("SA")) {
                return 5;
            }
            if (read.isReverseStrand() && read.mateIsReverseStrand() || !read.isReverseStrand() && !read.mateIsReverseStrand()) {
                return 2;
            }
            if (((read.getStart() < read.getMateStart() && read.isReverseStrand()) ||
                    (read.getStart() > read.getMateStart() && !read.isReverseStrand()))) {
                return 3;
            }
            return 0;
        } else {
            if (read.hasAttribute("SA")) {
                return 5;
            }
            return 0;
        }
    }

}