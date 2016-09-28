package org.broadinstitute.hellbender.tools.spark.sv;


import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.io.SyncFailedException;
import java.util.Map;

@CommandLineProgramProperties(summary="Compute metrics on the number of chimeric read pairs stratified by MAPQ",
        oneLineSummary="Compute metrics on the number of chimeric read pairs stratified by MAPQ",
        programGroup = SparkProgramGroup.class)
public class CalculateChimeraCountsByMapqSpark extends GATKSparkTool {
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

        final JavaPairRDD<String, ChimericPairInfo> stringChimericPairInfoJavaPairRDD =
                getReads().filter(r -> !r.isSecondaryAlignment() && !r.isSupplementaryAlignment() && !r.failsVendorQualityCheck())
                .filter(read -> isPartOfChimera(read) > 0)
                .mapToPair(read -> new Tuple2<>(read.getName(), read))
                .aggregateByKey(new ChimericPairInfo(),
                        CalculateChimeraCountsByMapqSpark::combine,
                        (ChimericPairInfo v1, ChimericPairInfo v2) -> {
                            if (! v1.initialized) {
                                return v2;
                            }
                            if (! v2.initialized) {
                                return v1;
                            }
                            ChimericPairInfo res = new ChimericPairInfo();
                            res.mapq1 = Math.min(v1.mapq1, v2.mapq1);
                            res.mapq2 = Math.min(v1.mapq2, v2.mapq2);
                            if (v1.chimeraType != v2.chimeraType) {
                                System.err.println("Got conflicting chimera types: " + v1.chimeraType + "-" + v2.chimeraType);
                            }
                            res.chimeraType = Math.min(v1.chimeraType, v2.chimeraType);
                            return res;
                        });
        final Map<Tuple2<Integer, Integer>, Long> chimericPairInfoLongMap =
                stringChimericPairInfoJavaPairRDD.map(v -> new Tuple2<>(Math.min(v._2().mapq1, v._2().mapq2), v._2().chimeraType)).countByValue();
        for (Tuple2<Integer, Integer> key : chimericPairInfoLongMap.keySet()) {
            System.out.println(key._1() + "\t" + key._2() + "\t" + chimericPairInfoLongMap.get(key));
        }

    }

    public static ChimericPairInfo combine(final ChimericPairInfo v1, final GATKRead v2) {

        if (v1.initialized) {
            v1.chimeraType = Math.max(v1.chimeraType, isPartOfChimera(v2));
            if (v2.isFirstOfPair()) {
                v1.mapq1 = Math.min(v1.mapq1, v2.getMappingQuality());
            } else {
                v1.mapq2 = Math.min(v1.mapq2, v2.getMappingQuality());
            }
            return v1;
        } else {
            v1.initialized = true;
            v1.chimeraType = isPartOfChimera(v2);
            if (v2.isFirstOfPair()) {
                v1.mapq1 = v2.getMappingQuality();
            } else {
                v1.mapq2 = v2.getMappingQuality();
            }
            return v1;

        }

    }

    static class ChimericPairInfo implements Serializable{
        private static final long serialVersionUID = 1L;

        boolean initialized = false;
        int mapq1 = Integer.MAX_VALUE;
        int mapq2 = Integer.MAX_VALUE;
        int chimeraType;

        @Override
        public String toString() {
            return "ChimericPairInfo{" +
                    "mapq1=" + mapq1 +
                    ", mapq2=" + mapq2 +
                    ", chimeraType=" + chimeraType +
                    '}';
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final ChimericPairInfo that = (ChimericPairInfo) o;

            if (mapq1 != that.mapq1) return false;
            if (mapq2 != that.mapq2) return false;
            return chimeraType == that.chimeraType;

        }

        @Override
        public int hashCode() {
            int result = mapq1;
            result = 31 * result + mapq2;
            result = 31 * result + chimeraType;
            return result;
        }
    }

    public static int isPartOfChimera(final GATKRead read) {
        if (read.isUnmapped()) {
            return 0;
        }
        if (read.isPaired() && !read.mateIsUnmapped()) {
            if (!read.getContig().equals(read.getMateContig())) {
                return 1;
            }
            if (read.isReverseStrand() && read.mateIsReverseStrand() || ! read.isReverseStrand() && ! read.mateIsReverseStrand()) {
                return 2;
            }
            if (((read.getStart() < read.getMateStart() && read.isReverseStrand()) ||
                            (read.getStart() > read.getMateStart() && !read.isReverseStrand())))  {
                return 3;
            }
            if (Math.abs(read.getFragmentLength()) > 100000) {
                return 4;
            }
            return 0;
        } else {
            if (read.hasAttribute("SA")) {
                return 5;
            }
        }
        return 0;
    }
}
