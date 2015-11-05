package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.Lists;
import org.apache.commons.collections4.iterators.IteratorIterable;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(summary = "Compares two BAMs", oneLineSummary = "Compares two BAMs",
        programGroup = SparkProgramGroup.class)
public class BqsrRecalBasesCompare extends GATKSparkTool  {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean requiresReads() { return true; }

        @Argument(doc="the second BAM", shortName = "I2", fullName = "input2", optional = false)
        protected String input2;

        @Override
        protected void runTool(final JavaSparkContext ctx) {
            JavaRDD<GATKRead> firstReads = getReads();
            ReadsSparkSource readsSource2 = new ReadsSparkSource(ctx);
            JavaRDD<GATKRead> secondReads = readsSource2.getParallelReads(input2, getIntervals());

            long firstBamSize = firstReads.count();
            long secondBamSize = secondReads.count();
            if (firstBamSize != secondBamSize) {
                throw new UserException("input bams have different numbers of total reads: "
                        + firstBamSize + "," + secondBamSize);
            }
            System.out.println("processing bams with " + firstBamSize + " reads");

            JavaPairRDD<String, Quals> firstQuals = firstReads.mapToPair(read -> {
                read.clearAttributes();
                return new Tuple2<>(read.toString(), new Quals(read.getBaseQualities()));
            });
            JavaPairRDD<String, Quals> secondQuals = secondReads.mapToPair(read -> {
                read.clearAttributes();
                return new Tuple2<>(read.toString(), new Quals(read.getBaseQualities()));
            });

            JavaPairRDD<String, Tuple2<Iterable<Quals>, Iterable<Quals>>> cogroup = firstQuals.cogroup(secondQuals);
            JavaRDD<Tuple2<Tuple2<Integer, Integer>, Integer>> map = cogroup.map(v1 -> {
                List<Quals> firstQuals1 = Lists.newArrayList(v1._2()._1());
                List<Quals> secondQuals1 = Lists.newArrayList(v1._2()._2());
                if (firstQuals1.size() != 1) {
                    throw new GATKException("expected only one read (first): " + firstQuals1.size());
                }
                if (secondQuals1.size() != 1) {
                    throw new GATKException("expected only one read (second): " + secondQuals1.size());
                }
                if (qFirst.size!= qSecond.size())
                for(int i = 0; i < qFirst.size()) {

                }
                for(byte b : firstQuals1.get(0).quals) {
                    ((int) b);
                }
                return null;
            });
        }
}

class Quals implements Serializable {
    private static final long serialVersionUID = 1L;

    public byte[] quals;

    public Quals(byte[] quals) {
        this.quals = quals;
    }
}
