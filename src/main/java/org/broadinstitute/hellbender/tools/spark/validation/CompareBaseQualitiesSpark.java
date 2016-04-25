package org.broadinstitute.hellbender.tools.spark.validation;

import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.TestSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.validation.CompareMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.transformers.BQSRReadTransformer.constructStaticQuantizedMapping;

@CommandLineProgramProperties(summary = "Compares two the quality scores of BAMs", oneLineSummary = "Diff qs of the BAMs",
        programGroup = TestSparkProgramGroup.class)

final public class CompareBaseQualitiesSpark extends GATKSparkTool  {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the second BAM", shortName = "I2", fullName = "input2", optional = false)
    protected String input2;

    @Argument(doc="summary output file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = true)
    protected String outputFilename = null;

    @Argument(doc="throw error on diff", shortName = "cd", fullName = "throwOnDiff", optional = true)
    protected boolean throwOnDiff = false;

    /**
     * Use static quantized quality scores to a given number of levels.
     */
    @Advanced
    @Argument(fullName="static_quantized_quals", shortName = "SQQ", doc = "Use static quantized quality scores to a given number of levels (with -"+ StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME+ ")", optional=true)
    public List<Integer> staticQuantizationQuals = new ArrayList<>();

    /**
     * Round down quantized only works with the static_quantized_quals option.  When roundDown = false, rounding is done in
     * probability space to the nearest bin.  When roundDown = true, the value is rounded to the nearest bin
     * that is smaller than the current bin.
     */
    @Advanced
    @Argument(fullName="round_down_quantized", shortName = "RDQ", doc = "Round quals down to nearest quantized qual", optional=true)
    public boolean roundDown = false;


    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final byte[] staticQuantizedMapping = constructStaticQuantizedMapping(staticQuantizationQuals, roundDown);

        JavaRDD<GATKRead> firstReads = getReads();
        ReadsSparkSource readsSource2 = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());
        JavaRDD<GATKRead> secondReads = readsSource2.getParallelReads(input2, referenceArguments.getReferenceFileName(), getIntervals(), bamPartitionSplitSize);

        long firstBamSize = firstReads.count();
        long secondBamSize = secondReads.count();
        // If a user gives two BAMs that aren't the same, we throw because this is invalid use of the tool.
        if (firstBamSize != secondBamSize) {
            throw new UserException("input bams have different numbers of total reads: "
                    + firstBamSize + "," + secondBamSize);
        }
        logger.info("Processing bams with " + firstBamSize + " reads");

        // To compare bases, we created two PairRDDs keyed by read name, strand, etc. where the value is the bases.
        // We cogroup the two RDDs, verify there is only one read per key per bam, and then diff the bases
        // position by position for the two reads and put the results in CompareMatrix.
        JavaPairRDD<String, Quals> firstQuals = firstReads.mapToPair(read -> {
            // For both BAMs, we clear the attributes because want to key by common unique per-read properties, start, contig, etc.
            return new Tuple2<>(getReadKey(read), new Quals(read.getBaseQualities()));
        });
        JavaPairRDD<String, Quals> secondQuals = secondReads.mapToPair(read -> {
            return new Tuple2<>(getReadKey(read), new Quals(read.getBaseQualities()));
        });

        JavaPairRDD<String, Tuple2<Iterable<Quals>, Iterable<Quals>>> cogroup = firstQuals.cogroup(secondQuals, getRecommendedNumReducers());
        final CompareMatrix finalMatrix = cogroup.map(v1 -> {
            List<Quals> lFirstQuals = Lists.newArrayList(v1._2()._1());
            List<Quals> lSecondQuals = Lists.newArrayList(v1._2()._2());
            CompareMatrix compareMatrix = new CompareMatrix(staticQuantizedMapping);
            if (lFirstQuals.size() != 1) {
                throw new GATKException("expected only one read per key in first bam: " + lFirstQuals.size());
            }
            if (lSecondQuals.size() != 1) {
                throw new GATKException("expected only one read per key in second bam: " + lSecondQuals.size());
            }

            compareMatrix.add(lFirstQuals.get(0).quals, lSecondQuals.get(0).quals);
            return compareMatrix;
        }).treeAggregate(new CompareMatrix(staticQuantizedMapping),
                CompareMatrix::add,
                CompareMatrix::add);


        finalMatrix.printOutResults(outputFilename);

        if (throwOnDiff && finalMatrix.hasNonDiagonalElements()) {
            throw new UserException("Quality scores from the two BAMs do not match");
        }
    }

    static String getReadKey(GATKRead r) {
        if (r.isUnmapped() || r.getCigar().isEmpty()) {
            return String.format("%s UNMAPPED", r.getName());
        } else {
            return String.format("%s %s:%d-%d", r.getName(), r.getContig(), r.getStart(), r.getEnd());

        }
    }

    private static class Quals implements Serializable {
        private static final long serialVersionUID = 1L;

        public byte[] quals;

        public Quals(byte[] quals) {
            this.quals = quals;
        }
    }

}
