package org.broadinstitute.hellbender.engine;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LocusWalkerUnitTest extends CommandLineProgramTest {

    @CommandLineProgramProperties(
            summary = "Dummy that reads file and counts how many pileup elements are transformed",
            oneLineSummary = "none",
            programGroup = TestProgramGroup.class
    )
    private static class TestTransformedLocusWalker extends LocusWalker {
        public int preTransformed = 0;
        public int postTransformed = 0;
        public int totalPileup = 0;

        @Override
        public ReadTransformer makePreReadFilterTransformer() {
            return new ReadTransformer() {
                private static final long serialVersionUID = 1L;
                @Override
                public GATKRead apply(GATKRead gatkRead) {
                    // test that the attribute is not set
                    Assert.assertFalse(gatkRead.hasAttribute("tr"));
                    gatkRead.setAttribute("tr", "PRE");

                    if (gatkRead.isReverseStrand()) {
                        preTransformed++;
                    }
                    return gatkRead;
                }
            };
        }

        @Override
        public CountingReadFilter makeReadFilter() {
            return new CountingReadFilter(new ReadFilter() {
                private static final long serialVersionUID = 1L;
                @Override
                public boolean test(GATKRead read) {
                    // test that the attribute is pre-transformed
                    Assert.assertEquals(read.getAttributeAsString("tr"), "PRE");
                    // test that the attribute is not set
                    Assert.assertFalse(read.hasAttribute("ft"));
                    // and set as filtered
                    read.setAttribute("ft", "yes");
                    return true;
                }
            });
        }

        @Override
        public ReadTransformer makePostReadFilterTransformer() {
            return new ReadTransformer() {
                private static final long serialVersionUID = 1L;
                @Override
                public GATKRead apply(GATKRead gatkRead) {
                    Assert.assertEquals(gatkRead.getAttributeAsString("tr"), "PRE");
                    gatkRead.setAttribute("tr", "POST");

                    if (gatkRead.mateIsUnmapped()) {
                        postTransformed++;
                    }
                    return gatkRead;
                }
            };
        }

        public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            alignmentContext.getBasePileup().getReads().forEach(read -> Assert.assertEquals(read.getAttributeAsString("tr"), "POST"));
            totalPileup++;
        }
    }

    @Test
    public void testTransformReads() throws IOException {

        final TestTransformedLocusWalker tool = new TestTransformedLocusWalker();

        final String[] args = {
                "-I", getTestDataDir()+ "/print_reads.sorted.bam",
                "-R", getTestDataDir()+ "/print_reads.fasta",
                "-L", "chr7:21-21"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.preTransformed, 4);
        Assert.assertEquals(tool.postTransformed, 1);
        Assert.assertEquals(tool.totalPileup, 1);
    }

    private static class TestEmitUncoveredLociTool extends LocusWalker {
        public int totalApplyCalls = 0;

        @Override
        public boolean emitEmptyLoci() {
            return true;
        }

        @Override
        public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            totalApplyCalls++;
        }
    }

    @Test
    public void testEmitUncoveredLoci() {
        // Most testing of this is in IntervalAlignmentContextIteratorUnitTest
        final TestEmitUncoveredLociTool tool = new TestEmitUncoveredLociTool();

        final String[] args = {
                "-I", getTestDataDir()+ "/print_reads.sorted.bam",
                "-R", getTestDataDir()+ "/print_reads.fasta",
                "-L", "chr7:21-30"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.totalApplyCalls, 10);
    }

}