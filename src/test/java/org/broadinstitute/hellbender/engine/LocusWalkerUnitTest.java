package org.broadinstitute.hellbender.engine;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
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

            // transforming name of unmapped mates
            return new ReadTransformer() {
                private static final long serialVersionUID = 1L;
                @Override
                public GATKRead apply(GATKRead gatkRead) {
                    if (gatkRead.isReverseStrand()) {
                        preTransformed++;
                    }
                    return gatkRead;
                }
            };
        }
        @Override
        public ReadTransformer makePostReadFilterTransformer() {

            // transforming name of unmapped mates
            return new ReadTransformer() {
                private static final long serialVersionUID = 1L;
                @Override
                public GATKRead apply(GATKRead gatkRead) {
                    if (gatkRead.mateIsUnmapped()) {
                        postTransformed++;
                    }
                    return gatkRead;
                }
            };
        }
        public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
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

}