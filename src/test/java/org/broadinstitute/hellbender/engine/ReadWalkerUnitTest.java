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
public class ReadWalkerUnitTest extends CommandLineProgramTest {

    @CommandLineProgramProperties(
            summary = "Dummy that reads file and count how many reads are transformed",
            oneLineSummary = "empty class",
            programGroup = TestProgramGroup.class
    )
    private static class TestTransformedReadWalker extends ReadWalker {
        public int preTransformed = 0;
        public int postTransformed = 0;
        public int totalReads = 0;

        @Override
        public ReadTransformer makePreReadFilterTransformer() {
            return new ReadTransformer() {
                private static final long serialVersionUID = 1L;
                @Override
                public GATKRead apply(GATKRead gatkRead) {
                    // test that the attribute is not set
                    Assert.assertFalse(gatkRead.hasAttribute("tr"));
                    // and set as pre-transformed
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
                    // test that the attribute is filtered
                    Assert.assertEquals(gatkRead.getAttributeAsString("ft"), "yes");
                    // test that the attribute is pre-transformed
                    Assert.assertEquals(gatkRead.getAttributeAsString("tr"), "PRE");
                    // and set as post-transformed
                    gatkRead.setAttribute("tr", "POST");

                    if (gatkRead.mateIsUnmapped()) {
                        postTransformed++;
                    }
                    return gatkRead;
                }
            };
        }


        @Override
        public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
            Assert.assertEquals(read.getAttributeAsString("tr"), "POST");
            totalReads++;
        }
    }

    @Test
    public void testTransformReads() throws IOException {

        final TestTransformedReadWalker tool = new TestTransformedReadWalker();

        final String[] args = {
                "-I", getTestDataDir()+ "/print_reads.sorted.bam",
                "-R", getTestDataDir()+ "/print_reads.fasta",
                "-L", "chr7:21-21"
        };

        tool.instanceMain(args);

        Assert.assertEquals(tool.preTransformed, 4);
        Assert.assertEquals(tool.postTransformed, 1);
        Assert.assertEquals(tool.totalReads, 5);
    }

}