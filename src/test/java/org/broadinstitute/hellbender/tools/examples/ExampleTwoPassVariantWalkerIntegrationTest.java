package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import static org.broadinstitute.hellbender.tools.examples.ExampleTwoPassVariantWalker.COPY_OF_QD_KEY_NAME;
import static org.broadinstitute.hellbender.tools.examples.ExampleTwoPassVariantWalker.QD_DISTANCE_FROM_MEAN;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.QUAL_BY_DEPTH_KEY;

public class ExampleTwoPassVariantWalkerIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test() throws IOException {
        final File outputVcf = File.createTempFile("output", "vcf");
        final String inputVcf = "src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/haploid-multisample.vcf";

        final String[] args = {
                "-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME, inputVcf,
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputVcf.getAbsolutePath()
        };

        runCommandLine(args);

        final FeatureDataSource<VariantContext> variantsBefore = new FeatureDataSource<>(inputVcf);
        final FeatureDataSource<VariantContext> variantsAfter = new FeatureDataSource<>(outputVcf);
        final Iterator<VariantContext> beforeIterator = variantsBefore.iterator();
        final Iterator<VariantContext> afterIterator = variantsAfter.iterator();


        while (beforeIterator.hasNext() && afterIterator.hasNext()){
            final VariantContext before = beforeIterator.next();
            final VariantContext after = afterIterator.next();

            Assert.assertEquals(before.getAttributeAsDouble(QUAL_BY_DEPTH_KEY, -1.0),
                    after.getAttributeAsDouble(COPY_OF_QD_KEY_NAME, -2.0));
            Assert.assertTrue(after.getAttributeAsDouble(QD_DISTANCE_FROM_MEAN, -1.0) > 0);
        }

        Assert.assertTrue(! afterIterator.hasNext() && ! beforeIterator.hasNext());
    }

}