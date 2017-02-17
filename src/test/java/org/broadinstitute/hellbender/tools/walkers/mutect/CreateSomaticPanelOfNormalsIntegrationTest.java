package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.testng.Assert.*;

/**
 * Created by David Benjamin on 2/17/17.
 */
public class CreateSomaticPanelOfNormalsIntegrationTest extends CommandLineProgramTest {

    private static final File PON_VCFS_DIR = new File(publicTestDir, "org/broadinstitute/hellbender/tools/mutect/createpon/");

    /**
     * In the following test, we have sample1.vcf:
     * 20	577548	.	C	G	100	PASS	SOMATIC;VAF=0.48275862069;DPR=29.0	GT	0/1
     * 20	1838610	.	T	A	100	PASS	SOMATIC;VAF=0.5;DPR=35.3333333333	GT	0/1
     * 20	2916255	.	G	C	100	PASS	SOMATIC;VAF=0.5;DPR=20.0	GT	0/1
     * 20	5270544	.	C	A	100	PASS	SOMATIC;VAF=0.5;DPR=38.0	GT	0/1
     * 20	5758517	.	G	T	100	PASS	SOMATIC;VAF=0.478260869565;DPR=25.0	GT	0/1
     * 20	6947936	.	A	G	100	PASS	SOMATIC;VAF=0.5;DPR=32.0	GT	0/1
     * 20	7492891	.	G	A	100	PASS	SOMATIC;VAF=0.481481481481;DPR=28.3333333333	GT	0/1
     * 20	8957515	.	T	C	100	PASS	SOMATIC;VAF=0.5;DPR=27.6666666667	GT	0/1
     *
     * and sample2.vcf:
     * 20	577548	.	C	G	100	PASS	SOMATIC;VAF=0.48275862069;DPR=29.0	GT	0/1
     * 20	1838610	.	T	A	100	PASS	SOMATIC;VAF=0.5;DPR=35.3333333333	GT	0/1
     * 20	2916255	.	G	C	100	PASS	SOMATIC;VAF=0.5;DPR=20.0	GT	0/1
     * 20	7492891	.	G	A	100	PASS	SOMATIC;VAF=0.481481481481;DPR=28.3333333333	GT	0/1
     * 20	8957515	.	T	C	100	PASS	SOMATIC;VAF=0.5;DPR=27.6666666667	GT	0/1
     * 20	9080929	.	C	A	100	PASS	SOMATIC;VAF=0.487179487179;DPR=42.3333333333	GT	0/1
     * with overlap:
     * 20	577548	.	C	G
     * 20	1838610	.	T	A
     * 20	2916255	.	G	C
     * 20	7492891	.	G	A
     * 20	8957515	.	T	C

     * @throws IOException
     */
    @Test
    public void test() throws IOException {
        final File vcf1 = new File(PON_VCFS_DIR, "sample1.vcf");
        final File vcf2 = new File(PON_VCFS_DIR, "sample2.vcf");

        final File vcfInputFile = createTempFile("vcfs", ".list");
        FileUtils.writeLines(vcfInputFile, Arrays.asList(vcf1.getAbsolutePath(), vcf2.getAbsolutePath()));

        final File outputVcf = createTempFile("pon", ".vcf");
        final String[] args = {
                "-" + CreateSomaticPanelOfNormals.INPUT_VCFS_LIST_SHORT_NAME, vcfInputFile.getAbsolutePath(),
                "-O", outputVcf.getAbsolutePath()
        };

        runCommandLine(args);

        final List<VariantContext> ponVariants =
                StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                .collect(Collectors.toList());

        Assert.assertEquals(ponVariants.size(), 5);
        final VariantContext vc1 = ponVariants.get(0);
        final VariantContext vc5 = ponVariants.get(4);
        Assert.assertEquals(vc1.getStart(), 577548);
        Assert.assertEquals(vc1.getNAlleles(), 2);
        Assert.assertTrue(vc1.getAlternateAllele(0).basesMatch("G"));
        Assert.assertEquals(vc5.getStart(), 8957515);
        Assert.assertEquals(vc5.getNAlleles(), 2);
        Assert.assertTrue(vc5.getAlternateAllele(0).basesMatch("C"));
    }

}