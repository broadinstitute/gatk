package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class SVVCFWriterUnitTest extends BaseTest{

    @Test
    public void testSortVariantsByCoordinate(){

        final String insOne = "AAA";new String(SVDiscoveryTestDataProvider.makeDummySequence(100, (byte)'A'));
        final String insTwo = "AAC";new String(SVDiscoveryTestDataProvider.makeDummySequence(100, (byte)'C'));

        final String contig = "21";
        final int pos = 100001;
        final int end = 100501;

        final VariantContext inversionOne = new VariantContextBuilder()
                .chr(contig).start(pos).stop(end)
                .alleles("G", SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_INV)
                .attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, insOne)
                .make();
        final VariantContext inversionTwo = new VariantContextBuilder()
                .chr(contig).start(pos).stop(end)
                .alleles("G", SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_INV)
                .attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, insTwo)
                .make();

        final VariantContext upstreamVariant = new VariantContextBuilder()
                .chr(contig).start(pos-50).stop(end)
                .alleles("T", SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DUP)
                .attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, insOne)
                .make();

        final VariantContext downstreamVariant = new VariantContextBuilder()
                .chr(contig).start(pos+20).stop(end+20)
                .alleles("C", SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_INS)
                .attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, insOne)
                .make();

        final File referenceDictionaryFile = new File(ReferenceUtils.getFastaDictionaryFileName(b37_reference_20_21));
        final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(referenceDictionaryFile);

        final List<VariantContext> sortedVariants = SVVCFWriter.sortVariantsByCoordinate(Arrays.asList(downstreamVariant, inversionTwo, inversionOne, upstreamVariant), dictionary);
        Assert.assertEquals(sortedVariants.size(), 4);
        Assert.assertEquals(sortedVariants.get(0), upstreamVariant);
        Assert.assertEquals(sortedVariants.get(1), inversionOne);
        Assert.assertEquals(sortedVariants.get(2), inversionTwo);
        Assert.assertEquals(sortedVariants.get(3), downstreamVariant);
    }
}
