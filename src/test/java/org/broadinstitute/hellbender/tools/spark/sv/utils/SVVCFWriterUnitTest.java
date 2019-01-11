package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIDHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.createBracketedSymbAlleleString;

public class SVVCFWriterUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testSortVariantsByCoordinate(){

        final String insOne = "AAA";new String(TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'A'));
        final String insTwo = "AAC";new String(TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(100, (byte)'C'));

        final String contig = "21";
        final int pos = 100001;
        final int end = 100501;

        final VariantContext inversionOne = new VariantContextBuilder()
                .chr(contig).start(pos).stop(end)
                .alleles("G", createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV))
                .attribute(GATKSVVCFConstants.INSERTED_SEQUENCE, insOne)
                .make();
        final VariantContext inversionTwo = new VariantContextBuilder()
                .chr(contig).start(pos).stop(end)
                .alleles("G", createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV))
                .attribute(GATKSVVCFConstants.INSERTED_SEQUENCE, insTwo)
                .make();

        final VariantContext upstreamVariant = new VariantContextBuilder()
                .chr(contig).start(pos-50).stop(end)
                .alleles("T", createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP))
                .attribute(GATKSVVCFConstants.INSERTED_SEQUENCE, insOne)
                .make();

        final VariantContext downstreamVariant = new VariantContextBuilder()
                .chr(contig).start(pos+20).stop(end+20)
                .alleles("C", createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS))
                .attribute(GATKSVVCFConstants.INSERTED_SEQUENCE, insOne)
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

    @Test(groups = "sv")
    public void testSetHeader() {
        SAMSequenceDictionary referenceSequenceDictionary = new ReferenceMultiSparkSource(
                b37_2bit_reference_20_21 , ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);
        final VCFHeader vcfHeader = SVVCFWriter.getVcfHeader(referenceSequenceDictionary);
        Assert.assertNotNull(vcfHeader.getSequenceDictionary());

        final List<String> refContigs = vcfHeader.getContigLines().stream().map(VCFContigHeaderLine::getID).collect(Collectors.toList());
        Assert.assertTrue(refContigs.size()==2);

        final List<String> headerKeys = vcfHeader.getIDHeaderLines().stream().map(VCFIDHeaderLine::getID).sorted().collect(Collectors.toList());
        Assert.assertTrue(headerKeys.remove(VCFConstants.END_KEY));
        Assert.assertTrue(headerKeys.removeAll(refContigs));
        Assert.assertEquals(headerKeys, GATKSVVCFUtilsUnitTest.expectedHeaderKeysInVCF);
    }
}
