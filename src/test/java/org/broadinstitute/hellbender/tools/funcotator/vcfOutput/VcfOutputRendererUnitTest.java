package org.broadinstitute.hellbender.tools.funcotator.vcfOutput;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Unit test class for the {@link VcfOutputRenderer} class.
 * Created by jonn on 9/1/17.
 */
public class VcfOutputRendererUnitTest extends GATKBaseTest {
    //==================================================================================================================
    // Static Variables:

    private static final String TEST_VCF = toolsTestDir + "/funcotator/hg38_trio.pik3ca.vcf";

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    //==================================================================================================================
    // Tests:

    /** Test that the exclusion list overrides the manually specified annotations */
    @Test
    public void testExclusionListOverridesManualDefaultAnnotations() {
        final Pair<VCFHeader, List<VariantContext>> entireInputVcf =  VariantContextTestUtils.readEntireVCFIntoMemory(TEST_VCF);
        final File outFile = createTempFile("vcf_output_renderer_exclusion", ".vcf");
        final VariantContextWriter vcfWriter = GATKVariantContextUtils.createVCFWriter(outFile.toPath(),null, false);

        final LinkedHashMap<String, String> dummyDefaults = new LinkedHashMap<>();
        dummyDefaults.put("FOO", "BAR");
        dummyDefaults.put("BAZ", "HUH?");

        final VcfOutputRenderer vcfOutputRenderer = new VcfOutputRenderer(vcfWriter,
            new ArrayList<>(), entireInputVcf.getLeft(), new LinkedHashMap<>(dummyDefaults),
            new LinkedHashMap<>(),
            new HashSet<>(), new HashSet<>(Arrays.asList("BAZ", "AC")), "Unknown");

        final VariantContext variant = entireInputVcf.getRight().get(0);

        final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(Collections.emptyList());
        vcfOutputRenderer.write(variant, funcotationMap);
        vcfOutputRenderer.close();

        // Check the output
        final Pair<VCFHeader, List<VariantContext>> tempVcf =  VariantContextTestUtils.readEntireVCFIntoMemory(outFile.getAbsolutePath());
        final VariantContext tempVariant = tempVcf.getRight().get(0);
        final String[] funcotatorKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(tempVcf.getLeft().getInfoHeaderLine("FUNCOTATION").getDescription());
        Assert.assertEquals(funcotatorKeys.length,1);
        Assert.assertEquals(funcotatorKeys[0],"FOO");
        final FuncotationMap tempFuncotationMap =
                FuncotationMap.createAsAllTableFuncotationsFromVcf(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, funcotatorKeys,
                        tempVariant.getAttributeAsString("FUNCOTATION", ""), tempVariant.getAlternateAllele(0), "TEST");
        Assert.assertTrue(tempFuncotationMap.get(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).get(0).hasField("FOO"));
        Assert.assertEquals(tempFuncotationMap.get(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).get(0).getField("FOO"), "BAR");
        Assert.assertFalse(tempFuncotationMap.get(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).get(0).hasField("BAZ"));

        // IMPORTANT: If the field is not a proper funcotation in VCFs, it will not be excluded.  I.e. if an input VCF has an excluded field, it will not be excluded.
        Assert.assertEquals(tempVariant.getAttribute("AC"), "1");
    }
}
