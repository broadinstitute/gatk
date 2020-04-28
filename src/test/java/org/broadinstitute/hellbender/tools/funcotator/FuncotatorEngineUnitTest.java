package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.DummyPlaceholderGatkTool;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class FuncotatorEngineUnitTest extends GATKBaseTest {
    final static private String INPUT_VCF = FuncotatorTestConstants.FUNCOTATOR_TEST_DIR + "/PIK3CA_SNPS_engine_test_chr3.vcf";
    private static final String DS_PIK3CA_DIR = largeFileTestDir + "funcotator" + File.separator + "small_ds_pik3ca" + File.separator;
    @DataProvider
    public Object[][] provideGt() {
        return new Object[][] {
                // ground truth gene name, hasClinvarAnnotation
                {new File(INPUT_VCF), Arrays.asList("PIK3CA", "PIK3CA", "PIK3CA"), new boolean[]{true, true, false}}
        };
    }
    @Test(dataProvider = "provideGt")
    public void testGetFuncotationFactoriesAndCreateFuncotationMapForVariant(final File vcfFile,
                                                                             final List<String> correspondingGeneName,
                                                                             final boolean[] hasClinvarHit) {

        final Pair<VCFHeader, List<VariantContext>> entireVcf = VariantContextTestUtils.readEntireVCFIntoMemory(vcfFile.getAbsolutePath());
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths("hg19", Collections.singletonList(DS_PIK3CA_DIR));

        final Pair<VCFHeader, List<VariantContext>> vcfFileContents = VariantContextTestUtils.readEntireVCFIntoMemory(vcfFile.getAbsolutePath());

        // Set up our arguments:
        final FuncotatorVariantArgumentCollection funcotatorArguments = new FuncotatorVariantArgumentCollection();
        funcotatorArguments.referenceVersion = BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19;
        funcotatorArguments.transcriptSelectionMode = TranscriptSelectionMode.CANONICAL;
        funcotatorArguments.lookaheadFeatureCachingInBp = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE;

        // Create the metadata directly from the input.
        final FuncotatorEngine funcotatorEngine =
                new FuncotatorEngine(
                        funcotatorArguments,
                        vcfFileContents.getLeft().getSequenceDictionary(),
                        VcfFuncotationMetadata.create(new ArrayList<>(entireVcf.getLeft().getInfoHeaderLines())),
                        DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                                configData,
                                new LinkedHashMap<>(),
                                TranscriptSelectionMode.CANONICAL,
                                new HashSet<>(),
                                new DummyPlaceholderGatkTool(),
                                FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE,
                                new FlankSettings(0, 0),
                                false,
                                FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT)
                );

        for (int i = 0; i < entireVcf.getRight().size(); i++) {
            final VariantContext vc = entireVcf.getRight().get(i);
            final SimpleInterval variantInterval = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());
            final ReferenceContext referenceContext = new ReferenceContext(ReferenceDataSource.of(Paths.get(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref())), variantInterval);
            final FeatureContext featureContext = FuncotatorTestUtils.createFeatureContext(funcotatorEngine.getFuncotationFactories(), "TEST", variantInterval,
                    0,0,0, null);
            final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForVariant(vc, referenceContext, featureContext);

            // Check that all of the transcripts at this location have the same gene name as the corresponding gene.
            //  The ground truth selected has the same gene name for all transcripts.
            //  Also, input VCF has no multiallelics.
            for (final String txId : funcotationMap.getTranscriptList()) {
                Assert.assertEquals(funcotationMap.getFieldValue(txId, "Gencode_19_hugoSymbol", vc.getAlternateAllele(0)), correspondingGeneName.get(i));
                Assert.assertTrue((funcotationMap.getFieldValue(txId, "dummy_ClinVar_VCF_ALLELEID", vc.getAlternateAllele(0)).isEmpty()) != hasClinvarHit[i]);
            }
        }
    }
}
