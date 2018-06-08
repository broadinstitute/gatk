package org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorReferenceTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Class to test the {@link VcfFuncotationFactory}.
 * Created by jonn on 3/26/18.
 */
public class VcfFuncotationFactoryUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final String FACTORY_NAME    = "TestFactory";
    private static final String FACTORY_VERSION = "TEST_VERSION";

    //==================================================================================================================
    // Private Members:

    private static final ReferenceDataSource CHR3_REF_DATA_SOURCE = ReferenceDataSource.of( new File(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()).toPath() );
    
    private static final LinkedHashMap<String, Object> FIELD_DEFAULT_MAP = new LinkedHashMap<>();
    
    static {
       FIELD_DEFAULT_MAP.put("ASP","false");
       FIELD_DEFAULT_MAP.put("ASS","false");
       FIELD_DEFAULT_MAP.put("CAF","");
       FIELD_DEFAULT_MAP.put("CDA","false");
       FIELD_DEFAULT_MAP.put("CFL","false");
       FIELD_DEFAULT_MAP.put("COMMON","");
       FIELD_DEFAULT_MAP.put("DSS","false");
       FIELD_DEFAULT_MAP.put("G5","false");
       FIELD_DEFAULT_MAP.put("G5A","false");
       FIELD_DEFAULT_MAP.put("GENEINFO","");
       FIELD_DEFAULT_MAP.put("GNO","false");
       FIELD_DEFAULT_MAP.put("HD","false");
       FIELD_DEFAULT_MAP.put("INT","false");
       FIELD_DEFAULT_MAP.put("KGPhase1","false");
       FIELD_DEFAULT_MAP.put("KGPhase3","false");
       FIELD_DEFAULT_MAP.put("LSD","false");
       FIELD_DEFAULT_MAP.put("MTP","false");
       FIELD_DEFAULT_MAP.put("MUT","false");
       FIELD_DEFAULT_MAP.put("NOC","false");
       FIELD_DEFAULT_MAP.put("NOV","false");
       FIELD_DEFAULT_MAP.put("NSF","false");
       FIELD_DEFAULT_MAP.put("NSM","false");
       FIELD_DEFAULT_MAP.put("NSN","false");
       FIELD_DEFAULT_MAP.put("OM","false");
       FIELD_DEFAULT_MAP.put("OTH","false");
       FIELD_DEFAULT_MAP.put("PM","false");
       FIELD_DEFAULT_MAP.put("PMC","false");
       FIELD_DEFAULT_MAP.put("R3","false");
       FIELD_DEFAULT_MAP.put("R5","false");
       FIELD_DEFAULT_MAP.put("REF","false");
       FIELD_DEFAULT_MAP.put("RS","");
       FIELD_DEFAULT_MAP.put("RSPOS","");
       FIELD_DEFAULT_MAP.put("RV","false");
       FIELD_DEFAULT_MAP.put("S3D","false");
       FIELD_DEFAULT_MAP.put("SAO","");
       FIELD_DEFAULT_MAP.put("SLO","false");
       FIELD_DEFAULT_MAP.put("SSR","");
       FIELD_DEFAULT_MAP.put("SYN","false");
       FIELD_DEFAULT_MAP.put("TOPMED","");
       FIELD_DEFAULT_MAP.put("TPA","false");
       FIELD_DEFAULT_MAP.put("U3","false");
       FIELD_DEFAULT_MAP.put("U5","false");
       FIELD_DEFAULT_MAP.put("VC","");
       FIELD_DEFAULT_MAP.put("VLD","false");
       FIELD_DEFAULT_MAP.put("VP","");
       FIELD_DEFAULT_MAP.put("WGT","");
       FIELD_DEFAULT_MAP.put("WTD","false");
       FIELD_DEFAULT_MAP.put("dbSNPBuildID","");
    }

    //==================================================================================================================
    // Helper Methods:

    private VariantContext createVariantContext(final String contig,
                                                final int start,
                                                final int end,
                                                final String refString,
                                                final String altString) {

        final Allele refAllele = Allele.create(refString, true);
        final Allele altAllele = Allele.create(altString);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );

        return variantContextBuilder.make();
    }

    private Object[] helpProvideForTestCreateFuncotations(final String contig,
                                                          final int start,
                                                          final int end,
                                                          final String refAlleleString,
                                                          final String altAlleleString,
                                                          final List<Funcotation> expected) {
        return new Object[] {
                createVariantContext(contig, start, end, refAlleleString, altAlleleString),
                new ReferenceContext( CHR3_REF_DATA_SOURCE, new SimpleInterval(contig, start, end)),
                expected
        };
    }


    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Iterator<Object[]> provideForTestGetName() {
        final ArrayList<Object[]> data = new ArrayList<>();

        for ( int i = 0; i < 20 ; ++i ) {
            data.add( new Object[] { RandomStringUtils.randomAlphanumeric(i) } );
        }

        return data.iterator();
    }

    @DataProvider
    private Object[][] provideForTestCreateFuncotationsOnVariant() {

        return new Object[][] {
                // Trivial Case: No overlapping features:
                helpProvideForTestCreateFuncotations("3", 61650, 61650, "T", "C",
                        Collections.singletonList(
                                new TableFuncotation(FIELD_DEFAULT_MAP.keySet().stream().map(s->FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                        FIELD_DEFAULT_MAP.values().stream().map(Object::toString).collect(Collectors.toList()),
                                        Allele.create("C"), FACTORY_NAME)
                        )
                ),
                // One overlapping VCF feature:
                helpProvideForTestCreateFuncotations("3", 61662, 61662, "T", "C",
                        Collections.singletonList(
                                new TableFuncotation(FIELD_DEFAULT_MAP.keySet().stream().map(s->FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                        Arrays.asList("true","false","0.9744,0.02556","false","false","1","false","true","false","","false","true","false","true","true","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","73009205","61662","false","false","0","true","0","false","0.954392,0.0456075","false","false","false","SNV","true","0x05010000000515043e000100","1","false","130"),
                                        Allele.create("C"), FACTORY_NAME)
                        )
                ),
                // No matching VCF features (three overlap by position only), since there are no indels in dbSNP (the test datasource), so the ground truth should be a default entry, which was constructed here manually:
                helpProvideForTestCreateFuncotations("3", 64157, 64166, "AGAAAGGTCA", "TCTTTCCAGT",
                        Collections.singletonList(new TableFuncotation(FIELD_DEFAULT_MAP.keySet().stream().map(s->FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                Arrays.asList("false","false","","false","false","","false","false","false","","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","","","false","false","","false","","false","","false","false","false","","false","","","false",""),
                                Allele.create("TCTTTCCAGT"), FACTORY_NAME))
                ),
        };
    }

    //==================================================================================================================
    // Tests:

    @Test
    public void testGetAnnotationFeatureClass() {
        final VcfFuncotationFactory vcfFuncotationFactory = new VcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3));
        Assert.assertEquals(vcfFuncotationFactory.getAnnotationFeatureClass(), VariantContext.class);
    }

    @Test
    public void testGetType() {
        final VcfFuncotationFactory vcfFuncotationFactory = new VcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3));
        Assert.assertEquals(vcfFuncotationFactory.getType(), FuncotatorArgumentDefinitions.DataSourceType.VCF);
    }

    @Test(dataProvider = "provideForTestGetName")
    public void testGetName(final String name) {
        final VcfFuncotationFactory vcfFuncotationFactory = new VcfFuncotationFactory(name, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3));
        Assert.assertEquals(vcfFuncotationFactory.getName(), name);
    }

    @Test
    public void testGetSupportedFuncotationFields() {
        final VcfFuncotationFactory vcfFuncotationFactory =
                new VcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH));

        final LinkedHashSet<String> expectedFieldNames = new LinkedHashSet<>();

        for ( final String field : FIELD_DEFAULT_MAP.keySet()) {
            expectedFieldNames.add(FACTORY_NAME + "_" + field);
        }

        Assert.assertEquals( vcfFuncotationFactory.getSupportedFuncotationFields(), expectedFieldNames );
    }

    @Test(dataProvider = "provideForTestCreateFuncotationsOnVariant")
    public void testCreateFuncotationsOnVariant(final VariantContext variant,
                                                final ReferenceContext referenceContext,
                                                final List<Funcotation> expected) {

        // Make our factory:
        final VcfFuncotationFactory vcfFuncotationFactory =
                new VcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH));

        // Create features from the file:
        final List<Feature> vcfFeatures;
        try (final VCFFileReader vcfReader = new VCFFileReader(IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH))) {
            vcfFeatures = vcfReader.query( variant.getContig(), variant.getStart(), variant.getEnd() ).stream().collect(Collectors.toList() );
        }

        // Check to see if we're in business:
        Assert.assertEquals(
                vcfFuncotationFactory.createFuncotationsOnVariant(
                        variant,
                        referenceContext,
                        vcfFeatures
                ),
                expected
        );

        Assert.assertEquals(
                vcfFuncotationFactory.createFuncotationsOnVariant(
                        variant,
                        referenceContext,
                        vcfFeatures,
                        Collections.emptyList()
                ),
                expected
        );


    }
}
