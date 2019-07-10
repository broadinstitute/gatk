package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ReferenceContextFuncotationFactoryUnitTest extends GATKBaseTest {

    private static final ReferenceDataSource refDataSourceHg19Ch3;

    static {
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForCreateFuncotations() {

//        variant, altAllele, gtfFeature, reference, transcript, version

        final SimpleInterval variantInterval =  new SimpleInterval("chr3", 178921515, 178921517);
        final Allele refAllele = Allele.create("GCA", true);
        final Allele altAllele = Allele.create("TTG");
        final VariantContext variant = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                variantInterval.getContig(),
                variantInterval.getStart(),
                variantInterval.getEnd(),
                Arrays.asList(refAllele, altAllele)
        ).make();

        final String versionString = "VERSION";
        final String dataSourceName = "TEST_GENCODE_NAME";

        // ======================

        return new Object[][] {
                {
                        variant,
                        new ReferenceContext( refDataSourceHg19Ch3, variantInterval ),
                        versionString,
                        dataSourceName,
                        TableFuncotation.create(
                                Collections.singletonList(ReferenceContextFuncotationFactory.REFERENCE_CONTEXT_FIELD_NAME),
                                Collections.singletonList("TATAAATAGTGCACTCAGAATAA"),
                                altAllele,
                                dataSourceName,
                                FuncotationMetadataUtils.createWithUnknownAttributes(Collections.singletonList(ReferenceContextFuncotationFactory.REFERENCE_CONTEXT_FIELD_NAME))
                        ),
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test ( dataProvider = "provideForCreateFuncotations")
    void testCreateFuncotations(final VariantContext variant,
                                final ReferenceContext reference,
                                final String version,
                                final String dataSourceName,
                                final Funcotation expectedReferenceContextFuncotation) {

        ReferenceContextFuncotationFactory referenceContextFuncotationFactory = new ReferenceContextFuncotationFactory(version, dataSourceName, FuncotatorArgumentDefinitions.REFERENCE_CONTEXT_WINDOW_SIZE_DEFAULT_VALUE);
        List<Funcotation> referenceContextFuncotations = referenceContextFuncotationFactory.createFuncotations(variant, reference, null);

        Assert.assertEquals(referenceContextFuncotations.size(), 1);
        Assert.assertEquals(referenceContextFuncotations.get(0), expectedReferenceContextFuncotation );
    }
}
