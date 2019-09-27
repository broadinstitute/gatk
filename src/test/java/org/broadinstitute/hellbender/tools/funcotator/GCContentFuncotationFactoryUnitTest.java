package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeatureBaseData;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfGeneFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfTranscriptFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

public class GCContentFuncotationFactoryUnitTest extends GATKBaseTest {

    private static final double doubleEqualsEpsilon = 0.000001;

    private static final ReferenceDataSource refDataSourceHg19Ch3;

    static {
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideDataForTestCalculateGcContent() {

        // Base Position:
        //
        //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        //0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
        //0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        final String seq = "AAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCC";
        final String contig = "test";

        return new Object[][] {
                // MNPs:
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100), 10, 10.0/11.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  1, 1.0/3.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  9, 9.0/19.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 10, 11.0/21.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 50, 50.0/100.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),   1, 0.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),   9, 4.0 / 14.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),  10, 5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),   1, 1 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),   9, 10.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),  10, 10.0/16.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),   1, 6.0/8.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),   9, 10.0/24.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  10, 11.0/26.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  50, 50.0/100.0 },

                // Insertions:
                { Allele.create("A", true), Allele.create("AG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("AGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("AGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("CG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("CGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("CGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100), 10, 1 },
                { Allele.create("T", true), Allele.create("TG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  1, 1.0/2.0 },
                { Allele.create("T", true), Allele.create("TGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  9, 9.0/18.0 },
                { Allele.create("T", true), Allele.create("TGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 10, 10.0/20.0 },
                { Allele.create("T", true), Allele.create("TGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 49, 49.0/98.0 },

                // Deletions:
                { Allele.create("AAAAA",  true), Allele.create("A"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),  10,  5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),  10, 10.0/15.0  },
                { Allele.create("TGGGGG", true), Allele.create("T"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  10, 10.0/25.0 },
                { Allele.create("TG", true),     Allele.create("T"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 51),  10, 10.0/21.0 },
        };
    }

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
                                Collections.singletonList(GCContentFuncotationFactory.GC_CONTENT_FIELD_NAME),
                                Collections.singletonList(String.valueOf(0.3399503722084367)),
                                altAllele,
                                dataSourceName,
                                FuncotationMetadataUtils.createWithUnknownAttributes(Collections.singletonList(GCContentFuncotationFactory.GC_CONTENT_FIELD_NAME))
                        ),
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideDataForTestCalculateGcContent")
    void testCalculateGcContent(final Allele refAllele,
                                final Allele altAllele,
                                final ReferenceContext referenceContext,
                                final int windowSize,
                                final double expected) {
        Assert.assertEquals( GCContentFuncotationFactory.calculateGcContent( refAllele, altAllele, referenceContext, windowSize ), expected, doubleEqualsEpsilon);
    }

    @Test ( dataProvider = "provideForCreateFuncotations")
    void testCreateFuncotations(final VariantContext variant,
                                final ReferenceContext reference,
                                final String version,
                                final String dataSourceName,
                                final Funcotation expectedGCContentFuncotation) {

        GCContentFuncotationFactory gcContentFuncotationFactory = new GCContentFuncotationFactory(version, dataSourceName, FuncotatorArgumentDefinitions.GC_CONTENT_WINDOW_SIZE_DEFAULT_VALUE);
        List<Funcotation> gcContentFuncotations = gcContentFuncotationFactory.createFuncotations(variant, reference, null);

        Assert.assertEquals(gcContentFuncotations.size(), 1);
        Assert.assertEquals(gcContentFuncotations.get(0), expectedGCContentFuncotation );
    }
}
