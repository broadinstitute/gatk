package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;

/**
 *  * Unit test class for the {@link GencodeFuncotation} class.
 * Created by jonn on 9/1/17.
 */
public class GencodeFuncotationUnitTest extends GATKBaseTest {
    
    //==================================================================================================================
    // Helper Methods:
    private static GencodeFuncotation createGencodeFuncotation(final String hugoSymbol, final String ncbiBuild,
                                                               final String chromosome, final int start, final int end,
                                                               final GencodeFuncotation.VariantClassification variantClassification,
                                                               final GencodeFuncotation.VariantClassification secondaryVariantClassification,
                                                               final GencodeFuncotation.VariantType variantType,
                                                               final String refAllele,
                                                               final String tumorSeqAllele2, final String genomeChange,
                                                               final String annotationTranscript, final String transcriptStrand,
                                                               final Integer transcriptExon,
                                                               final Integer transcriptStartPos,
                                                               final Integer transcriptEndPos,
                                                               final String cDnaChange, final String codonChange,
                                                               final String proteinChange, final Double gcContent,
                                                               final String referenceContext,
                                                               final String geneTranscriptType,
                                                               final List<String> otherTranscripts) {

        final GencodeFuncotation gencodeFuncotation = new GencodeFuncotation();

        gencodeFuncotation.setVersion("TEST_VERSION");
        gencodeFuncotation.setDataSourceName(GencodeFuncotationFactory.DEFAULT_NAME);

        gencodeFuncotation.setHugoSymbol( hugoSymbol );
        gencodeFuncotation.setNcbiBuild( ncbiBuild );
        gencodeFuncotation.setChromosome( chromosome );
        gencodeFuncotation.setStart( start );
        gencodeFuncotation.setEnd( end );
        gencodeFuncotation.setVariantClassification( variantClassification );
        gencodeFuncotation.setSecondaryVariantClassification(secondaryVariantClassification);
        gencodeFuncotation.setVariantType( variantType );
        gencodeFuncotation.setRefAllele( refAllele );
        gencodeFuncotation.setTumorSeqAllele2( tumorSeqAllele2 );

        gencodeFuncotation.setGenomeChange( genomeChange );
        gencodeFuncotation.setAnnotationTranscript( annotationTranscript );
        gencodeFuncotation.setTranscriptStrand( transcriptStrand );
        gencodeFuncotation.setTranscriptExonNumber( transcriptExon );
        gencodeFuncotation.setTranscriptStartPos( transcriptStartPos );
        gencodeFuncotation.setTranscriptEndPos( transcriptEndPos );
        gencodeFuncotation.setcDnaChange( cDnaChange );
        gencodeFuncotation.setCodonChange( codonChange );
        gencodeFuncotation.setProteinChange( proteinChange );
        gencodeFuncotation.setGcContent( gcContent );
        gencodeFuncotation.setReferenceContext( referenceContext );
        gencodeFuncotation.setGeneTranscriptType( geneTranscriptType );
        gencodeFuncotation.setOtherTranscripts( otherTranscripts );

        return gencodeFuncotation;
    }

    private static GencodeFuncotation setFuncotationFieldOverride( final GencodeFuncotation gencodeFuncotation,
                                                                   final String field,
                                                                   final String overrideValue) {

        final GencodeFuncotation otherGencodeFuncotation = new GencodeFuncotation(gencodeFuncotation);

        otherGencodeFuncotation.setFieldSerializationOverrideValue( field, overrideValue );

        return otherGencodeFuncotation;
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestGetFieldNames() {
        //final GencodeFuncotation gencodeFuncotation, final LinkedHashSet<String> expected
        return new Object[][] {
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "big_%20_changes", "T1",
                                "3'", 1, 1, 1, "A", "ATC", "Lys", 1.0, null, "HPF000", Arrays.asList("ONE", "TWO", "THREE")),
                        new LinkedHashSet<>(
                                Arrays.asList("Gencode_TEST_VERSION_hugoSymbol",
                                        "Gencode_TEST_VERSION_ncbiBuild",
                                        "Gencode_TEST_VERSION_chromosome",
                                        "Gencode_TEST_VERSION_start",
                                        "Gencode_TEST_VERSION_end",
                                        "Gencode_TEST_VERSION_variantClassification",
                                        "Gencode_TEST_VERSION_secondaryVariantClassification",
                                        "Gencode_TEST_VERSION_variantType",
                                        "Gencode_TEST_VERSION_refAllele",
                                        "Gencode_TEST_VERSION_tumorSeqAllele1",
                                        "Gencode_TEST_VERSION_tumorSeqAllele2",
                                        "Gencode_TEST_VERSION_genomeChange",
                                        "Gencode_TEST_VERSION_annotationTranscript",
                                        "Gencode_TEST_VERSION_transcriptStrand",
                                        "Gencode_TEST_VERSION_transcriptExon",
                                        "Gencode_TEST_VERSION_transcriptPos",
                                        "Gencode_TEST_VERSION_cDnaChange",
                                        "Gencode_TEST_VERSION_codonChange",
                                        "Gencode_TEST_VERSION_proteinChange",
                                        "Gencode_TEST_VERSION_gcContent",
                                        "Gencode_TEST_VERSION_referenceContext",
                                        "Gencode_TEST_VERSION_geneTranscriptType",
                                        "Gencode_TEST_VERSION_otherTranscripts")
                        )
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetField() {
        //final GencodeFuncotation gencodeFuncotation, final String fieldName, final String expected

        final GencodeFuncotation gencodeFuncotation = createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                "A", "T", "big_%20_changes", "T1",
                "3'", 1, 1, 1,"A", "ATC", "Lys", 1.0, "ATGCGCAT", "HPF000", Arrays.asList("ONE", "TWO", "THREE"));

        return new Object[][] {
                {
                        gencodeFuncotation,
                        "Gencode_TEST_VERSION_hugoSymbol",
                        "TESTGENE"
                },
                {
                        gencodeFuncotation,
                        "hugoSymbol",
                        "TESTGENE"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_hugoSymbol", "GARBAGEDAY"),
                        "hugoSymbol",
                        "GARBAGEDAY"
                },
                {
                        gencodeFuncotation,
                        "ncbiBuild",
                        "BUILD1"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_ncbiBuild", "OVERRIDE"),
                        "ncbiBuild",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "chromosome",
                        "chr1"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_chromosome", "OVERRIDE"),
                        "chromosome",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "start",
                        "1"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_start", "OVERRIDE"),
                        "start",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "end",
                        "100"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_end", "OVERRIDE"),
                        "end",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "variantClassification",
                        GencodeFuncotation.VariantClassification.NONSENSE.toString()
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_variantClassification", "OVERRIDE"),
                        "variantClassification",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "secondaryVariantClassification",
                        GencodeFuncotation.VariantClassification.INTRON.toString()
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_secondaryVariantClassification", "OVERRIDE"),
                        "secondaryVariantClassification",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "variantType",
                        GencodeFuncotation.VariantType.SNP.toString()
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_variantType", "OVERRIDE"),
                        "variantType",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "refAllele",
                        "A"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_refAllele", "OVERRIDE"),
                        "refAllele",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "tumorSeqAllele1",
                        "A"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_tumorSeqAllele1", "OVERRIDE"),
                        "tumorSeqAllele1",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "tumorSeqAllele2",
                        "T"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_tumorSeqAllele2", "OVERRIDE"),
                        "tumorSeqAllele2",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "genomeChange",
                        "big_%20_changes"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_genomeChange", "OVERRIDE"),
                        "genomeChange",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "annotationTranscript",
                        "T1"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_annotationTranscript", "OVERRIDE"),
                        "annotationTranscript",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "transcriptStrand",
                        "3'"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_transcriptStrand", "OVERRIDE"),
                        "transcriptStrand",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "transcriptExon",
                        "1"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_transcriptExon", "OVERRIDE"),
                        "transcriptExon",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "transcriptPos",
                        "1"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_transcriptPos", "OVERRIDE"),
                        "transcriptPos",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "cDnaChange",
                        "A"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_cDnaChange", "OVERRIDE"),
                        "cDnaChange",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "codonChange",
                        "ATC"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_codonChange", "OVERRIDE"),
                        "codonChange",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "proteinChange",
                        "Lys"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_proteinChange", "OVERRIDE"),
                        "proteinChange",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "gcContent",
                        "1.0"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_gcContent", "OVERRIDE"),
                        "gcContent",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "referenceContext",
                        "ATGCGCAT"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_referenceContext", "OVERRIDE"),
                        "referenceContext",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "geneTranscriptType",
                        "HPF000"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_geneTranscriptType", "OVERRIDE"),
                        "geneTranscriptType",
                        "OVERRIDE"
                },
                {
                        gencodeFuncotation,
                        "otherTranscripts",
                        "ONE" + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER + "TWO" + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER + "THREE"
                },
                {
                        setFuncotationFieldOverride(gencodeFuncotation, "Gencode_TEST_VERSION_otherTranscripts", "OVERRIDE"),
                        "otherTranscripts",
                        "OVERRIDE"
                }
        };
    }

    @DataProvider
    Object[][] provideForTestGetFieldFail() {
        //final GencodeFuncotation gencodeFuncotation, final String fieldName, final String expected

        final GencodeFuncotation gencodeFuncotation = createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                "A", "T", "big_%20_changes", "T1",
                "3'", 1, 1, 1, "A", "ATC", "Lys", 1.0, "ATGCGCAT", "HPF000", Arrays.asList("ONE", "TWO", "THREE"));


        return new Object[][] {
                {
                    gencodeFuncotation,
                    "TESTFIELD_OMICRON"
                },
                {
                        gencodeFuncotation,
                        "hugoSymmbol"
                },
                {
                        gencodeFuncotation,
                        "GENCODE_hugoSymbol"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestGetFieldNames")
    public void testGetFieldNames(final GencodeFuncotation gencodeFuncotation, final LinkedHashSet<String> expected) {
        Assert.assertEquals(gencodeFuncotation.getFieldNames(), expected);
    }

    @Test(dataProvider = "provideForTestGetField")
    public void testGetField(final GencodeFuncotation gencodeFuncotation, final String fieldName, final String expected) {
        Assert.assertEquals(gencodeFuncotation.getField(fieldName), expected);
    }

    @Test(dataProvider = "provideForTestGetFieldFail", expectedExceptions = GATKException.class)
    public void testGetFieldFail(final GencodeFuncotation gencodeFuncotation, final String fieldName) {
        gencodeFuncotation.getField(fieldName);
    }
}
