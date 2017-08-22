package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 *  * Unit test class for the {@link GencodeFuncotation} class.
 * Created by jonn on 9/1/17.
 */
public class GencodeFuncotationUnitTest extends BaseTest {
    
    //==================================================================================================================
    // Helper Methods:
    private static GencodeFuncotation createFuncotation(final String hugoSymbol, final String ncbiBuild,
                                                        final String chromosome, final int start, final int end,
                                                        final GencodeFuncotation.VariantClassification variantClassification,
                                                        final GencodeFuncotation.VariantClassification secondaryVariantClassification,
                                                        final GencodeFuncotation.VariantType variantType,
                                                        final String refAllele, final String tumorSeqAllele1,
                                                        final String tumorSeqAllele2, final String genomeChange,
                                                        final String annotationTranscript, final String transcriptStrand,
                                                        final Integer transcriptExon, final Integer transcriptPos,
                                                        final String cDnaChange, final String codonChange,
                                                        final String proteinChange, final List<String> otherTranscripts) {

        final GencodeFuncotation gencodeFuncotation = new GencodeFuncotation();

        gencodeFuncotation.setHugoSymbol( hugoSymbol );
        gencodeFuncotation.setNcbiBuild( ncbiBuild );
        gencodeFuncotation.setChromosome( chromosome );
        gencodeFuncotation.setStart( start );
        gencodeFuncotation.setEnd( end );
        gencodeFuncotation.setVariantClassification( variantClassification );
        gencodeFuncotation.setSecondaryVariantClassification(secondaryVariantClassification);
        gencodeFuncotation.setVariantType( variantType );
        gencodeFuncotation.setRefAllele( refAllele );
        gencodeFuncotation.setTumorSeqAllele1( tumorSeqAllele1 );
        gencodeFuncotation.setTumorSeqAllele2( tumorSeqAllele2 );

        gencodeFuncotation.setGenomeChange( genomeChange );
        gencodeFuncotation.setAnnotationTranscript( annotationTranscript );
        gencodeFuncotation.setTranscriptStrand( transcriptStrand );
        gencodeFuncotation.setTranscriptExon( transcriptExon );
        gencodeFuncotation.setTranscriptPos( transcriptPos );
        gencodeFuncotation.setcDnaChange( cDnaChange );
        gencodeFuncotation.setCodonChange( codonChange );
        gencodeFuncotation.setProteinChange( proteinChange );
        gencodeFuncotation.setOtherTranscripts( otherTranscripts );

        return gencodeFuncotation;
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] createGencodeFuncotationsAndStringSerializations() {

        final String D = "|";

//        final String hugoSymbol, final String ncbiBuild,
//        final String chromosome, final int start, final int end,
//        final GencodeFuncotation.VariantClassification variantClassification,
//        final GencodeFuncotation.VariantType variantType,
//        final String refAllele, final String tumorSeqAllele1,
//        final String tumorSeqAllele2, final String genomeChange,
//        final String annotationTranscript, final String transcriptStrand,
//        final int transcriptExon, final int transcriptPos,
//        final String cDnaChange, final String codonChange,
//        final String proteinChange, final List<String> otherTranscripts

        return new Object[][] {
                {
                    createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                            GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                            "A", "T", "T", "big changes", "T1",
                            "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),
                        
                    "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                            GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                            "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                            "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", null),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D
                },
                {
                        createFuncotation(null, "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", null, "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", null, 50, 60,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + D + 50 + D + 60 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                null, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, null, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, null,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "G", "C", "C", null, "T1",
                                "3'", 1, 1, "A", "ACC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "G" + D + "C" + D + "C" + D + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ACC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", null,
                                null, 1, 1, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + D +
                                D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", null, null, "A", "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + D + D + "A" + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, null,  "ATC", "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + D + "ATC" + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", null, "Lys", Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + D + "Lys" + D + "ONE;TWO;THREE"
                },
                {
                        createFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", null, Arrays.asList("ONE", "TWO", "THREE")),

                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + D + "ONE;TWO;THREE"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test
    void testGetSerializedFieldNames() {
        final List<String> serializedFieldNames = GencodeFuncotation.getSerializedFieldNames();

        final List<String> expectedFieldNames = Arrays.asList(
                "hugoSymbol",
                "ncbiBuild",
                "chromosome",
                "start",
                "end",
                "variantClassification",
                "secondaryVariantClassification",
                "variantType",
                "refAllele",
                "tumorSeqAllele1",
                "tumorSeqAllele2",
                "genomeChange",
                "annotationTranscript",
                "transcriptStrand",
                "transcriptExon",
                "transcriptPos",
                "cDnaChange",
                "codonChange",
                "proteinChange",
                "otherTranscripts");

        Assert.assertEquals(serializedFieldNames, expectedFieldNames);
    }

    @Test(dataProvider = "createGencodeFuncotationsAndStringSerializations")
    void testSerializeToVcfString(final GencodeFuncotation gencodeFuncotation, final String expected) {
        Assert.assertEquals(gencodeFuncotation.serializeToVcfString(), expected);
    }
}
