package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
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
                                                               final String refAllele, final String tumorSeqAllele1,
                                                               final String tumorSeqAllele2, final String genomeChange,
                                                               final String annotationTranscript, final String transcriptStrand,
                                                               final Integer transcriptExon, final Integer transcriptPos,
                                                               final String cDnaChange, final String codonChange,
                                                               final String proteinChange, final Double gcContent,
                                                               final List<String> otherTranscripts) {

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
        gencodeFuncotation.setTranscriptExonNumber( transcriptExon );
        gencodeFuncotation.setTranscriptPos( transcriptPos );
        gencodeFuncotation.setcDnaChange( cDnaChange );
        gencodeFuncotation.setCodonChange( codonChange );
        gencodeFuncotation.setProteinChange( proteinChange );
        gencodeFuncotation.setGcContent( gcContent );
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
    Object[][] provideDataForTestSerializationOverrides() {

        final GencodeFuncotation gencodeFuncotation = createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                "A", "T", "T", "big changes", "T1",
                "3'", 1, 1, "A", "ATC", "Lys", 1.0, Arrays.asList("ONE", "TWO", "THREE"));

        final String overrideVal = "GARBAGEDAY!";

        final String D = "|";

        return new Object[][] {
                { setFuncotationFieldOverride(gencodeFuncotation, "hugoSymbol", overrideVal),
                        overrideVal + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "ncbiBuild", overrideVal),
                        "TESTGENE" + D + overrideVal + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "chromosome", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + overrideVal + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "start", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + overrideVal + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "end", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + overrideVal + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "variantClassification", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                overrideVal + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "secondaryVariantClassification", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + overrideVal + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "variantType", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + overrideVal + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "refAllele", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                overrideVal + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "tumorSeqAllele1", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + overrideVal + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "tumorSeqAllele2", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + overrideVal + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "genomeChange", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + overrideVal + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "annotationTranscript", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + overrideVal + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "transcriptStrand", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                overrideVal + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "transcriptExon", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + overrideVal + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "transcriptPos", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + overrideVal + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "cDnaChange", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + overrideVal + D + "ATC" + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "codonChange", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + overrideVal + D + "Lys" + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "proteinChange", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + overrideVal + D + "1.0" + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "gcContent", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + overrideVal + D + "ONE;TWO;THREE"
                },
                { setFuncotationFieldOverride(gencodeFuncotation, "otherTranscripts", overrideVal),
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "1.0" + D + overrideVal
                },
        };
    }

    @DataProvider
    Object[][] createGencodeFuncotationsAndStringSerializations() {

        final String D = "|";

        return new Object[][] {
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", null, Arrays.asList("ONE", "TWO", "THREE")),
                        "AEIOUY!!!!!!!",
                        "AEIOUY!!!!!!!TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", null, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + D + "ONE;TWO;THREE"
                },
                {
                    createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                            GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                            "A", "T", "T", "big changes", "T1",
                            "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        "",
                    "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                            GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                            "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                            "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, null),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D
                },
                {
                        createGencodeFuncotation(null, "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", null, "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", null, 50, 60,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + D + 50 + D + 60 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                null, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, null, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, null,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "G", "C", "C", null, "T1",
                                "3'", 1, 1, "A", "ACC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "G" + D + "C" + D + "C" + D + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ACC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", null,
                                null, 1, 1, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + D +
                                D + "1" + D + 1 + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", null, null, "A", "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + D + D + "A" + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, null,  "ATC", "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + D + "ATC" + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", null, "Lys", 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + D + "Lys" + D + "0.332" + D + "ONE;TWO;THREE"
                },
                {
                        createGencodeFuncotation("TESTGENE", "BUILD1", "chr1", 1, 100,
                                GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.SNP,
                                "A", "T", "T", "big changes", "T1",
                                "3'", 1, 1, "A", "ATC", null, 0.332, Arrays.asList("ONE", "TWO", "THREE")),
                        null,
                        "TESTGENE" + D + "BUILD1" + D + "chr1" + D + 1 + D + 100 + D +
                                GencodeFuncotation.VariantClassification.NONSENSE + D + GencodeFuncotation.VariantClassification.INTRON + D + GencodeFuncotation.VariantType.SNP + D +
                                "A" + D + "T" + D + "T" + D + "big changes" + D + "T1" + D +
                                "3'" + D + "1" + D + 1 + D + "A" + D + "ATC" + D + D + "0.332" + D + "ONE;TWO;THREE"
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
                "gcContent",
                "otherTranscripts");

        Assert.assertEquals(serializedFieldNames, expectedFieldNames);
    }

    @Test(dataProvider = "provideDataForTestSerializationOverrides")
    void testSerializationOverrides(final GencodeFuncotation gencodeFuncotation, final String expected) {
        Assert.assertEquals(gencodeFuncotation.serializeToVcfString(), expected);
    }

    @Test(dataProvider = "createGencodeFuncotationsAndStringSerializations")
    void testSerializeToVcfString(final GencodeFuncotation gencodeFuncotation, final String manualAnnotationString, final String expected) {
        Assert.assertEquals(gencodeFuncotation.serializeToVcfString(manualAnnotationString), expected);
    }
}
