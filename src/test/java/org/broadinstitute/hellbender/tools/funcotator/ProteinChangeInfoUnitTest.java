package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.stream.Collectors;

/**
 * Unit tests for the {@link ProteinChangeInfo} class.
 * Created by jonn on 10/23/18.
 */
public class ProteinChangeInfoUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestCreateProteinChangeInfo() {

        // PIK3CA is on chr3 and is transcribed in the FORWARD direction:
        final ReferenceDataSource pik3caTranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
        final String              pik3caFullTranscriptName   = pik3caTranscriptDataSource.getSequenceDictionary().getSequences().stream().filter(s -> s.getSequenceName().startsWith(FuncotatorTestConstants.PIK3CA_TRANSCRIPT)).map(SAMSequenceRecord::getSequenceName).collect(Collectors.joining());
        final ReferenceSequence   pik3caReferenceSequence    = pik3caTranscriptDataSource.queryAndPrefetch(pik3caFullTranscriptName, 158, 3364);

        // MUC16 is on chr19 and is transcribed in the REVERSE direction:
        // (The magic numbers here are the coding region for this transcript for MUC16)
        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
        final String              muc16FullTranscriptName   = muc16TranscriptDataSource.getSequenceDictionary().getSequences().stream().filter(s -> s.getSequenceName().startsWith(FuncotatorTestConstants.MUC16_TRANSCRIPT)).map(SAMSequenceRecord::getSequenceName).collect(Collectors.joining());
        final ReferenceSequence   muc16ReferenceSequence    = muc16TranscriptDataSource.queryAndPrefetch(muc16FullTranscriptName, 205, 43728);

        return new Object[][] {
                // ==========================================================
                // Standard Code:
                // ---------------------------
                // MNPs:
                //    length 1:
                {
                        Allele.create("A", true),
                        Allele.create("T"),
                        1,
                        1,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(1,1,"M", "L")
                },
                //    length >=3:
                {
                        Allele.create("CACT", true),
                        Allele.create("GGGA"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(543,544,"IT", "MG")
                },

                // Insertions:
                //     FS
                {
                        Allele.create("C", true),
                        Allele.create("CG"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,544,"T", "")
                },
                //     Non-FS between codons
                //         + strand
                {
                        Allele.create("C", true),
                        Allele.create("CGGA"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(543,544,"", "G")
                },
                //         - strand
                {
                        Allele.create("T", true),
                        Allele.create("TTGG"),
                        10083,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, false,
                        ProteinChangeInfo.create(3361,3362,"", "W")
                },
                //     Non-FS within codons
                {
                        Allele.create("A", true),
                        Allele.create("ACTG"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,545,"", "A")
                },
                {
                        Allele.create("A", true),
                        Allele.create("ATCG"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,544,"T", "IA")
                },

                // Deletions:
                //     FS
                {
                        Allele.create("CA", true),
                        Allele.create("C"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,544,"T", "")
                },
                //     Non-FS between codons
                //         + strand
                {
                        Allele.create("CACT", true),
                        Allele.create("C"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,544,"T", "")
                },
                //         - strand
                {
                        Allele.create("TCTG", true),
                        Allele.create("T"),
                        10083,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, false,
                        ProteinChangeInfo.create(3362,3362,"L", "")
                },
                {
                        Allele.create("TCTGAGC", true),
                        Allele.create("T"),
                        10083,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, false,
                        ProteinChangeInfo.create(3362,3362,"LS", "")
                },
                //     Non-FS within codons
                //         + strand
                {
                        Allele.create("ACTG", true),
                        Allele.create("A"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,545,"TE", "K")
                },
                {
                        Allele.create("CTGA", true),
                        Allele.create("C"),
                        1631,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(545,545,"E", "")
                },
                {
                        Allele.create("CTGAGCA", true),
                        Allele.create("C"),
                        1631,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(545,545,"EQ", "")
                },
                {
                        Allele.create("ACTGAGC", true),
                        Allele.create("A"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, false,
                        ProteinChangeInfo.create(544,546,"TEQ", "K")
                },
                //         - strand
                {
                        Allele.create("CTCT", true),
                        Allele.create("C"),
                        10082,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, false,
                        ProteinChangeInfo.create(3362,3362,"L", "")
                },
                {
                        Allele.create("TCTC", true),
                        Allele.create("T"),
                        10081,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, false,
                        ProteinChangeInfo.create(3361,3361,"S", "")
                },
                {
                        Allele.create("CTCTGAG", true),
                        Allele.create("C"),
                        10082,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, false,
                        ProteinChangeInfo.create(3362,3362,"LS", "")
                },

                // ==========================================================
                // Mitochondrial Code:
                // ---------------------------
                // MNPs:
                //    length 1:
                {
                        Allele.create("A", true),
                        Allele.create("T"),
                        1,
                        1,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(1,1,"M", "L")
                },
                //    length >=3:
                {
                        Allele.create("CACT", true),
                        Allele.create("GGGA"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(543,544,"IT", "MG")
                },

                // Insertions:
                //     FS
                {
                        Allele.create("C", true),
                        Allele.create("CG"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,544,"T", "")
                },
                //     Non-FS between codons
                //         + strand
                {
                        Allele.create("C", true),
                        Allele.create("CGGA"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(543,544,"", "G")
                },
                //         - strand
                {
                        Allele.create("T", true),
                        Allele.create("TTGG"),
                        10083,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, true,
                        ProteinChangeInfo.create(3361,3362,"", "W")
                },
                //     Non-FS within codons
                {
                        Allele.create("A", true),
                        Allele.create("ACTG"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,545,"", "A")
                },
                {
                        Allele.create("A", true),
                        Allele.create("ATCG"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,544,"T", "IA")
                },

                // Deletions:
                //     FS
                {
                        Allele.create("CA", true),
                        Allele.create("C"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,544,"T", "")
                },
                //     Non-FS between codons
                //         + strand
                {
                        Allele.create("CACT", true),
                        Allele.create("C"),
                        1629,
                        1627,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,544,"T", "")
                },
                //         - strand
                {
                        Allele.create("TCTG", true),
                        Allele.create("T"),
                        10083,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, true,
                        ProteinChangeInfo.create(3362,3362,"L", "")
                },
                {
                        Allele.create("TCTGAGC", true),
                        Allele.create("T"),
                        10083,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, true,
                        ProteinChangeInfo.create(3362,3362,"LS", "")
                },
                //     Non-FS within codons
                //         + strand
                {
                        Allele.create("ACTG", true),
                        Allele.create("A"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,545,"TE", "K")
                },
                {
                        Allele.create("CTGA", true),
                        Allele.create("C"),
                        1631,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(545,545,"E", "")
                },
                {
                        Allele.create("CTGAGCA", true),
                        Allele.create("C"),
                        1631,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(545,545,"EQ", "")
                },
                {
                        Allele.create("ACTGAGC", true),
                        Allele.create("A"),
                        1630,
                        1630,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(544,546,"TEQ", "K")
                },
                //         - strand
                {
                        Allele.create("CTCT", true),
                        Allele.create("C"),
                        10082,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, true,
                        ProteinChangeInfo.create(3362,3362,"L", "")
                },
                {
                        Allele.create("TCTC", true),
                        Allele.create("T"),
                        10081,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, true,
                        ProteinChangeInfo.create(3361,3361,"S", "")
                },
                {
                        Allele.create("CTCTGAG", true),
                        Allele.create("C"),
                        10082,
                        10081,
                        muc16ReferenceSequence.getBaseString(),
                        Strand.NEGATIVE, true,
                        ProteinChangeInfo.create(3362,3362,"LS", "")
                },

                // ==========================================================
                // Mitochondrial Code specific differences:
                // -----------------------------
                {
                        Allele.create("T", true),
                        Allele.create("A"),
                        207,
                        205,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(69,69,"I", "M")
                },
                {
                        Allele.create("A", true),
                        Allele.create("C"),
                        3207,
                        3205,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(1069,1069,"W", "C")
                },
                {
                        Allele.create("G", true),
                        Allele.create("A"),
                        3171,
                        3169,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(1057,1057,"W", "W")
                },
                {
                        Allele.create("G", true),
                        Allele.create("A"),
                        3171,
                        3169,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(1057,1057,"W", "W")
                },
                {
                        Allele.create("AA", true),
                        Allele.create("GG"),
                        3089,
                        3088,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(1030,1030,"K", "*")
                },
                {
                        Allele.create("A", true),
                        Allele.create("G"),
                        3089,
                        3088,
                        pik3caReferenceSequence.getBaseString(),
                        Strand.POSITIVE, true,
                        ProteinChangeInfo.create(1030,1030,"K", "*")
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestCreateProteinChangeInfo")
    void testCreateProteinChangeInfo( final Allele refAllele,
                                      final Allele altAllele,
                                      final int codingSequenceAlleleStart,
                                      final int alignedCodingSequenceAlleleStart,
                                      final String codingSequence,
                                      final Strand strand,
                                      final boolean isMitochondria,
                                      final ProteinChangeInfo expected ) {

        Assert.assertEquals(
                ProteinChangeInfo.create(
                        refAllele,
                        altAllele,
                        codingSequenceAlleleStart,
                        alignedCodingSequenceAlleleStart,
                        codingSequence,
                        strand,
                        isMitochondria),
                expected
        );
    }

}
