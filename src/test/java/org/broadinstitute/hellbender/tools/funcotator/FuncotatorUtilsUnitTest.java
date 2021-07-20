package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.collections.MapUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A unit test suite for the {@link FuncotatorUtils} class.
 * Created by jonn on 9/1/17.
 */
public class FuncotatorUtilsUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Static Variables:
    private static final Path TEST_REFERENCE = IOUtils.getPath(hg19MiniReference);
    private static final String TEST_REFERENCE_CONTIG = "1";
    private static final int TEST_REFERENCE_START = 12000;
    private static final int TEST_REFERENCE_END = 16000;

    private static final ReferenceDataSource refDataSourceHg19Ch3;

    private static String hg19Chr3Ref;

    // Initialization of static variables:
    static {

        hg19Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref();

        refDataSourceHg19Ch3 = ReferenceDataSource.of(IOUtils.getPath(hg19Chr3Ref));
    }

    //==================================================================================================================
    // Helper Methods:

    /**
     * Prints the bases of the {@link FuncotatorUtilsUnitTest#TEST_REFERENCE} from {@link FuncotatorUtilsUnitTest#TEST_REFERENCE_START} to {@link FuncotatorUtilsUnitTest#TEST_REFERENCE_END}
     * The print out has each base numbered - the results must be interpreted vertically (i.e. read the numbers from top to bottom to get the index of the base).
     * For example, the first 21 bases are as follows (with labels for the example):
     *
     *   Base 5       Base 19
     *     |             |
     * 000000000000000000000
     * 000000000000000000000
     * 000000000011111111112
     * 123456789012345678901
     * TCATCTGCAGGTGTCTGACTT
     *
     */
    private void printReferenceBases() {
        printReferenceBases(TEST_REFERENCE, TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END);
    }

    /**
     * Writes the bases of the {@link FuncotatorUtilsUnitTest#TEST_REFERENCE} from {@link FuncotatorUtilsUnitTest#TEST_REFERENCE_START} to {@link FuncotatorUtilsUnitTest#TEST_REFERENCE_END} to a file.
     * The print out has each base numbered - the results must be interpreted vertically (i.e. read the numbers from top to bottom to get the index of the base).
     * For example, the first 21 bases are as follows (with labels for the example):
     *
     *   Base 5       Base 19
     *     |             |
     * 000000000000000000000
     * 000000000000000000000
     * 000000000011111111112
     * 012345678901234567890
     * TCATCTGCAGGTGTCTGACTT
     *
     * @param contig Contig from which to print bases.
     * @param start Start point in contig from which to print bases.
     * @param end End point in contig to which to print bases.
     */
    private void printReferenceBases(final Path refFile, final String contig, final int start, final int end) {
        final ReferenceContext ref = new ReferenceContext(new ReferenceFileSource(refFile), new SimpleInterval(contig, start, end));

        final int numReferenceRows = (int)Math.ceil(Math.log10( end - start ) + 1);
        final ArrayList<StringBuilder> referenceIndexStringBuilders = new ArrayList<>();

        for ( int j = 0; j < numReferenceRows ; ++j ) {
            final StringBuilder sb = new StringBuilder();
            for ( int i = 0; i < ref.getBases().length; ++i ) {
                sb.append((int)(i / Math.pow(10, j)) % 10);
            }
            referenceIndexStringBuilders.add(sb);
        }


        try (Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("refBases_"+contig+"_"+start+"-"+end+".txt")))) {
            writer.write("Location: " + contig + ":" + start + ":" + end + "\n");
            writer.write("=================================================================================\n");
            for ( final StringBuilder sb : referenceIndexStringBuilders ) {
                writer.write( sb.toString() );
                writer.write( "\n" );
            }
            writer.write( new String(ref.getBases()) + "\n\n" );
        }
        catch ( final IOException ex ) {
            throw new GATKException("Could not create an output file!", ex);
        }
    }

    private static Object[] helpCreateDataForTestGetBasesInWindowAroundReferenceAllele(final String refAlleleBases,
                                                                                       final String strand,
                                                                                       final int windowSizeInBases,
                                                                                       final int startPos,
                                                                                       final int endPos,
                                                                                       final String expected) {
        return new Object[] {
            Allele.create(refAlleleBases, true),
            new ReferenceContext( refDataSourceHg19Ch3, new SimpleInterval("chr3", startPos, endPos) ),
            Strand.decode(strand),
            windowSizeInBases,
            expected
        };
    }

    private static Object[] helpProvideForCreateIntronicCDnaString(final int start, final List<Locatable> exonList, final String ref, final String alt, final Integer index, final Integer dist) {
        if ( index == null ) {
            return new Object[]{ start, exonList, ref, alt, "NA" };
        }
        else {
            return new Object[]{ start, exonList, ref, alt, "c.e" + (index+1) + (dist > 0 ? '+' : '-') + Math.abs(dist) + ref + '>' + alt };
        }
    }


    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideReferenceAndExonListAndExpected() {

        return new Object[][] {
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.emptyList(),
                        Strand.POSITIVE,
                        ""
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550)),
                        Strand.POSITIVE,
                        "GCAGAGACGGGAGGGGCAGAGCCGCAGGCACAGCCAAGAGGGCTGAAGAAA"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Arrays.asList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550),
                                new SimpleInterval("1", TEST_REFERENCE_START + 551, TEST_REFERENCE_START + 600)
                        ),
                        Strand.POSITIVE,
                        "GCAGAGACGGGAGGGGCAGAGCCGCAGGCACAGCCAAGAGGGCTGAAGAAATGGTAGAACGGAGCAGCTGGTGATGTGTGGGCCCACCGGCCCCAGGCTCC"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Arrays.asList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 500),
                                new SimpleInterval("1", TEST_REFERENCE_START + 501, TEST_REFERENCE_START + 501),
                                new SimpleInterval("1", TEST_REFERENCE_START + 502, TEST_REFERENCE_START + 502),
                                new SimpleInterval("1", TEST_REFERENCE_START + 503, TEST_REFERENCE_START + 503),
                                new SimpleInterval("1", TEST_REFERENCE_START + 504, TEST_REFERENCE_START + 504),
                                new SimpleInterval("1", TEST_REFERENCE_START + 505, TEST_REFERENCE_START + 505),
                                new SimpleInterval("1", TEST_REFERENCE_START + 506, TEST_REFERENCE_START + 506),
                                new SimpleInterval("1", TEST_REFERENCE_START + 507, TEST_REFERENCE_START + 507),
                                new SimpleInterval("1", TEST_REFERENCE_START + 508, TEST_REFERENCE_START + 508),
                                new SimpleInterval("1", TEST_REFERENCE_START + 509, TEST_REFERENCE_START + 509)
                        ),
                        Strand.POSITIVE,
                        "GCAGAGACGG"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 500)
                        ),
                        Strand.POSITIVE,
                        "G"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, 1, 10)),
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 1)
                        ),
                        Strand.POSITIVE,
                        FuncotatorConstants.MASKED_ANY_BASE_STRING
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.emptyList(),
                        Strand.NEGATIVE,
                        ""
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550)),
                        Collections.singletonList(new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550)),
                        Strand.NEGATIVE,
                        "TTTCTTCAGCCCTCTTGGCTGTGCCTGCGGCTCTGCCCCTCCCGTCTCTGC"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 600)),
                        Arrays.asList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550),
                                new SimpleInterval("1", TEST_REFERENCE_START + 551, TEST_REFERENCE_START + 600)
                        ),
                        Strand.NEGATIVE,
                        "GGAGCCTGGGGCCGGTGGGCCCACACATCACCAGCTGCTCCGTTCTACCATTTCTTCAGCCCTCTTGGCTGTGCCTGCGGCTCTGCCCCTCCCGTCTCTGC"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Arrays.asList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 509, TEST_REFERENCE_START + 509),
                                new SimpleInterval("1", TEST_REFERENCE_START + 508, TEST_REFERENCE_START + 508),
                                new SimpleInterval("1", TEST_REFERENCE_START + 507, TEST_REFERENCE_START + 507),
                                new SimpleInterval("1", TEST_REFERENCE_START + 506, TEST_REFERENCE_START + 506),
                                new SimpleInterval("1", TEST_REFERENCE_START + 505, TEST_REFERENCE_START + 505),
                                new SimpleInterval("1", TEST_REFERENCE_START + 504, TEST_REFERENCE_START + 504),
                                new SimpleInterval("1", TEST_REFERENCE_START + 503, TEST_REFERENCE_START + 503),
                                new SimpleInterval("1", TEST_REFERENCE_START + 502, TEST_REFERENCE_START + 502),
                                new SimpleInterval("1", TEST_REFERENCE_START + 501, TEST_REFERENCE_START + 501),
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 500)
                        ),
                        Strand.NEGATIVE,
                        "TGGTCAGCCACTGCAGCC"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Arrays.asList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 500),
                                new SimpleInterval("1", TEST_REFERENCE_START + 501, TEST_REFERENCE_START + 501),
                                new SimpleInterval("1", TEST_REFERENCE_START + 502, TEST_REFERENCE_START + 502),
                                new SimpleInterval("1", TEST_REFERENCE_START + 503, TEST_REFERENCE_START + 503),
                                new SimpleInterval("1", TEST_REFERENCE_START + 504, TEST_REFERENCE_START + 504),
                                new SimpleInterval("1", TEST_REFERENCE_START + 505, TEST_REFERENCE_START + 505),
                                new SimpleInterval("1", TEST_REFERENCE_START + 506, TEST_REFERENCE_START + 506),
                                new SimpleInterval("1", TEST_REFERENCE_START + 507, TEST_REFERENCE_START + 507),
                                new SimpleInterval("1", TEST_REFERENCE_START + 508, TEST_REFERENCE_START + 508),
                                new SimpleInterval("1", TEST_REFERENCE_START + 509, TEST_REFERENCE_START + 509)
                        ),
                        Strand.NEGATIVE,
                        "TGGTCAGCCACTGCAGCC"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START + 500, TEST_REFERENCE_END + 500)),
                        Collections.singletonList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 500)
                        ),
                        Strand.NEGATIVE,
                        "C"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 500)
                        ),
                        Strand.NEGATIVE,
                        "T"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("1", TEST_REFERENCE_START, TEST_REFERENCE_START)
                        ),
                        Strand.NEGATIVE,
                        "C"
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, 1, 10)),
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 1)
                        ),
                        Strand.NEGATIVE,
                        FuncotatorConstants.MASKED_ANY_BASE_STRING
                },
        };
    }

    @DataProvider
    Object[][] provideReferenceAndExonListForGatkExceptions() {

        return new Object[][] {
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("2", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550)
                        ),
                        Strand.POSITIVE
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("2", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550)
                        ),
                        Strand.NEGATIVE
                },
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("2", TEST_REFERENCE_START + 500, TEST_REFERENCE_START + 550)
                        ),
                        Strand.NONE
                },
        };
    }

    @DataProvider
    Object[][] provideReferenceAndExonListForIllegalArgumentExceptions() {

        return new Object[][] {
                {
                        new ReferenceContext(new ReferenceFileSource(TEST_REFERENCE), new SimpleInterval(TEST_REFERENCE_CONTIG, TEST_REFERENCE_START, TEST_REFERENCE_END)),
                        Collections.singletonList(
                                new SimpleInterval("1", TEST_REFERENCE_START + 500, TEST_REFERENCE_START)
                        ),
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetIndelAdjustedAlleleChangeStartPosition() {

        return new Object[][] {
                {
                        new VariantContextBuilder().chr("1").start(1).stop(1).alleles(Arrays.asList(Allele.create("A", true), Allele.create("T"))).make(),
                        1
                },
                {
                        new VariantContextBuilder().chr("1").start(1).stop(1).alleles(Arrays.asList(Allele.create("A", true), Allele.create("TG"))).make(),
                        1
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1234).alleles(Arrays.asList(Allele.create("A", true), Allele.create("AT"))).make(),
                        1235
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1234).alleles(Arrays.asList(Allele.create("A", true), Allele.create("ATG"))).make(),
                        1235
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1234).alleles(Arrays.asList(Allele.create("A", true), Allele.create("ATGTC"))).make(),
                        1235
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1236).alleles(Arrays.asList(Allele.create("ATG", true), Allele.create("AT"))).make(),
                        1236
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1236).alleles(Arrays.asList(Allele.create("ATG", true), Allele.create("A"))).make(),
                        1235
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1237).alleles(Arrays.asList(Allele.create("ATGT", true), Allele.create("ATGTC"))).make(),
                        1238
                },
                {
                        new VariantContextBuilder().chr("1").start(1234).stop(1238).alleles(Arrays.asList(Allele.create("ATGTC", true), Allele.create("ATGT"))).make(),
                        1238
                },

        };
    }

    @DataProvider
    Object[][] provideDataForTestGetNonOverlappingAltAlleleBaseString() {
        return new Object[][] {
                { Allele.create("A", true),          Allele.create("A"),          false, "" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("AAAAAAAAAA"), false, "" },

                { Allele.create("AAAAAAAAAA", true), Allele.create("TAAAAAAAAA"), false, "T" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("AAAAAAAAAT"), false, "T" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("ATTTTTTTTT"), false, "TTTTTTTTT" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("TTTTTTTTTA"), false, "TTTTTTTTT" },

                { Allele.create("AAAAA", true), Allele.create("TAAAAAAAAA"),      false, "TAAAA" },
                { Allele.create("AAAAA", true), Allele.create("AAAAAAAAAT"),      false, "AAAAT" },
                { Allele.create("AAAAA", true), Allele.create("ATTTTTTTTT"),      false, "TTTTTTTTT" },
                { Allele.create("AAAAA", true), Allele.create("TTTTTTTTTA"),      false, "TTTTTTTTT" },

                { Allele.create("AAAAAAAAAA", true), Allele.create("TAAAA"),      false, "T" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("AAAAT"),      false, "T" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("ATTTT"),      false, "TTTT" },
                { Allele.create("AAAAAAAAAA", true), Allele.create("TTTTA"),      false, "TTTT" },

//                { Allele.create("A", true),          Allele.create("A"),          true,  "" },
//                { Allele.create("AAAAAAAAAA", true), Allele.create("AAAAAAAAAA"), true,  "" },
//
//                { Allele.create("AAAAAAAAAA", true), Allele.create("TAAAAAAAAA"), true,  "T" },
//                { Allele.create("AAAAAAAAAA", true), Allele.create("AAAAAAAAAT"), true,  "T" },
//                { Allele.create("AAAAAAAAAA", true), Allele.create("ATTTTTTTTT"), true,  "TTTTTTTTT" },
//                { Allele.create("AAAAAAAAAA", true), Allele.create("TTTTTTTTTA"), true,  "TTTTTTTTT" },
        };
    }

    @DataProvider
    Object[][] provideDataForGetStartPositionInTranscript() {

        final List<? extends Locatable> exons_forward = Arrays.asList(
                new SimpleInterval("chr1", 10,19),
                new SimpleInterval("chr1", 30,39),
                new SimpleInterval("chr1", 50,59),
                new SimpleInterval("chr1", 70,79),
                new SimpleInterval("chr1", 90,99)
        );

        final List<? extends Locatable> exons_backward = Arrays.asList(
                new SimpleInterval("chr1", 90,99),
                new SimpleInterval("chr1", 70,79),
                new SimpleInterval("chr1", 50,59),
                new SimpleInterval("chr1", 30,39),
                new SimpleInterval("chr1", 10,19)
        );

        // And a spot-check test from real data:
        final List<? extends Locatable> spot_check_exons = Arrays.asList(
                new SimpleInterval("3", 178916614, 178916965),
                new SimpleInterval("3", 178917478, 178917687),
                new SimpleInterval("3", 178919078, 178919328),
                new SimpleInterval("3", 178921332, 178921577),
                new SimpleInterval("3", 178922291, 178922376),
                new SimpleInterval("3", 178927383, 178927488),
                new SimpleInterval("3", 178927974, 178928126),
                new SimpleInterval("3", 178928219, 178928353),
                new SimpleInterval("3", 178935998, 178936122)
        );

        return new Object[][] {
                { new SimpleInterval("chr1", 1, 1),     exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 25, 67),   exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 105, 392), exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 10, 10),   exons_forward, Strand.POSITIVE,  1 },
                { new SimpleInterval("chr1", 99, 99),   exons_forward, Strand.POSITIVE, 50 },
                { new SimpleInterval("chr1", 50, 67),   exons_forward, Strand.POSITIVE, 21 },
                { new SimpleInterval("chr1", 67, 75),   exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 91, 97),   exons_forward, Strand.POSITIVE, 42 },

                { new SimpleInterval("chr1", 1, 1),     exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 25, 67),   exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 105, 392), exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 10, 10),   exons_backward, Strand.NEGATIVE, 50 },
                { new SimpleInterval("chr1", 99, 99),   exons_backward, Strand.NEGATIVE,  1 },
                { new SimpleInterval("chr1", 50, 67),   exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 67, 75),   exons_backward, Strand.NEGATIVE, 15 },

                // Spot check:
                { new SimpleInterval("3", 178936090, 178936096), spot_check_exons, Strand.POSITIVE, 1632 }
        };
    }

    @DataProvider
    Object[][] providePositionAndExpectedAlignedPosition() {
        return new Object[][] {
                {   1,   1},
                {   2,   1},
                {   3,   1},
                {   4,   4},
                {   5,   4},
                {   6,   4},
                { 324, 322},
                { 325, 325},
                {1635,1633},
                {1636,1636},
                {1637,1636},

                // 0 and negative positions:
                {   0,   -2},
                {  -1,   -2},
                {  -2,   -2},
                {  -3,   -5},
                {  -4,   -5},
                {  -5,   -5},
        };
    }

    @DataProvider
    Object[][] provideDataForTestGetAlignedEndPositionOneArg() {
        return new Object[][] {
                {1,  3},
                {2,  3},
                {3,  3},
                {4,  6},
                {5,  6},
                {6,  6},
                {7,  9},
                {8,  9},
                {9,  9},
                {10,12},
                {11,12},
                {12,12},
        };
    }

    @DataProvider
    Object[][] provideDataForGetAlternateSequence() {

        // TODO: Add - strand tests!

        return new Object[][] {
                {
                    "01234567890A1234567890123456789", 12, Allele.create((byte)'A'), Allele.create((byte)'A'), "01234567890A1234567890123456789"
                },
                {
                    "01234567890A1234567890123456789", 11, Allele.create((byte)'A'), Allele.create((byte)'A'), "0123456789AA1234567890123456789"
                },
                {
                    "01234567890A1234567890123456789", 12, Allele.create((byte)'A'), Allele.create("ATGCATGC".getBytes()), "01234567890ATGCATGC1234567890123456789"
                },
                {
                    "AAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCC", 12, Allele.create("GGGGCCC".getBytes()), Allele.create((byte)'T'), "AAAAATTTTTGTCCAAAAATTTTTGGGGGCCCCC"
                },
                {
                    "A", 1, Allele.create((byte)'A'), Allele.create("ATGCATGC".getBytes()), "ATGCATGC"
                },
                {
                    "BA", 2, Allele.create((byte)'A'), Allele.create("ATGCATGC".getBytes()), "BATGCATGC"
                },
                {
                    "AB", 1, Allele.create((byte)'A'), Allele.create("ATGCATGC".getBytes()), "ATGCATGCB"
                },
                {
                    "ATGCATGC", 2, Allele.create((byte)'T'), Allele.create((byte)'G'), "AGGCATGC"
                },
        };
    }

    @DataProvider
    Object[][] provideDataForGetEukaryoticAminoAcidByCodon() {
        return new Object[][] {
                {null, null},
                {"", null},
                {"XQZ", null},
                {"ATG", AminoAcid.METHIONINE},
                {"CCA", AminoAcid.PROLINE},
                {"CCC", AminoAcid.PROLINE},
                {"CCG", AminoAcid.PROLINE},
                {"CCT", AminoAcid.PROLINE},

                // IUPAC base decoding:
                {"GCN", AminoAcid.ALANINE},
                {"CGN", AminoAcid.ARGANINE},
                {"AGR", AminoAcid.ARGANINE},
                {"CGY", AminoAcid.ARGANINE},
                {"MGR", AminoAcid.ARGANINE},
                {"AAY", AminoAcid.ASPARAGINE},
                {"GAY", AminoAcid.ASPARTIC_ACID},
                {"TGY", AminoAcid.CYSTEINE},
                {"CAR", AminoAcid.GLUTAMINE},
                {"GAR", AminoAcid.GLUTAMIC_ACID},
                {"GGN", AminoAcid.GLYCINE},
                {"CAY", AminoAcid.HISTIDINE},
                {"ATH", AminoAcid.ISOLEUCINE},
                {"CTN", AminoAcid.LEUCINE},
                {"TTR", AminoAcid.LEUCINE},
                {"CTY", AminoAcid.LEUCINE},
                {"YTR", AminoAcid.LEUCINE},
                {"AAR", AminoAcid.LYSINE},
                {"TTY", AminoAcid.PHENYLALANINE},
                {"CCN", AminoAcid.PROLINE},
                {"TCN", AminoAcid.SERINE},
                {"AGY", AminoAcid.SERINE},
                {"ACN", AminoAcid.THREONINE},
                {"TAY", AminoAcid.TYROSINE},
                {"GTN", AminoAcid.VALINE},
                {"TRA", AminoAcid.STOP_CODON},
                {"TAR", AminoAcid.STOP_CODON},
        };
    }

    @DataProvider
    Object[][] provideDataForGetMitochondrialAminoAcidByCodon() {
        return new Object[][]{
                // Bad codon values:
                {null,  false, FuncotatorUtils.Genus.UNSPECIFIED, null},
                {"",    false, FuncotatorUtils.Genus.UNSPECIFIED, null},
                {"XQZ", false, FuncotatorUtils.Genus.UNSPECIFIED, null},
                {null,  false, FuncotatorUtils.Genus.UNSPECIFIED, null},
                {"",    false, FuncotatorUtils.Genus.UNSPECIFIED, null},
                {"XQZ", false, FuncotatorUtils.Genus.UNSPECIFIED, null},

                // In sequence (non-starting) codon values:
                {"ATG", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.METHIONINE},
                {"CCA", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.PROLINE},
                {"CCC", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.PROLINE},
                {"CCG", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.PROLINE},
                {"CCT", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.PROLINE},
                {"ATT", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ISOLEUCINE},
                {"ATA", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.METHIONINE},
                {"AGA", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.STOP_CODON},
                {"AGG", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.STOP_CODON},
                {"TGA", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.TRYPTOPHAN},

                // Start codon special cases:
                // Bos:
                {"ATA", true, FuncotatorUtils.Genus.BOS,       AminoAcid.METHIONINE},
                {"ATA", true, FuncotatorUtils.Genus.HOMO,      AminoAcid.METHIONINE},
                {"ATT", true, FuncotatorUtils.Genus.HOMO,      AminoAcid.METHIONINE},
                {"ATT", true, FuncotatorUtils.Genus.MUS,       AminoAcid.METHIONINE},
                {"ATC", true, FuncotatorUtils.Genus.MUS,       AminoAcid.METHIONINE},
                {"ATY", true, FuncotatorUtils.Genus.MUS,       AminoAcid.METHIONINE},
                {"GTG", true, FuncotatorUtils.Genus.CORTURNIX, AminoAcid.METHIONINE},
                {"GTG", true, FuncotatorUtils.Genus.GALLUS,    AminoAcid.METHIONINE},

                // IUPAC base decoding:
                {"GCN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ALANINE},
                {"CGN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ARGANINE},
                {"CGY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ARGANINE},
                {"MGR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ARGANINE},
                {"AAY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ASPARAGINE},
                {"GAY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ASPARTIC_ACID},
                {"TGY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.CYSTEINE},
                {"CAR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.GLUTAMINE},
                {"GAR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.GLUTAMIC_ACID},
                {"GGN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.GLYCINE},
                {"CAY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.HISTIDINE},
                {"ATH", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.ISOLEUCINE},
                {"CTN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.LEUCINE},
                {"TTR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.LEUCINE},
                {"CTY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.LEUCINE},
                {"YTR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.LEUCINE},
                {"AAR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.LYSINE},
                {"TTY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.PHENYLALANINE},
                {"CCN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.PROLINE},
                {"TCN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.SERINE},
                {"AGY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.SERINE},
                {"ACN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.THREONINE},
                {"TAY", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.TYROSINE},
                {"GTN", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.VALINE},
                {"TRA", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.STOP_CODON},
                {"TAR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.STOP_CODON},
                {"AGR", false, FuncotatorUtils.Genus.UNSPECIFIED, AminoAcid.STOP_CODON},
        };
    }

    @DataProvider
    Object[][] provideDataForGetCodonChangeString()  {

//        final String refAllele,
//        final String altAllele,
//        final int alleleStart,
//        final int codingSequenceAlleleStart,
//        final int alignedCodingSequenceAlleleStart,
//        final String alignedCodingSeqRefAllele,
//        final String alignedCodingSeqAltAllele,
//        final alignedAlternateAllele,
//        final int alignedRefAlleleStop,
//        final String contig
//        final Strand strand
//        final Locatable startCodon,
//        final ReferenceSequence codingSequence,
//        final String expected

        // PIK3CA is on chr3 and is transcribed in the FORWARD direction:
        final ReferenceDataSource pik3caTranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
        final String              pik3caFullTranscriptName   = pik3caTranscriptDataSource.getSequenceDictionary().getSequences().stream().filter(s -> s.getSequenceName().startsWith(FuncotatorTestConstants.PIK3CA_TRANSCRIPT)).map(SAMSequenceRecord::getSequenceName).collect(Collectors.joining());
        final ReferenceSequence   pik3caCodingSequence       = pik3caTranscriptDataSource.queryAndPrefetch(pik3caFullTranscriptName, 158, 3364);

        // MUC16 is on chr19 and is transcribed in the REVERSE direction:
        // (The magic numbers here are the coding region for this transscript for MUC16)
        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
        final String              muc16FullTranscriptName   = muc16TranscriptDataSource.getSequenceDictionary().getSequences().stream().filter(s -> s.getSequenceName().startsWith(FuncotatorTestConstants.MUC16_TRANSCRIPT)).map(SAMSequenceRecord::getSequenceName).collect(Collectors.joining());
        final ReferenceSequence   muc16CodingSequence    = muc16TranscriptDataSource.queryAndPrefetch(muc16FullTranscriptName, 205, 43728);

        final SimpleInterval pik3caStart = new SimpleInterval("chr3", 178916614, 178916616);
        final SimpleInterval muc16Start = new SimpleInterval("chr19", 9091812, 9091814);

        return new Object[][] {
                // ONP
                //     alignedCodingSequenceAlleleStart == alignedReferenceAlleleStop   *
                //     alignedCodingSequenceAlleleStart != alignedReferenceAlleleStop
                { "G", "T", 178921515, 997, 997, "GCA", "TCA", "TCA", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-999)Gca>Tca" },
                { "C", "T", 178921516, 998, 997, "GCA", "GTA", "GTA", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-999)gCa>gTa" },
                { "A", "T", 178921517, 999, 997, "GCA", "GCT", "GCT", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-999)gcA>gcT" },
                // INDEL
                //     StartCodonIndel + strand
                { "A", "AT", 178916614, 1, 1, "ATG", "ATTG", "ATTG", 3, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "" },
                //     FrameShift                                                       *
                { "G", "GC", 178921515, 997, 997, "GCA", "GCCA", "GCCA", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-999)gcafs" },
                //     In-Frame Insertion
                //          StartCodon - strand
                { "A", "TGTA", 9091814, 1, 1, "ATG", "TGTATG", "TGTATG", 3, muc16Start.getContig(), Strand.NEGATIVE, muc16CodingSequence, muc16Start, "c.(1-3)atg>TGTatg" },
                //              + Strand
                { "A", "AGCG", 178921515, 999, 997, "GCA", "GCAGCG", "GCAGCG", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(1000-1002)ctc>GCGctc" },
                //              - Strand
                { "T", "TTCG", 9091794, 21, 19, "CTT", "CTCGTT", "CTCGTT", 21, muc16Start.getContig(), Strand.NEGATIVE, muc16CodingSequence, muc16Start, "c.(19-21)ctt>ctCGTt" },
                //          Within a codon:
                //              + Strand
                { "G", "GCAG", 178921515, 997, 997, "GCA", "GCAGCA", "GCAGCA", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-999)gca>gCAGca" },
                //              - Strand
                { "C", "CG", 9091793, 22, 22, "CCT", "CGCT", "CGCT", 24, muc16Start.getContig(), Strand.NEGATIVE, muc16CodingSequence, muc16Start, "c.(22-24)cctfs" },
                //     In-Frame Deletion
                //          Between Codons
                //              + Strand
                { "TGCA", "T", 178921514, 996, 994, "AGTGCA", "AGT", "AGT", 996, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-999)gcadel" },
                //              - Strand
                { "TCCT", "T", 9091794, 21, 19, "CTTCCT", "CTT", "CTT", 21, muc16Start.getContig(), Strand.NEGATIVE, muc16CodingSequence, muc16Start, "c.(19-24)cttcct>ctt" },
                //          Within a codon
                { "GCAC", "G", 178921515, 997, 997, "GCACTC", "GTC", "GTC", 999, pik3caStart.getContig(), Strand.POSITIVE, pik3caCodingSequence, pik3caStart, "c.(997-1002)gcactc>gtc" },
        };
    }

    @DataProvider
    Object[][] provideDataForIsInFrameWithEndOfRegion() {
        return new Object[][] {

                // Starting position checks:
                { 1, 1, false },
                { 1, 2, false },
                { 1, 3, true  },
                { 1, 4, false },
                { 1, 5, false },
                { 1, 6, true  },
                { 1, 7, false },
                { 1, 8, false },
                { 1, 9, true  },
                { 1, 10, false },

                // Middle position checks:
                { 1, 10, false },
                { 2, 10, true },
                { 3, 10, false },
                { 4, 10, false },
                { 5, 10, true },
                { 6, 10, false },
                { 7, 10, false },
                { 8, 10, true },
                { 9, 10, false },
                { 10, 10, false },
                { 56, 473, false },
                { 57, 473, true },
                { 58, 473, false },

                // Last position should always be out of frame.
                { 10, 10, false },
                { 11, 11, false },
                { 12, 12, false },
                { 100, 100, false },
                { 1000, 1000, false },

                // Positions before start of the region.
                // Include these to see if given positions directly flanking the region
                // are in frame.
                { 1, 1, false },
                { 0, 1, false },
                { -1, 1, true  },
                { -2, 1, false  },
                { -3, 1, false  },
                { -4, 1, true  },

                { 1, 2, false },
                { 0, 2, true },
                { -1, 2, false  },
                { -2, 2, false  },
                { -3, 2, true  },
                { -4, 2, false  },
        };
    }
    
    @DataProvider
    Object[][] provideStringDataForGetAminoAcidByLetter() {
        return new Object[][] {
                { "A", AminoAcid.ALANINE },
                { "R", AminoAcid.ARGANINE },
                { "N", AminoAcid.ASPARAGINE },
                { "D", AminoAcid.ASPARTIC_ACID },
                { "C", AminoAcid.CYSTEINE },
                { "E", AminoAcid.GLUTAMIC_ACID },
                { "Q", AminoAcid.GLUTAMINE },
                { "G", AminoAcid.GLYCINE },
                { "H", AminoAcid.HISTIDINE },
                { "I", AminoAcid.ISOLEUCINE },
                { "L", AminoAcid.LEUCINE },
                { "K", AminoAcid.LYSINE },
                { "M", AminoAcid.METHIONINE },
                { "F", AminoAcid.PHENYLALANINE },
                { "P", AminoAcid.PROLINE },
                { "S", AminoAcid.SERINE },
                { "*", AminoAcid.STOP_CODON },
                { "T", AminoAcid.THREONINE },
                { "W", AminoAcid.TRYPTOPHAN },
                { "Y", AminoAcid.TYROSINE },
                { "V", AminoAcid.VALINE },
                { "ahuewifaef", null },
                { "", null },
                { "X", null },
                { "x", null },
                { "7", null },
        };
    }

    @DataProvider
    Object[][] provideCharDataForGetAminoAcidByLetter() {
        return new Object[][] {
                { 'A', AminoAcid.ALANINE },
                { 'R', AminoAcid.ARGANINE },
                { 'N', AminoAcid.ASPARAGINE },
                { 'D', AminoAcid.ASPARTIC_ACID },
                { 'C', AminoAcid.CYSTEINE },
                { 'E', AminoAcid.GLUTAMIC_ACID },
                { 'Q', AminoAcid.GLUTAMINE },
                { 'G', AminoAcid.GLYCINE },
                { 'H', AminoAcid.HISTIDINE },
                { 'I', AminoAcid.ISOLEUCINE },
                { 'L', AminoAcid.LEUCINE },
                { 'K', AminoAcid.LYSINE },
                { 'M', AminoAcid.METHIONINE },
                { 'F', AminoAcid.PHENYLALANINE },
                { 'P', AminoAcid.PROLINE },
                { 'S', AminoAcid.SERINE },
                { '*', AminoAcid.STOP_CODON },
                { 'T', AminoAcid.THREONINE },
                { 'W', AminoAcid.TRYPTOPHAN },
                { 'Y', AminoAcid.TYROSINE },
                { 'V', AminoAcid.VALINE },
                { 'a', null },
                { '\0', null },
                { 'X', null },
                { 'x', null },
                { '7', null },
        };
    }

    @DataProvider
    Object[][] provideDataForGetProteinChangePosition() {
        return new Object[][] {
                { 1, 1 },
                { 2, 1 },
                { 3, 1 },
                { 4, 2 },
                { 5, 2 },
                { 6, 2 },
                { 300, 100 },
                { 301, 101 },
                { 302, 101 },
                { 303, 101 },
                { 304, 102 },
        };
    }

    @DataProvider
    Object[][] provideDataForRenderProteinChangeString() {

        final SimpleInterval startCodon = new SimpleInterval("TEST", 1,3);

        // NOTE: The ref and alt alleles only matter for length - content is not checked.

        // Arrays have the following fields:
        // protChangeStartPos,protChangeEndPos,refAminoAcidSeq,altAminoAcidSeq,refAllele,altAllele,startCodon,refAlleleStartPosition,expected

        return new Object[][] {
                // SNPs:
                //    No protein change:
                { 1,     1, "N", "N", "C", "A", startCodon, 1, "p.N1N" },
                //    With protein change:
                { 1,     1, "N", "K", "C", "A", null,       1, "p.N1K" },
                { 1,     1, "N", "K", "C", "A", startCodon, 1, "p.N1K" },
                { 5,     5, "F", "L", "C", "A", startCodon, 1, "p.F5L" },
                { 10,   10, "A", "V", "C", "G", startCodon, 1, "p.A10V" },
                { 100, 100, "Y", "P", "A", "T", startCodon, 1, "p.Y100P" },

                // DNPs:
                //    No protein change:
                { 2,     2, "K",  "K",  "CG", "AT", null,     1, "p.K2K" },
                { 1,     2, "NK", "NK", "CG", "AT", null,     1, "p.1_2NK>NK" },
                //    1 protein change:
                { 1,     1, "N", "G", "CG", "AT", null,       1, "p.N1G" },
                { 1,     1, "K", "T", "CG", "AT", startCodon, 1, "p.K1T" },
                { 5,     5, "V", "I", "TT", "CG", startCodon, 1, "p.V5I" },
                { 10,   10, "R", "*", "AA", "GT", startCodon, 1, "p.R10*" },
                { 100, 100, "*", "L", "AG", "TC", startCodon, 1, "p.*100L" },
                //    2 protein changes:
                { 1,     2, "NK", "TG", "CG", "AT", null,       1, "p.1_2NK>TG" },
                { 1,     2, "NK", "TG", "CG", "AT", startCodon, 1, "p.1_2NK>TG" },
                { 5,     6, "FV", "LI", "TT", "CG", startCodon, 1, "p.5_6FV>LI" },
                { 10,   11, "RS", "W*", "AA", "GT", startCodon, 1, "p.10_11RS>W*" },
                { 100, 101, "Q*", "FL", "AG", "TC", startCodon, 1, "p.100_101Q*>FL" },

                // TNPs:
                //    No protein change:
                { 2,     2, "R",  "R",  "CGA", "AGG", null,       1, "p.R2R" },
                //    1 protein change:
                { 1,     1, "I", "G", "ATT", "GGA", null,         1, "p.I1G" },
                { 5,     5, "V", "T", "GTT", "ACC", startCodon,   1, "p.V5T" },
                //    2 protein changes:
                { 101, 102, "LS", "FR", "ATC", "TCG", null,       1, "p.101_102LS>FR" },
                { 201, 202, "VK", "GR", "TGA", "GGC", startCodon, 1, "p.201_202VK>GR" },

                // MNPs:
                //    No protein change:
                { 2,     5, "RQV*",  "RQV*",  "CGCCAGGTGTAA", "CGTCAAGTTTGA", null, 1, "p.2_5RQV*>RQV*" },
                { 10,   32, "RACDEFGHIKLMNPQRSTVWYR",
                            "RACDEFGHIKLMNPQRSTVWYR",
                            "CGCGCGTGCGATGAATTTGGCCATATTAAACTGATGAACCCGCAGCGCAGCACCGTGTGGTATCGC",
                            "AGAGCTTGTGACGAGTTCGGACACATCAAGTTCATGAATCCCCAAAGAAGAACAGTATGGTACAGA",
                            startCodon, 1, "p.10_32RACDEFGHIKLMNPQRSTVWYR>RACDEFGHIKLMNPQRSTVWYR" },
                //    2 protein changes:
                { 17,   18, "F*", "LT", "CTAG", "GACT", null, 1, "p.17_18F*>LT" },
                //    10 protein changes:
                { 51,   60, "ACDEFGHIKL", "MNPQRSTVWY",
                        "GCGTGCGATGAATTTGGCCATATTAAACTG", "ATGAACCCGCAGCGCAGCACCGTGTGGTAT",
                        null, 1, "p.51_60ACDEFGHIKL>MNPQRSTVWY" },

                // FS INSs:
                //    Start codon overlap:
                { 1,     1, "M", "M", "A", "ATG",  startCodon, 1, "" },
                { 1,     1, "M", "M", "A", "AT",   startCodon, 1, "" },
                { 1,     1, "M", "M", "T", "TG",   startCodon, 2, "" },
                { 1,     1, "M", "M", "T", "TGC",  startCodon, 2, "" },
                //    No start codon overlap:
                { 2,     2, "K", "K", "T", "TG",   startCodon, 3,  "p.K2fs" },
                { 4,     4, "G", "G", "C", "CA",   startCodon, 10, "p.G4fs" },
                { 4,     4, "G", "G", "C", "CAC",  startCodon, 10, "p.G4fs" },
                { 4,     4, "Q", "Q", "A", "ATG",  startCodon, 10, "p.Q4fs" },
                { 4,     4, "V", "V", "T", "TG",   startCodon, 10, "p.V4fs" },

                // FS DELs:
                //    Start codon overlap:
                { 1,     1, "M", "M", "ATG", "AT", startCodon, 1, "" },
                { 1,     1, "M", "M", "AT", "A",   startCodon, 1, "" },
                { 1,     1, "M", "M", "TG", "T",   startCodon, 2, "" },
                { 1,     1, "M", "M", "CA", "C",   new SimpleInterval("TEST", 100, 102), 99, "" },
                //    No start codon overlap:
                {  4,    4, "M", "M", "ATG", "AT",    startCodon, 10, "p.M4fs" },
                {  5,    5, "V", "V", "AT", "A",      startCodon, 12, "p.V5fs" },
                { 10,   10, "W", "W", "TG", "T",      startCodon, 30, "p.W10fs" },
                {  9,    9, "GS", "GS", "CTGGG", "C", startCodon, 27, "p.GS9fs" },

                // IF INSs:
                //    Between codons:
                //       3 bases:
                { 544, 545, "", "G", "T", "TGGA", startCodon, 178936090, "p.544_545insG" },
                { 545, 546, "", "W", "G", "GTGG", startCodon, 178936093, "p.545_546insW" },
                //       6 bases:
                { 544, 545, "", "GC", "T", "TGGACTT", startCodon, 178936090, "p.544_545insGC" },
                { 545, 546, "", "W*", "G", "GTGGTGA", startCodon, 178936093, "p.545_546insW*" },
                //       9 bases:
                { 544, 545, "", "GCA", "T", "TGGACTTATT", startCodon, 178936090, "p.544_545insGCA" },
                { 545, 546, "", "WV*", "G", "GTGGGTATGA", startCodon, 178936093, "p.545_546insWV*" },
                //    Within a Codon:
                //       3 bases:
                { 544, 545, "", "A", "A", "ACTG", startCodon, 178936088, "p.544_545insA" },
                { 544, 545, "", "S", "C", "CCAG", startCodon, 178936089, "p.544_545insS" },
                { 544, 545, "T", "IA", "A", "ATCG", startCodon, 178936088, "p.544_545T>IA" },
                //       6 bases:
                { 545, 545, "E", "DV*", "G", "GATGTCT", startCodon, 178936091, "p.E545DV*" },
                { 545, 545, "E", "DLQ", "A", "ACCTGCA", startCodon, 178936092, "p.E545DLQ" },
                { 545, 546, "", "QE", "A", "AGCAGGA", startCodon, 178936092, "p.545_546insQE" },

                // IF DELs:
                //    Complete codon deletions:
                //       3 bases:
                { 544, 544, "T", "", "CACT", "C", startCodon, 1631, "p.T544del" },
                { 545, 545, "E", "", "TGAG", "T", startCodon, 1634, "p.E545del" },
                //       6 bases:
                { 544, 544, "TE", "", "CACTGAG", "C", startCodon, 1631, "p.TE544del" },
                { 545, 545, "EQ", "", "TGAGCAG", "T", startCodon, 1634, "p.EQ545del" },
                //       9 bases:
                { 544, 544, "TEQ", "", "CACTGAGCAG", "C", startCodon, 1631, "p.TEQ544del" },
                { 545, 545, "EQE", "", "TGAGCAGGAG", "T", startCodon, 1634, "p.EQE545del" },
                //    Within a Codon:
                //       3 bases:
                { 544, 545, "TE", "K", "ACTG", "A", startCodon, 1632, "p.544_545TE>K" },
                { 545, 545, "E",  "E", "CTGA", "C", startCodon, 1633, "p.E545del" },
                //       6 bases:
                { 544, 546, "TEQ", "K", "ACTGAGC", "A", startCodon, 1632, "p.544_546TEQ>K" },
                { 545, 545, "EQ", "", "CTGAGCA", "C", startCodon, 1633, "p.EQ545del" },
                //       9 bases:
                { 545, 545, "EQE", "", "ACTGAGCAGG", "A", startCodon, 1632, "p.EQE545del" },
                { 545, 545, "EQE", "", "CTGAGCAGGA", "C", startCodon, 1633, "p.EQE545del" },
        };
    }

    @DataProvider
    Object[][] provideForGetClosestExonIndex() {

        final List<Locatable> exonList = new ArrayList<>();
        exonList.add( new SimpleInterval("testContig", 20, 30) );
        exonList.add( new SimpleInterval("testContig", 40, 50) );
        exonList.add( new SimpleInterval("testContig", 60, 70) );
        exonList.add( new SimpleInterval("testContig", 80, 90) );
        exonList.add( new SimpleInterval("testContig", 10000, 15000) );
        exonList.add( new SimpleInterval("testContig", 10000000, 15000000) );

        return new Object[][] {
                // No exons:
                {  8, Collections.emptyList(), -1 },
                // Start before first exon start:
                {  8, exonList, 0 },
                // Start inside first exon:
                { 25, exonList, 0 },
                // Start after first exon end:
                { 31, exonList, 0 },
                // Start after exon start, just before a "break point":
                { 35, exonList, 0 },
                // Start between exons, just after a "break point":
                { 36, exonList, 1 },
                // Start between exons:
                { 500, exonList, 3 },
                // Start after last exon:
                { 50000000, exonList, 5 },
        };
    }

    @DataProvider
    Object[][] provideForCreateIntronicCDnaString() {

        final List<Locatable> exonList = new ArrayList<>();
        exonList.add( new SimpleInterval("testContig", 20, 30) );
        exonList.add( new SimpleInterval("testContig", 40, 50) );
        exonList.add( new SimpleInterval("testContig", 60, 70) );
        exonList.add( new SimpleInterval("testContig", 80, 90) );
        exonList.add( new SimpleInterval("testContig", 10000, 15000) );
        exonList.add( new SimpleInterval("testContig", 10000000, 15000000) );

        return new Object[][] {
                // No exons:
                helpProvideForCreateIntronicCDnaString(  8, Collections.emptyList(), "A", "C", null, null ),
                // Start before first exon start:
                helpProvideForCreateIntronicCDnaString(  8, exonList,"A", "T", 0, -12 ),
                // Start inside first exon:
                helpProvideForCreateIntronicCDnaString( 25, exonList,"C", "T", 0, 5 ),
                // Start after first exon end:
                helpProvideForCreateIntronicCDnaString( 31, exonList, "G", "T", 0, 1 ),
                // Start after exon start, just before a "break point":
                helpProvideForCreateIntronicCDnaString( 35, exonList, "ATG", "GCT", 0, 5),
                // Start between exons, just after a "break point":
                helpProvideForCreateIntronicCDnaString( 36, exonList, "ATG", "GCT", 1, -4),
                // Start between exons:
                helpProvideForCreateIntronicCDnaString( 500, exonList, "TTTTTTTTT", "AAAAAAAAA", 3, 410),
                // Start after last exon:
                helpProvideForCreateIntronicCDnaString( 50000000, exonList, "A", "GTG",5, 35000000 ),

        };
    }

    @DataProvider
    Object[][] provideForTestGetCodingSequenceChangeString() {

//        final int codingSequenceAlleleStart,
//        final String referenceAllele,
//        final String alternateAllele,
//        final Strand strand,
//        final String expected

        // + Strand tests from PIK3CA
        // - Strand tests from MUC16

        // TODO: Add tests for the deletion case requiring exon start/end/allele position.

        return new Object[][] {
                // SNP
                { 1632, "A", "T", Strand.POSITIVE, "c.1632A>T" },
                // ONP (|ONP| > 1)
                { 1633, "CT", "GG", Strand.POSITIVE, "c.1633_1634CT>GG" },
                // Insertion (+ Strand)
                { 1634, "T", "TGG", Strand.POSITIVE, "c.1634_1635insGG" },
                // Insertion (- Strand)
                { 43303, "A", "TCA", Strand.NEGATIVE, "c.43302_43303insTC" },
                // Deletion (+ Strand)
                { 1634, "TGAG", "T", Strand.POSITIVE, "c.1635_1637delGAG" },
                // Deletion (- Strand)
                { 43138, "TCT", "T", Strand.NEGATIVE, "c.43138_43139delTC" },
        };
    }

    @DataProvider
    Object[][] provideForTestIsIndelBetweenCodons() {

        // Arrays have the following fields:
        // codingSequenceAlleleStart, alignedCodingSequenceAlleleStart, refAllele, strand, expected

        return new Object[][] {
                // + Strand:
                { 4,  4, "A",   Strand.POSITIVE, false },
                { 5,  4, "A",   Strand.POSITIVE, false },
                { 6,  4, "A",   Strand.POSITIVE, true },
                { 4,  4, "AC",  Strand.POSITIVE, false },
                { 5,  4, "AC",  Strand.POSITIVE, true },
                { 6,  4, "AC",  Strand.POSITIVE, false },
                { 4,  4, "ACT", Strand.POSITIVE, true },
                { 5,  4, "ACT", Strand.POSITIVE, false },
                { 6,  4, "ACT", Strand.POSITIVE, false },

                // - Strand:
                { 4,  4, "A",   Strand.NEGATIVE, true },
                { 5,  4, "A",   Strand.NEGATIVE, false },
                { 6,  4, "A",   Strand.NEGATIVE, false },
                { 4,  4, "AC",  Strand.NEGATIVE, true },
                { 5,  4, "AC",  Strand.NEGATIVE, false },
                { 6,  4, "AC",  Strand.NEGATIVE, false },
                { 4,  4, "ACT", Strand.NEGATIVE, true },
                { 5,  4, "ACT", Strand.NEGATIVE, false },
                { 6,  4, "ACT", Strand.NEGATIVE, false },
        };
    }

    @DataProvider
    Object[][] provideDataForGetProteinChangeEndPosition() {
        return new Object[][] {
                {1,  3, 1},
                {1,  6, 2},
                {1,  9, 3},
                {1, 12, 4},
                {1, 15, 5},
                {1, 18, 6},
                {1, 21, 7},
                {1, 24, 8},
        };
    }

    @DataProvider
    Object[][] provideDataForGetAlignedCodingSequenceAllele() {

        final String seq = "ATGAAAGGGGTGCCTATGCTAGATAGACAGATAGTGTGTGTGTGTGTGCGCGCGCGCGCGCGTTGTTAG";

        //CTA ACA ACG CGC GCG CGC GCG CAC ACA CAC ACA CAC TAT CTG TCT ATC TAG CAT AGG CAC CCC TTT CAT

//        final Allele refAllele,
//        final Integer refAlleleStart,

        return new Object[][] {
                { seq,  1, 3,  Allele.create("ATG", true), 1, Strand.POSITIVE, "ATG" },
                { seq,  4, 6,  Allele.create("AAA", true), 4, Strand.POSITIVE, "AAA" },
                { seq,  7, 9,  Allele.create("GGG", true), 7, Strand.POSITIVE, "GGG" },
                { seq, 10, 12, Allele.create("GTG", true), 10, Strand.POSITIVE, "GTG" },
                { seq, 13, 15, Allele.create("CCT", true), 13, Strand.POSITIVE, "CCT" },
                { seq, 16, 18, Allele.create("ATG", true), 16, Strand.POSITIVE, "ATG" },
                { seq, 19, 21, Allele.create("CTA", true), 19, Strand.POSITIVE, "CTA" },
                { seq,  1,  6, Allele.create("ATGAAA", true), 1, Strand.POSITIVE, "ATGAAA" },
                { seq,  4,  9, Allele.create("AAAGGG", true), 4, Strand.POSITIVE, "AAAGGG" },
                { seq,  7, 12, Allele.create("GGGGTG", true), 7, Strand.POSITIVE, "GGGGTG" },
                { seq, 10, 15, Allele.create("GTGCCT", true), 10, Strand.POSITIVE, "GTGCCT" },
                { seq, 13, 18, Allele.create("CCTATG", true), 13, Strand.POSITIVE, "CCTATG" },
                { seq, 16, 21, Allele.create("ATGCTA", true), 16, Strand.POSITIVE, "ATGCTA" },
                { seq, 19, 24, Allele.create("CTAGAT", true), 19, Strand.POSITIVE, "CTAGAT" },
                { seq, 1, seq.length(), Allele.create(seq, true), 1, Strand.POSITIVE, seq },

                { seq,  1, 3,  Allele.create("CTA", true), 1, Strand.NEGATIVE, "CTA" },
                { seq,  4, 6,  Allele.create("ACA", true), 4, Strand.NEGATIVE, "ACA" },
                { seq,  7, 9,  Allele.create("ACG", true), 7, Strand.NEGATIVE, "ACG" },
                { seq, 10, 12, Allele.create("CGC", true), 10, Strand.NEGATIVE, "CGC" },
                { seq, 13, 15, Allele.create("GCG", true), 13, Strand.NEGATIVE, "GCG" },
                { seq, 16, 18, Allele.create("CGC", true), 16, Strand.NEGATIVE, "CGC" },
                { seq, 19, 21, Allele.create("GCG", true), 19, Strand.NEGATIVE, "GCG" },
                { seq,  1,  6, Allele.create("CTAACA", true), 1, Strand.NEGATIVE, "CTAACA" },
                { seq,  4,  9, Allele.create("ACAACG", true), 4, Strand.NEGATIVE, "ACAACG" },
                { seq,  7, 12, Allele.create("ACGCGC", true), 7, Strand.NEGATIVE, "ACGCGC" },
                { seq, 10, 15, Allele.create("CGCGCG", true), 10, Strand.NEGATIVE, "CGCGCG" },
                { seq, 13, 18, Allele.create("GCGCGC", true), 13, Strand.NEGATIVE, "GCGCGC" },
                { seq, 16, 21, Allele.create("CGCGCG", true), 16, Strand.NEGATIVE, "CGCGCG" },
                { seq, 19, 24, Allele.create("GCGCAC", true), 19, Strand.NEGATIVE, "GCGCAC" },
                { seq, 1, seq.length(), Allele.create(ReadUtils.getBasesReverseComplement( seq.getBytes() ), true), 1, Strand.NEGATIVE, ReadUtils.getBasesReverseComplement( seq.getBytes() ) },
        };
    }

    @DataProvider
    Object[][] provideDataForTestGetAlignedRefAllele() {

//        final String referenceSnippet,
//        final int referencePadding,
//        final Allele refAllele,
//        final int codingSequenceRefAlleleStart,
//        final int alignedRefAlleleStart
//        final Strand strand
//        expected

        // TODO: Make tests with Strand.NEGATIVE!!!!  (issue #5351 - https://github.com/broadinstitute/gatk/issues/5351)
        // TODO: Make tests with alt allele longer/shorter than ref allele!  (issue #5351 - https://github.com/broadinstitute/gatk/issues/5351)

//                                                 11111111112222222222333333333344444444445555555555666666
//                                       012345678901234567890123456789012345678901234567890123456789012345
        final String referenceSnippet = "AAATTTGGGCCCATGATATAGGCGCCGTAGCAGTAGATAGCCCCCCAACCGGGGCCCGGGTTTAAA";

        return new Object[][] {

                // alignedRefAlleleStart <= 0
                {referenceSnippet, 10, Allele.create("CC", true), 1, -1, Strand.POSITIVE, "GCCCAT"},
                {referenceSnippet, 11, Allele.create("CA", true), 1, -1, Strand.POSITIVE, "CCCATG"},
                {referenceSnippet, 12, Allele.create("AT", true), 1, -1, Strand.POSITIVE, "CCATGA"},

                // (codingSequenceRefAlleleStart <= 0) && (alignedRefAlleleStart <= 0)
                {referenceSnippet, 10, Allele.create("CC", true), 0, -2, Strand.POSITIVE, "GCCCAT"},
                {referenceSnippet, 11, Allele.create("CA", true), 0, -2, Strand.POSITIVE, "CCCATG"},
                {referenceSnippet, 12, Allele.create("AT", true), 0, -2, Strand.POSITIVE, "CCATGA"},

                {referenceSnippet, 10, Allele.create("CCA", true),  -1, -2, Strand.POSITIVE, "CCCATG"},
                {referenceSnippet, 11, Allele.create("CAT", true),  -1, -2, Strand.POSITIVE, "CCATGA"},
                {referenceSnippet, 12, Allele.create("ATG", true),  -1, -2, Strand.POSITIVE, "CATGAT"},

                // Allele same as reference sequence:
                {referenceSnippet, 0, Allele.create("AAA", true), 1, 1, Strand.POSITIVE, "AAA"},
                {referenceSnippet, 1, Allele.create("AA", true), 2, 1, Strand.POSITIVE, "AAA"},
                {referenceSnippet, 2, Allele.create("A", true), 3, 1, Strand.POSITIVE, "AAA"},

                {referenceSnippet, 6, Allele.create("GGG", true), 1, 1, Strand.POSITIVE, "GGG"},
                {referenceSnippet, 7, Allele.create("GG", true), 2, 1, Strand.POSITIVE, "GGG"},
                {referenceSnippet, 8, Allele.create("G", true), 3, 1, Strand.POSITIVE, "GGG"},

                {referenceSnippet, 6, Allele.create("GGG", true), 60, 60, Strand.POSITIVE, "GGG"},
                {referenceSnippet, 7, Allele.create("GG", true), 61, 60, Strand.POSITIVE, "GGG"},
                {referenceSnippet, 8, Allele.create("G", true), 62, 60, Strand.POSITIVE, "GGG"},

                {referenceSnippet, 17, Allele.create("ATAG", true), 18, 17, Strand.POSITIVE, "TATAGG"},

                {referenceSnippet, 4, Allele.create(referenceSnippet.substring(4), true), 18, 17, Strand.POSITIVE, referenceSnippet.substring(3)},

                // Allele diffferent from reference sequence:
                {referenceSnippet, 0, Allele.create("TAA", true), 1, 1, Strand.POSITIVE, "TAA"},
                {referenceSnippet, 1, Allele.create("AG", true), 2, 1, Strand.POSITIVE, "AAG"},
                {referenceSnippet, 2, Allele.create("T", true), 3, 1, Strand.POSITIVE, "AAT"},

                {referenceSnippet, 6, Allele.create("GCC", true), 1, 1, Strand.POSITIVE, "GCC"},
                {referenceSnippet, 7, Allele.create("AG", true), 2, 1, Strand.POSITIVE, "GAG"},
                {referenceSnippet, 8, Allele.create("A", true), 3, 1, Strand.POSITIVE, "GGA"},

                {referenceSnippet, 6, Allele.create("GCC", true), 60, 60, Strand.POSITIVE, "GCC"},
                {referenceSnippet, 7, Allele.create("AG", true), 61, 60, Strand.POSITIVE, "GAG"},
                {referenceSnippet, 8, Allele.create("A", true), 62, 60, Strand.POSITIVE, "GGA"},

                {referenceSnippet, 17, Allele.create("ATAT", true), 18, 17, Strand.POSITIVE, "TATATG"},

                {referenceSnippet, 4, Allele.create(
                        new String(new char[referenceSnippet.length() - 4]).replace("\0", "A"), true),
                        18, 17, Strand.POSITIVE,
                        "T" + new String(new char[referenceSnippet.length() - 4]).replace("\0", "A")},
        };
    }

    @DataProvider
    Object[][] provideDataForGetTranscriptAlleleStartPosition() {

        final VariantContext variant =
                new VariantContextBuilder(
                        "",
                        "1",
                        100,
                        100,
                        Arrays.asList( Allele.create("A", true), Allele.create("T") )
                ).make();

        return new Object[][] {
                {
                    variant,
                    Collections.singletonList( new SimpleInterval("1", 100, 100) ),
                    Strand.POSITIVE,
                    1
                },
                {
                    variant,
                    Collections.singletonList( new SimpleInterval("1", 100, 100) ),
                    Strand.NEGATIVE,
                    1
                },
                {
                    variant,
                    Collections.singletonList( new SimpleInterval("1", 50, 200) ),
                    Strand.POSITIVE,
                    51
                },
                {
                    variant,
                    Collections.singletonList( new SimpleInterval("1", 50, 200) ),
                    Strand.NEGATIVE,
                    101
                },
                {
                    variant,
                    Arrays.asList( new SimpleInterval("1", 50, 100), new SimpleInterval("1", 101, 200) ),
                    Strand.POSITIVE,
                    51
                },
                {
                    variant,
                    Arrays.asList( new SimpleInterval("1", 50, 100), new SimpleInterval("1", 101, 200) ),
                    Strand.NEGATIVE,
                    101
                },
                {
                    variant,
                    Arrays.asList(
                            new SimpleInterval("1",   1, 10),
                            new SimpleInterval("1",  11, 20),
                            new SimpleInterval("1",  21, 30),
                            new SimpleInterval("1",  31, 40),
                            new SimpleInterval("1",  41, 50),
                            new SimpleInterval("1",  51, 60),
                            new SimpleInterval("1",  61, 70),
                            new SimpleInterval("1",  71, 80),
                            new SimpleInterval("1",  81, 90),
                            new SimpleInterval("1",  91, 100),
                            new SimpleInterval("1", 101, 110),
                            new SimpleInterval("1", 111, 120),
                            new SimpleInterval("1", 121, 130),
                            new SimpleInterval("1", 131, 140),
                            new SimpleInterval("1", 141, 150),
                            new SimpleInterval("1", 151, 160),
                            new SimpleInterval("1", 161, 170),
                            new SimpleInterval("1", 171, 180),
                            new SimpleInterval("1", 181, 190),
                            new SimpleInterval("1", 191, 200)
                    ),
                    Strand.POSITIVE,
                    100
                },
                {
                    variant,
                        Arrays.asList(
                                new SimpleInterval("1",   1, 10),
                                new SimpleInterval("1",  11, 20),
                                new SimpleInterval("1",  21, 30),
                                new SimpleInterval("1",  31, 40),
                                new SimpleInterval("1",  41, 50),
                                new SimpleInterval("1",  51, 60),
                                new SimpleInterval("1",  61, 70),
                                new SimpleInterval("1",  71, 80),
                                new SimpleInterval("1",  81, 90),
                                new SimpleInterval("1",  91, 100),
                                new SimpleInterval("1", 101, 110),
                                new SimpleInterval("1", 111, 120),
                                new SimpleInterval("1", 121, 130),
                                new SimpleInterval("1", 131, 140),
                                new SimpleInterval("1", 141, 150),
                                new SimpleInterval("1", 151, 160),
                                new SimpleInterval("1", 161, 170),
                                new SimpleInterval("1", 171, 180),
                                new SimpleInterval("1", 181, 190),
                                new SimpleInterval("1", 191, 200)
                        ),
                    Strand.NEGATIVE,
                    101
                },
        };
    }

    @DataProvider
    Object[][] provideDataForTestCreateSpliceSiteCodonChange() {

        return new Object[][] {
                {1000, 5, 1000, 1500, Strand.POSITIVE, 0, "c.e5-0"},
                {1000, 4, 1, 1500, Strand.POSITIVE,    0, "c.e4+500"},
                {1000, 3, 500, 1500, Strand.POSITIVE,  0, "c.e3-500"},

                {1000, 5, 1000, 1500, Strand.NEGATIVE, 0, "c.e5+0"},
                {1000, 4, 1, 1500, Strand.NEGATIVE,    0, "c.e4-500"},
                {1000, 3, 500, 1500, Strand.NEGATIVE,  0, "c.e3+500"},

                {1000, 5, 1500, 500, Strand.NEGATIVE,  0, "c.e5+500"},

                {1000, 5, 1000, 1500, Strand.POSITIVE, 1, "c.e5+1"},
                {1000, 4, 1, 1500, Strand.POSITIVE,    2, "c.e4+502"},
                {1000, 3, 500, 1500, Strand.POSITIVE,  3, "c.e3-497"},

                {1000, 5, 1000, 1500, Strand.NEGATIVE, 4, "c.e5+4"},
                {1000, 4, 1, 1500, Strand.NEGATIVE,    5, "c.e4-495"},
                {1000, 3, 500, 1500, Strand.NEGATIVE,  6, "c.e3+506"},

                {1000, 5, 1500, 500, Strand.NEGATIVE,  7, "c.e5+507"},

                {1000, 5, 1000, 1500, Strand.POSITIVE, -1, "c.e5-1"},
                {1000, 4, 1, 1500, Strand.POSITIVE,    -2, "c.e4+498"},
                {1000, 3, 500, 1500, Strand.POSITIVE,  -3, "c.e3-503"},

                {1000, 5, 1000, 1500, Strand.NEGATIVE, -4, "c.e5-4"},
                {1000, 4, 1, 1500, Strand.NEGATIVE,    -5, "c.e4-505"},
                {1000, 3, 500, 1500, Strand.NEGATIVE,  -6, "c.e3+494"},

                {1000, 5, 1500, 500, Strand.NEGATIVE,  -7, "c.e5+493"},
        };
    }

    @DataProvider
    Object[][] provideDataForTestGetOverlappingExonPositions() {
//        refAllele, altAllele, contig, start, stop, strand, exonPositionList

        final String contig = "chr3";

        final List<SimpleInterval> exonPositionList = new ArrayList<>(10);

        final SimpleInterval interval1 = new SimpleInterval(contig, 100, 199);
        final SimpleInterval interval2 = new SimpleInterval(contig, 300, 399);
        final SimpleInterval interval3 = new SimpleInterval(contig, 500, 599);
        final SimpleInterval interval4 = new SimpleInterval(contig, 700, 799);
        final SimpleInterval interval5 = new SimpleInterval(contig, 900, 999);

        exonPositionList.add( interval1 );
        exonPositionList.add( interval2 );
        exonPositionList.add( interval3 );
        exonPositionList.add( interval4 );
        exonPositionList.add( interval5 );

        return new Object[][] {
                // + Strand:
                { Allele.create("A", true),          Allele.create("T"),          contig,  50,  50, Strand.POSITIVE, exonPositionList, null },
                { Allele.create("A", true),          Allele.create("T"),          contig, 299, 299, Strand.POSITIVE, exonPositionList, null },
                { Allele.create("A", true),          Allele.create("T"),          contig, 400, 400, Strand.POSITIVE, exonPositionList, null },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig,  40,  49, Strand.POSITIVE, exonPositionList, null },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 290, 299, Strand.POSITIVE, exonPositionList, null },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 400, 409, Strand.POSITIVE, exonPositionList, null },

                { Allele.create("A", true),          Allele.create("T"),          contig, 100, 100, Strand.POSITIVE, exonPositionList, interval1 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 150, 150, Strand.POSITIVE, exonPositionList, interval1 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 199, 199, Strand.POSITIVE, exonPositionList, interval1 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig,  91, 100, Strand.POSITIVE, exonPositionList, interval1 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 120, 129, Strand.POSITIVE, exonPositionList, interval1 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 199, 208, Strand.POSITIVE, exonPositionList, interval1 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 300, 300, Strand.POSITIVE, exonPositionList, interval2 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 350, 350, Strand.POSITIVE, exonPositionList, interval2 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 399, 399, Strand.POSITIVE, exonPositionList, interval2 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 291, 300, Strand.POSITIVE, exonPositionList, interval2 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 320, 329, Strand.POSITIVE, exonPositionList, interval2 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 399, 408, Strand.POSITIVE, exonPositionList, interval2 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 500, 500, Strand.POSITIVE, exonPositionList, interval3 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 550, 550, Strand.POSITIVE, exonPositionList, interval3 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 599, 599, Strand.POSITIVE, exonPositionList, interval3 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 491, 500, Strand.POSITIVE, exonPositionList, interval3 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 520, 529, Strand.POSITIVE, exonPositionList, interval3 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 599, 608, Strand.POSITIVE, exonPositionList, interval3 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 700, 700, Strand.POSITIVE, exonPositionList, interval4 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 750, 750, Strand.POSITIVE, exonPositionList, interval4 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 799, 799, Strand.POSITIVE, exonPositionList, interval4 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 691, 700, Strand.POSITIVE, exonPositionList, interval4 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 720, 729, Strand.POSITIVE, exonPositionList, interval4 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 799, 808, Strand.POSITIVE, exonPositionList, interval4 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 900,  900, Strand.POSITIVE, exonPositionList, interval5 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 950,  950, Strand.POSITIVE, exonPositionList, interval5 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 999,  999, Strand.POSITIVE, exonPositionList, interval5 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 891,  900, Strand.POSITIVE, exonPositionList, interval5 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 920,  929, Strand.POSITIVE, exonPositionList, interval5 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 999, 1008, Strand.POSITIVE, exonPositionList, interval5 },
                
                // - Strand:
                { Allele.create("A", true),          Allele.create("T"),          contig,  50,  50, Strand.NEGATIVE, exonPositionList, null },
                { Allele.create("A", true),          Allele.create("T"),          contig, 299, 299, Strand.NEGATIVE, exonPositionList, null },
                { Allele.create("A", true),          Allele.create("T"),          contig, 400, 400, Strand.NEGATIVE, exonPositionList, null },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig,  40,  49, Strand.NEGATIVE, exonPositionList, null },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 290, 299, Strand.NEGATIVE, exonPositionList, null },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 400, 409, Strand.NEGATIVE, exonPositionList, null },

                { Allele.create("A", true),          Allele.create("T"),          contig, 100, 100, Strand.NEGATIVE, exonPositionList, interval1 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 150, 150, Strand.NEGATIVE, exonPositionList, interval1 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 199, 199, Strand.NEGATIVE, exonPositionList, interval1 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig,  91, 100, Strand.NEGATIVE, exonPositionList, interval1 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 120, 129, Strand.NEGATIVE, exonPositionList, interval1 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 199, 208, Strand.NEGATIVE, exonPositionList, interval1 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 300, 300, Strand.NEGATIVE, exonPositionList, interval2 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 350, 350, Strand.NEGATIVE, exonPositionList, interval2 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 399, 399, Strand.NEGATIVE, exonPositionList, interval2 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 291, 300, Strand.NEGATIVE, exonPositionList, interval2 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 320, 329, Strand.NEGATIVE, exonPositionList, interval2 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 399, 408, Strand.NEGATIVE, exonPositionList, interval2 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 500, 500, Strand.NEGATIVE, exonPositionList, interval3 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 550, 550, Strand.NEGATIVE, exonPositionList, interval3 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 599, 599, Strand.NEGATIVE, exonPositionList, interval3 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 491, 500, Strand.NEGATIVE, exonPositionList, interval3 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 520, 529, Strand.NEGATIVE, exonPositionList, interval3 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 599, 608, Strand.NEGATIVE, exonPositionList, interval3 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 700, 700, Strand.NEGATIVE, exonPositionList, interval4 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 750, 750, Strand.NEGATIVE, exonPositionList, interval4 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 799, 799, Strand.NEGATIVE, exonPositionList, interval4 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 691, 700, Strand.NEGATIVE, exonPositionList, interval4 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 720, 729, Strand.NEGATIVE, exonPositionList, interval4 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 799, 808, Strand.NEGATIVE, exonPositionList, interval4 },

                { Allele.create("A", true),          Allele.create("T"),          contig, 900,  900, Strand.NEGATIVE, exonPositionList, interval5 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 950,  950, Strand.NEGATIVE, exonPositionList, interval5 },
                { Allele.create("A", true),          Allele.create("T"),          contig, 999,  999, Strand.NEGATIVE, exonPositionList, interval5 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 891,  900, Strand.NEGATIVE, exonPositionList, interval5 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 920,  929, Strand.NEGATIVE, exonPositionList, interval5 },
                { Allele.create("ATGCATGCAT", true), Allele.create("GCATGCATGC"), contig, 999, 1008, Strand.NEGATIVE, exonPositionList, interval5 },

        };
    }

    @DataProvider
    Object[][] provideDataForTestAssertValidStrand_ValidStrands() {
        return new Object[][] {
                { Strand.POSITIVE },
                { Strand.NEGATIVE },
        };
    }

    @DataProvider
    Object[][] provideDataForTestAssertValidStrand_InvalidStrands() {
        return new Object[][] {
                { null },
                { Strand.NONE }
        };
    }

    @DataProvider
    Object[][] provideDataForTestIsPositionInFrame() {
        return new  Object[][] {
                {   1, true },
                {   2, false },
                {   3, false },
                {   4, true },
                {   5, false },
                {   6, false },
                {   7, true },
                {   8, false },
                {   9, false },
                {  10, true },
                { 289, true },
                { 290, false },
                { 291, false },
                { 292, true },
                { 293, false },
                { 294, false },

        };
    }

    @DataProvider
    Object[][] provideDataForTestGetBasesInWindowAroundReferenceAllele() {

        final int offset = 100000000;

        return new Object[][] {
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  1, 500000 + offset, 500000 + offset,          "TAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  2, 500000 + offset, 500000 + offset,         "GTACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  3, 500000 + offset, 500000 + offset,        "AGTACAT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  4, 500000 + offset, 500000 + offset,       "TAGTACATT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  5, 500000 + offset, 500000 + offset,      "TTAGTACATTA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  6, 500000 + offset, 500000 + offset,     "GTTAGTACATTAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  7, 500000 + offset, 500000 + offset,    "AGTTAGTACATTAAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  8, 500000 + offset, 500000 + offset,   "TAGTTAGTACATTAACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+",  9, 500000 + offset, 500000 + offset,  "CTAGTTAGTACATTAACAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "+", 10, 500000 + offset, 500000 + offset, "ACTAGTTAGTACATTAACAAT"),

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  1, 500000 + offset, 500003 + offset,          "TACATT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  2, 500000 + offset, 500003 + offset,         "GTACATTA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  3, 500000 + offset, 500003 + offset,        "AGTACATTAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  4, 500000 + offset, 500003 + offset,       "TAGTACATTAAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  5, 500000 + offset, 500003 + offset,      "TTAGTACATTAACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  6, 500000 + offset, 500003 + offset,     "GTTAGTACATTAACAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  7, 500000 + offset, 500003 + offset,    "AGTTAGTACATTAACAAT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  8, 500000 + offset, 500003 + offset,   "TAGTTAGTACATTAACAATG"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+",  9, 500000 + offset, 500003 + offset,  "CTAGTTAGTACATTAACAATGA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT",  "+", 10, 500000 + offset, 500003 + offset, "ACTAGTTAGTACATTAACAATGAC"),

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  1, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(         "TAC".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  2, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(        "GTACA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  3, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(       "AGTACAT".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  4, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(      "TAGTACATT".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  5, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(     "TTAGTACATTA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  6, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(    "GTTAGTACATTAA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  7, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(   "AGTTAGTACATTAAC".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  8, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(  "TAGTTAGTACATTAACA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-",  9, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement( "CTAGTTAGTACATTAACAA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "-", 10, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement("ACTAGTTAGTACATTAACAAT".getBytes() )) ) ,

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  1, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(          "TACATT".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  2, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(         "GTACATTA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  3, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(        "AGTACATTAA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  4, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(       "TAGTACATTAAC".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  5, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(      "TTAGTACATTAACA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  6, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(     "GTTAGTACATTAACAA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  7, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(    "AGTTAGTACATTAACAAT".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  8, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(   "TAGTTAGTACATTAACAATG".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-",  9, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement(  "CTAGTTAGTACATTAACAATGA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "-", 10, 500000 + offset, 500003 + offset, new String( BaseUtils.simpleReverseComplement( "ACTAGTTAGTACATTAACAATGAC".getBytes() )) ),
        };
    }

    @DataProvider
    private Object[][] provideDataForTestCreateReferenceSnippet() {

        // NOTE: Genome positions start at 1, not 0.
        //                                                                                                                          1
        //                       0        1         2         3         4    *    5         6         7         8         9         0
        //                       1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        final String sequence = "CGTCGACGGAACAAAAGTAGACCATCCCTCTTGGTAAGTACGTCTTCATACTCTACAAATACCCATAGCACAATTCGGAGCCCAACGCCCGACGGGTCAT";
        //   REVERSE COMPLEMENT: ATGACCCGTCGGGCGTTGGGCTCCGAATTGTGCTATGGGTATTTGTAGAGTATGAAGACGTACTTACCAAGAGGGATGGTCTACTTTTGTTCCGTCGACG
        //                       1                                                      *
        //                       0         9         8         7         6         5         4         3         2         1
        //                       0987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321

        final String contig = "specialTestContig";

        final int startPos = 45;
        final int endPos = 45;

        // NOTE: Variants are always described in the + direction.
        //       The start/end positions are for the REFERENCE allele ONLY.

        return new Object[][] {
                // SNP in + direction
                {
                        Allele.create("T", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTA"
                },
                // DNP in + direction
                {
                        Allele.create("TT", true),
                        Allele.create("CC"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos+1),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTAC"
                },
                // TNP in + direction
                {
                        Allele.create("TTC", true),
                        Allele.create("CGG"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos+2),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTACA"
                },
                // Insertion - 1 base in + direction
                {
                        Allele.create("T", true),
                        Allele.create("TC"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "AAGTACGTCTTCATACTCTA"
                },
                // Insertion - 2 bases in + direction
                {
                        Allele.create("T", true),
                        Allele.create("TCC"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "AAGTACGTCTTCATACTCTA"
                },
                // Insertion - 3 bases in + direction
                {
                        Allele.create("T", true),
                        Allele.create("TCCA"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "AAGTACGTCTTCATACTCTA"
                },
                // Deletion - 1 base in + direction
                {
                        Allele.create("CT", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-1,  endPos),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTA"
                },
                // Deletion - 2 bases in + direction
                {
                        Allele.create("CTT", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-1,  endPos+1),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTAC"
                },
                // Deletion - 3 bases in + direction
                {
                        Allele.create("CTTC", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig, startPos - 1, endPos + 2),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTACA"
                },

                // ================================================================================

                // SNP in - direction
                {
                        Allele.create("T", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTA"
                },
                // DNP in - direction
                {
                        Allele.create("CT", true),
                        Allele.create("TA"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-1,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTAC"
                },
                // TNP in - direction
                {
                        Allele.create("TCT", true),
                        Allele.create("ATG"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-2,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTACC"
                },
                // Insertion - 1 base in - direction
                {
                        Allele.create("T", true),
                        Allele.create("TC"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTT"
                },
                // Insertion - 2 bases in - direction
                {
                        Allele.create("T", true),
                        Allele.create("TCC"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTT"
                },
                // Insertion - 3 bases in - direction
                {
                        Allele.create("T", true),
                        Allele.create("TCCA"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTT"
                },


        //                       0        1         2         3         4    *    5         6         7         8         9         0
        //                       1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
      //final String sequence = "CGTCGACGGAACAAAAGTAGACCATCCCTCTTGGTAAGTACGTCTTCATACTCTACAAATACCCATAGCACAATTCGGAGCCCAACGCCCGACGGGTCAT";
        //   REVERSE COMPLEMENT: ATGACCCGTCGGGCGTTGGGCTCCGAATTGTGCTATGGGTATTTGTAGAGTATGAAGACGTACTTACCAAGAGGGATGGTCTACTTTTGTTCCGTCGACG
        //                       1                                                      *
        //                       0         9         8         7         6         5         4         3         2         1
        //                       0987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321

                // Deletion - 1 base in - direction
                {
                        Allele.create("CT", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-1,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTA"
                },
                // Deletion - 2 bases in - direction
                {
                        Allele.create("TCT", true),
                        Allele.create("T"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-2,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTAC"
                },
                // Deletion - 3 bases in - direction
                {
                        Allele.create("GTCT", true),
                        Allele.create("G"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-3,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTACC"
                },
                // Deletion - 4 bases in - direction
                {
                        Allele.create("CGTCT", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(sequence, contig,  startPos-4,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTACCA"
                },

                // ================================================================================

                // Special test case 1:
                {
                        Allele.create("CT", true),
                        Allele.create("C"),
                        FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(
                                //         0        1     *   2         3
                                //         123456789012345678901234567890
                                "ACCCAGCCCACTCACCTTTCTCTCCTGGAA",
                                contig,
                                16,
                                17 ),
                        Strand.NEGATIVE,
                        "CAGGAGAGAAAGGTGAGTGGG"
                }
        };
    }

    @DataProvider
    Object[][] provideForTestGetStrandCorrectedAllele() {
        return new Object[][] {
                { Allele.create("A"),       Strand.POSITIVE, Allele.create("A") },
                { Allele.create("AA"),      Strand.POSITIVE, Allele.create("AA") },
                { Allele.create("AAT"),     Strand.POSITIVE, Allele.create("AAT") },
                { Allele.create("AATT"),    Strand.POSITIVE, Allele.create("AATT") },
                { Allele.create("AATTG"),   Strand.POSITIVE, Allele.create("AATTG") },
                { Allele.create("AATTGC"),  Strand.POSITIVE, Allele.create("AATTGC") },
                { Allele.create("AATTGCG"), Strand.POSITIVE, Allele.create("AATTGCG") },
                { Allele.create("A"),       Strand.NEGATIVE, Allele.create("T") },
                { Allele.create("AA"),      Strand.NEGATIVE, Allele.create("TT") },
                { Allele.create("AAT"),     Strand.NEGATIVE, Allele.create("ATT") },
                { Allele.create("AATT"),    Strand.NEGATIVE, Allele.create("AATT") },
                { Allele.create("AATTG"),   Strand.NEGATIVE, Allele.create("CAATT") },
                { Allele.create("AATTGC"),  Strand.NEGATIVE, Allele.create("GCAATT") },
                { Allele.create("AATTGCG"), Strand.NEGATIVE, Allele.create("CGCAATT") },
        };
    }

    //==================================================================================================================
    // Tests:

//    @Test(dataProvider = "provideReferenceAndExonListAndExpected")
//    void testGetCodingSequence(final ReferenceContext reference, final List<Locatable> exonList, final Strand strand, final String expected) {
//        final String codingSequence = FuncotatorUtils.getCodingSequence(reference, exonList, strand);
//        Assert.assertEquals( codingSequence, expected );
//    }

//    @Test(dataProvider = "provideReferenceAndExonListForGatkExceptions",
//            expectedExceptions = GATKException.class)
//    void testGetCodingSequenceWithGatkExceptions(final ReferenceContext reference, final Strand strand, final List<Locatable> exonList) {
//        FuncotatorUtils.getCodingSequence(reference, exonList, strand);
//    }

//    @Test(dataProvider = "provideReferenceAndExonListForIllegalArgumentExceptions",
//            expectedExceptions = IllegalArgumentException.class)
//    void testGetCodingSequenceWithIllegalArgumentExceptions(final ReferenceContext reference, final List<Locatable> exonList) {
//        FuncotatorUtils.getCodingSequence(reference, exonList);
//    }

    @Test(dataProvider = "provideForTestGetIndelAdjustedAlleleChangeStartPosition")
    void testGetIndelAdjustedAlleleChangeStartPosition(final VariantContext variant, final int expected) {
         Assert.assertEquals( FuncotatorUtils.getIndelAdjustedAlleleChangeStartPosition(variant), expected );
    }

    @Test(dataProvider = "provideDataForTestGetNonOverlappingAltAlleleBaseString")
    void testGetNonOverlappingAltAlleleBaseString(final Allele refAllele, final Allele altAllele, final boolean copyRefBasesWhenAltIsPastEnd, final String expected) {
        Assert.assertEquals( FuncotatorUtils.getNonOverlappingAltAlleleBaseString(refAllele, altAllele, copyRefBasesWhenAltIsPastEnd), expected);
    }

    @Test(dataProvider = "provideDataForGetStartPositionInTranscript")
    void testGetStartPositionInTranscript(final Locatable variant, final List<? extends Locatable> transcript, final Strand strand, final int expected) {
        Assert.assertEquals( FuncotatorUtils.getStartPositionInTranscript(variant, transcript, strand), expected );
    }

    @Test(dataProvider = "providePositionAndExpectedAlignedPosition")
    void testGetAlignedPosition(final int pos, final int expected) {
        Assert.assertEquals(FuncotatorUtils.getAlignedPosition(pos), expected);
    }

    @Test(dataProvider = "provideDataForTestGetAlignedEndPositionOneArg")
    void testGetAlignedEndPositionOneArg(final int alleleEndPos, final int expected) {
        Assert.assertEquals(FuncotatorUtils.getAlignedEndPosition(alleleEndPos), expected);
    }

    @Test(dataProvider = "provideDataForGetAlternateSequence")
    void testGetAlternateSequence(final String refCodingSeq, final int startPos, final Allele refAllele, final Allele altAllele, final String expected) {
        Assert.assertEquals(FuncotatorUtils.getAlternateSequence(new StrandCorrectedReferenceBases(refCodingSeq,Strand.POSITIVE), startPos, refAllele, altAllele, Strand.POSITIVE), expected);
    }

    @Test(dataProvider = "provideDataForGetEukaryoticAminoAcidByCodon")
    void testGetEukaryoticAminoAcidByCodon(final String codon, final AminoAcid expected) {
        Assert.assertEquals(FuncotatorUtils.getEukaryoticAminoAcidByCodon(codon), expected);
    }

    @Test(dataProvider = "provideDataForGetMitochondrialAminoAcidByCodon")
    void testGetMitochondrialAminoAcidByCodon(final String codon,
                                              final boolean isFirst,
                                              final FuncotatorUtils.Genus genus,
                                              final AminoAcid expected) {
        Assert.assertEquals(FuncotatorUtils.getMitochondrialAminoAcidByCodon(codon, isFirst, genus), expected);

        // Make sure to test the other case as well:
        if ( !isFirst || (genus == FuncotatorUtils.Genus.UNSPECIFIED) ) {
            Assert.assertEquals(FuncotatorUtils.getMitochondrialAminoAcidByCodon(codon), expected);
        }
    }

    @Test(dataProvider = "provideDataForGetCodonChangeString")
    void testGetCodonChangeString( final String refAllele,
                                   final String altAllele,
                                   final int alleleStart,
                                   final int codingSequenceAlleleStart,
                                   final int alignedCodingSequenceAlleleStart,
                                   final String alignedCodingSeqRefAllele,
                                   final String alignedCodingSeqAltAllele,
                                   final String alignedAlternateAllele,
                                   final int alignedRefAlleleStop,
                                   final String contig,
                                   final Strand strand,
                                   final ReferenceSequence codingSequence,
                                   final Locatable startCodon,
                                   final String expected ) {

        final SequenceComparison seqComp = new SequenceComparison();
        seqComp.setReferenceAllele(refAllele);
        seqComp.setAlternateAllele(altAllele);
        seqComp.setAlleleStart(alleleStart);
        seqComp.setCodingSequenceAlleleStart(codingSequenceAlleleStart);
        seqComp.setAlignedCodingSequenceAlleleStart(alignedCodingSequenceAlleleStart);
        seqComp.setAlignedCodingSequenceReferenceAllele(alignedCodingSeqRefAllele);
        seqComp.setAlignedCodingSequenceAlternateAllele(alignedCodingSeqAltAllele);
        seqComp.setAlignedAlternateAllele(alignedAlternateAllele);
        seqComp.setAlignedReferenceAlleleStop(alignedRefAlleleStop);
        seqComp.setContig(contig);
        seqComp.setStrand(strand);
        seqComp.setTranscriptCodingSequence(codingSequence);

        Assert.assertEquals(FuncotatorUtils.getCodonChangeString(seqComp, startCodon), expected);
    }

    @Test
    void testGetAminoAcidNames() {
        Assert.assertEquals(FuncotatorUtils.getAminoAcidNames(),
                new String[]{
                        "Alanine",
                        "Arganine",
                        "Asparagine",
                        "Aspartic acid",
                        "Cysteine",
                        "Glutamic acid",
                        "Glutamine",
                        "Glycine",
                        "Histidine",
                        "Isoleucine",
                        "Leucine",
                        "Lysine",
                        "Methionine",
                        "Phenylalanine",
                        "Proline",
                        "Serine",
                        "Stop codon",
                        "Threonine",
                        "Tryptophan",
                        "Tyrosine",
                        "Valine",
                        "Undecodable Amino Acid"
                }
        );
    }

    @Test
    void testGetAminoAcidCodes() {
        Assert.assertEquals(FuncotatorUtils.getAminoAcidCodes(),
                new String[] {
                        "Ala",
                        "Arg",
                        "Asn",
                        "Asp",
                        "Cys",
                        "Glu",
                        "Gln",
                        "Gly",
                        "His",
                        "Ile",
                        "Leu",
                        "Lys",
                        "Met",
                        "Phe",
                        "Pro",
                        "Ser",
                        "Stop",
                        "Thr",
                        "Trp",
                        "Tyr",
                        "Val",
                        "UNDECODABLE"
                }
        );
    }

    @Test (dataProvider = "provideDataForIsInFrameWithEndOfRegion")
    void testIsInFrameWithEndOfRegion(final int pos, final int length, final boolean expected) {
        Assert.assertEquals( FuncotatorUtils.isInFrameWithEndOfRegion(pos, length), expected );
    }

    @Test (dataProvider = "provideStringDataForGetAminoAcidByLetter")
    void testGetAminoAcidByLetter(final String letter, final AminoAcid expected) {
        Assert.assertEquals( FuncotatorUtils.getAminoAcidByLetter(letter), expected );
    }

    @Test (dataProvider = "provideCharDataForGetAminoAcidByLetter")
    void testGetAminoAcidByLetter(final char letter, final AminoAcid expected) {
        Assert.assertEquals( FuncotatorUtils.getAminoAcidByLetter(letter), expected );
    }

    @Test (dataProvider = "provideDataForGetProteinChangePosition")
    void testGetProteinChangePosition(final Integer alignedCodingSequenceStartPos, final int expected) {
        Assert.assertEquals( FuncotatorUtils.getProteinChangePosition(alignedCodingSequenceStartPos) , expected );
    }

    @Test (dataProvider = "provideDataForRenderProteinChangeString")
    void testRenderProteinChangeString(final int protChangeStartPos,
                                    final int protChangeEndPos,
                                    final String refAminoAcidSeq,
                                    final String altAminoAcidSeq,
                                    final String refAllele,
                                    final String altAllele,
                                    final Locatable startCodon,
                                    final int refAlleleStartPosition,
                                    final String expected) {

        final SequenceComparison seqComp = new SequenceComparison();

        seqComp.setProteinChangeInfo(
                ProteinChangeInfo.create(
                        protChangeStartPos,
                        protChangeEndPos,
                        refAminoAcidSeq,
                        altAminoAcidSeq
                )
        );

        seqComp.setAlleleStart(refAlleleStartPosition);

        seqComp.setReferenceAllele(refAllele);
        seqComp.setAlternateAllele(altAllele);
        if (startCodon != null) {
            seqComp.setContig(startCodon.getContig());
        }
        else {
            seqComp.setContig("TESTCONTIG");
        }

        Assert.assertEquals( FuncotatorUtils.renderProteinChangeString(seqComp, startCodon), expected );
    }

    @Test(dataProvider = "provideForGetClosestExonIndex")
    void testGetClosestExonIndex( final int variantStartPos,
                             final List<? extends Locatable> exonList,
                             final int expected ) {
        Assert.assertEquals( FuncotatorUtils.getClosestExonIndex(variantStartPos, exonList), expected );
    }

    @Test(dataProvider = "provideForCreateIntronicCDnaString")
    void testCreateIntronicCDnaString(final int variantStartPos,
                                      final List<? extends Locatable> exonList,
                                      final String strandCorrectedRefAllele,
                                      final String strandCorrectedAltAllele,
                                      final String expected ) {
        Assert.assertEquals( FuncotatorUtils.createIntronicCDnaString(variantStartPos, exonList, strandCorrectedRefAllele, strandCorrectedAltAllele), expected );
    }

    @Test(dataProvider = "provideForTestGetCodingSequenceChangeString")
    void testGetCodingSequenceChangeString(final int codingSequenceAlleleStart,
                                           final String referenceAllele,
                                           final String alternateAllele,
                                           final Strand strand,
                                           final String expected) {

        Assert.assertEquals( FuncotatorUtils.getCodingSequenceChangeString(codingSequenceAlleleStart, referenceAllele, alternateAllele, strand, null, null, null), expected );
    }

    @Test(dataProvider = "provideForTestIsIndelBetweenCodons")
    void testIsIndelBetweenCodons(final int codingSequenceAlleleStart,
                                  final int alignedCodingSequenceAlleleStart,
                                  final String refAllele,
                                  final Strand strand,
                                  final boolean expected) {

        Assert.assertEquals(
                FuncotatorUtils.isIndelBetweenCodons(codingSequenceAlleleStart, alignedCodingSequenceAlleleStart, refAllele, strand),
                expected
        );
    }

    @Test (dataProvider = "provideDataForGetProteinChangeEndPosition")
    void testGetProteinChangeEndPosition(final Integer proteinChangeStartPosition, final Integer alignedAlternateAlleleLength, final int expected) {
        Assert.assertEquals( FuncotatorUtils.getProteinChangeEndPosition(proteinChangeStartPosition, alignedAlternateAlleleLength) , expected );
    }

    @Test (dataProvider = "provideDataForGetAlignedCodingSequenceAllele")
    void testGetAlignedCodingSequenceAllele(  final String codingSequence,
                                final Integer alignedAlleleStart,
                                final Integer alignedAlleleStop,
                                final Allele refAllele,
                                final Integer refAlleleStart,
                                final Strand strand,
                                final String expected) {
        final String alignedRefAllele = FuncotatorUtils.getAlignedCodingSequenceAllele(codingSequence, alignedAlleleStart, alignedAlleleStop, refAllele, refAlleleStart, strand);
        Assert.assertEquals( alignedRefAllele, expected );
    }

    @Test(dataProvider = "provideDataForTestGetAlignedRefAllele")
    void testGetAlignedRefAllele( final String referenceSnippet,
                                  final int referencePadding,
                                  final Allele refAllele,
                                  final int codingSequenceRefAlleleStart,
                                  final int alignedRefAlleleStart,
                                  final Strand strand,
                                  final String expected) {

        // Make a Dummy Locatable for Logging:
        final Locatable dummyLocatableForLogging = new SimpleInterval("ReferenceSnippet", 1, 100);

        final Allele altAllele = Allele.create(Utils.dupChar('A', refAllele.length()));

        Assert.assertEquals(
                FuncotatorUtils.getAlignedRefAllele(
                        new StrandCorrectedReferenceBases(referenceSnippet, strand),
                        referencePadding,
                        refAllele,
                        altAllele,
                        codingSequenceRefAlleleStart,
                        alignedRefAlleleStart,
                        strand,
                        dummyLocatableForLogging),
                expected
        );
    }

    @Test (dataProvider = "provideDataForGetTranscriptAlleleStartPosition")
    void testGetTranscriptAlleleStartPosition(final VariantContext variant, final List<Locatable> exons, final Strand strand, final int expected) {
        Assert.assertEquals( FuncotatorUtils.getTranscriptAlleleStartPosition(variant, exons, strand), expected );
    }

    @Test (dataProvider = "provideDataForTestCreateSpliceSiteCodonChange")
    void testCreateSpliceSiteCodonChange(final int variantStart,
                                         final int exonNumber,
                                         final int exonStart,
                                         final int exonEnd,
                                         final Strand strand,
                                         final int offsetIndelAdjustment,
                                         final String expected) {

        Assert.assertEquals( FuncotatorUtils.createSpliceSiteCodonChange(variantStart, exonNumber, exonStart, exonEnd, strand, offsetIndelAdjustment), expected );
    }

    @Test (dataProvider = "provideDataForTestGetOverlappingExonPositions")
    void testGetOverlappingExonPositions(final Allele refAllele,
                                         final Allele altAllele,
                                         final String contig,
                                         final int start,
                                         final int stop,
                                         final Strand strand,
                                         final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                         final SimpleInterval expected) {
        Assert.assertEquals( FuncotatorUtils.getOverlappingExonPositions(refAllele, altAllele, contig, start, stop, strand, exonPositionList), expected);
    }

    @Test (dataProvider = "provideDataForTestAssertValidStrand_ValidStrands")
    void testAssertValidStrand_ValidStrands( final Strand strand ) {
        FuncotatorUtils.assertValidStrand( strand );
    }

    @Test (dataProvider = "provideDataForTestAssertValidStrand_InvalidStrands",
           expectedExceptions = {GATKException.class, IllegalArgumentException.class})
    void testAssertValidStrand_InvalidStrands( final Strand strand ) {
        FuncotatorUtils.assertValidStrand( strand );
    }

    @Test (dataProvider = "provideDataForTestIsPositionInFrame")
    void testIsPositionInFrame( final int position, final boolean expected ) {
        Assert.assertEquals( FuncotatorUtils.isPositionInFrame(position), expected );
    }

    @Test (dataProvider = "provideDataForTestGetBasesInWindowAroundReferenceAllele")
    void testGetBasesInWindowAroundReferenceAllele(final Allele refAllele,
                                                   final ReferenceContext referenceContext,
                                                   final Strand strand,
                                                   final int referenceWindow,
                                                   final String expected) {

        final StrandCorrectedReferenceBases basesInWindow = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(refAllele, referenceContext, strand, referenceWindow);
        Assert.assertEquals( basesInWindow, new StrandCorrectedReferenceBases(expected, strand) );
    }

    @Test(dataProvider = "provideDataForTestCreateReferenceSnippet")
    void testCreateReferenceSnippet(final Allele refAllele, final Allele altAllele, final ReferenceContext reference, final Strand strand, final String expected ) {
        final int referenceWindow = 10;
        Assert.assertEquals( FuncotatorUtils.createReferenceSnippet(refAllele, altAllele, reference, strand, referenceWindow), new StrandCorrectedReferenceBases(expected, strand));
    }

    @Test(enabled = false)
    public void testSequenceDictionaryMD5Sums() {

        final String refDir  = "/Users/jonn/Development/references" + File.separator;

        final Map<String, List<String>> nameFilenameMap = new LinkedHashMap<>();

        nameFilenameMap.put("BROAD REF", Arrays.asList("Homo_sapiens_assembly19.fasta", ""));
        nameFilenameMap.put("B37",       Arrays.asList("human_g1k_v37.fasta", ""));
        nameFilenameMap.put("GRCh37",    Arrays.asList("GRCh37.p13.genome.fasta", "chr"));
        nameFilenameMap.put("HG19",      Arrays.asList("ucsc.hg19.fasta", "chr"));

        final int chr3Length = 198022430;
        final int chrYLength = 59373566;

        for ( final Map.Entry<String, List<String>> entry : nameFilenameMap.entrySet() ) {

            final String fileName = entry.getValue().get(0);
            final String contigDecorator = entry.getValue().get(1);

            final ReferenceDataSource referenceDataSource = ReferenceDataSource.of(IOUtils.getPath(refDir + fileName), true);

            System.out.println();
            System.out.println("================================================================================");
            System.out.println("|                             " + entry.getKey() + " (" + fileName + ")");
            System.out.println("================================================================================");

            Assert.assertEquals(referenceDataSource.queryAndPrefetch(contigDecorator + "3", 1, chr3Length).getBases().length, chr3Length);
            Assert.assertEquals(referenceDataSource.queryAndPrefetch(contigDecorator + "Y", 1, chrYLength).getBases().length, chrYLength);

            System.out.println("CHR 3: " + Utils.calcMD5(referenceDataSource.queryAndPrefetch(contigDecorator + "3", 1, chr3Length).getBases()));
            System.out.println("CHR Y: " + Utils.calcMD5(referenceDataSource.queryAndPrefetch(contigDecorator + "Y", 1, chrYLength).getBases()));
        }
    }

    @Test(dataProvider = "provideForTestGetStrandCorrectedAllele")
    public void testGetStrandCorrectedAllele(final Allele allele, final Strand strand, final Allele expected) {
         Assert.assertEquals( FuncotatorUtils.getStrandCorrectedAllele(allele, strand), expected );
    }

    @DataProvider
    public Object[][] provideIsGencodeFuncotation() {
        return new Object[][] {
                {TableFuncotation.create(Collections.singletonList("FIELD"), Collections.singletonList("VALUE"), Allele.create("A"), "TEST", null), false},
                {new GencodeFuncotationBuilder().setAnnotationTranscript("TXID").build(), true}
        };
    }

    @Test(dataProvider = "provideIsGencodeFuncotation")
    public void testIsGencodeFuncotation(final Funcotation f, final boolean gt) {
        Assert.assertEquals(FuncotatorUtils.isGencodeFuncotation(f), gt);
    }

    @DataProvider
    public Object[][] provideAreGencodeFuncotations() {
        return new Object[][] {
                {Collections.singletonList(TableFuncotation.create(Collections.singletonList("FIELD"), Collections.singletonList("VALUE"), Allele.create("A"), "TEST", null)), false},
                {Collections.singletonList(new GencodeFuncotationBuilder().setAnnotationTranscript("TXID").build()), true},
                {Arrays.asList(
                        TableFuncotation.create(Collections.singletonList("FIELD"), Collections.singletonList("VALUE"), Allele.create("A"), "TEST", null),
                        TableFuncotation.create(Collections.singletonList("FIELD2"), Collections.singletonList("VALUE1"), Allele.create("A"), "TEST", null)
                    ), false
                },
                {Arrays.asList(
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID1").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID2").build()
                    ), true
                },
                {Arrays.asList(
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID1").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID2").build(),
                        TableFuncotation.create(Collections.singletonList("FIELD"), Collections.singletonList("VALUE"), Allele.create("A"), "TEST", null)
                    ), true
                },
                {Arrays.asList(
                        TableFuncotation.create(Collections.singletonList("FIELD"), Collections.singletonList("VALUE"), Allele.create("A"), "TEST", null),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID1").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID2").build()
                    ), true
                }
        };
    }

    @Test(dataProvider = "provideAreGencodeFuncotations")
    public void testAreGencodeFuncotations(final List<Funcotation> f, final boolean gt) {
        Assert.assertEquals(FuncotatorUtils.areAnyGencodeFuncotation(f), gt);
    }

    @DataProvider
    public Object[][] provideCreateFuncotations() {
        final Map<String, String> attributes1 = ImmutableMap.of("FOOFIELD", "FOO", "BAZFIELD", "BAZ");
        final List<VCFInfoHeaderLine> attributes1AsVcfHeaderLine = attributes1.keySet().stream()
                .map(k -> new VCFInfoHeaderLine(k, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Description here"))
                .collect(Collectors.toList());

        final Map<String, String> attributes2 = ImmutableMap.of("FOOFIELD", "FOO,FOOINDEL", "BAZFIELD", "BAZ,BAZINDEL");
        final List<VCFInfoHeaderLine> attributes2AsVcfHeaderLine = attributes2.keySet().stream()
                .map(k -> new VCFInfoHeaderLine(k, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Description here"))
                .collect(Collectors.toList());

        final List<VCFInfoHeaderLine> attributes2AsVcfHeaderLineWithExtras = attributes2.keySet().stream()
                .map(k -> new VCFInfoHeaderLine(k, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Description here"))
                .collect(Collectors.toList());
        attributes2AsVcfHeaderLineWithExtras.add(new VCFInfoHeaderLine("EXTRA", VCFHeaderLineCount.A, VCFHeaderLineType.String, "Description here"));


        return new Object[][] {
            { new VariantContextBuilder(
                    FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                    "chr3",
                    1000000,
                    1000000,
                    Arrays.asList(Allele.create("A", true), Allele.create("C")))
                    .attributes(attributes1)
                    .make(),
                    VcfFuncotationMetadata.create(attributes1AsVcfHeaderLine),
                    "TEST1"
            }, { new VariantContextBuilder(
                    FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                    "chr3",
                    1000000,
                    1000000,
                    Arrays.asList(Allele.create("A", true), Allele.create("C"), Allele.create("ATT")))
                    .attributes(attributes1)
                    .make(),
                    VcfFuncotationMetadata.create(attributes2AsVcfHeaderLine),
                    "TEST1"
            }, { new VariantContextBuilder(
                    FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                    "chr3",
                    1000000,
                    1000000,
                    Arrays.asList(Allele.create("A", true), Allele.create("C"), Allele.create("ATT")))
                    .attributes(attributes1)
                    .make(),
                    VcfFuncotationMetadata.create(attributes2AsVcfHeaderLineWithExtras),
                    "TEST1"
            }
        };
    }

    @Test(dataProvider = "provideCreateFuncotations")
    public void testCreateFuncotations(final VariantContext vc, final FuncotationMetadata metadata, final String datasourceName) {
        final List<Funcotation> funcotations = FuncotatorUtils.createFuncotations(vc, metadata, datasourceName);

        Assert.assertTrue(funcotations.stream().allMatch(f -> f.getDataSourceName().equals(datasourceName)));
        Assert.assertEquals(funcotations.stream().map(Funcotation::getAltAllele).collect(Collectors.toSet()), new HashSet<>(vc.getAlternateAlleles()));
        Assert.assertEquals(funcotations.stream().map(Funcotation::getMetadata).collect(Collectors.toSet()), new HashSet<>(Collections.singletonList(metadata)));
    }

    @DataProvider
    public Object[][] provideMafSanitizing() {
        return new Object[][] {
                {"\tHAHAHAH\t", "_%09_HAHAHAH_%09_"},
                {"\tHAHAHAH\n", "_%09_HAHAHAH_%0A_"},
                {"FOO", "FOO"}
        };
    }

    @Test(dataProvider = "provideMafSanitizing")
    public void testSanitizeFuncotationFieldForMaf(final String individualFuncotationField, final String gt) {
        final String guess = FuncotatorUtils.sanitizeFuncotationFieldForMaf(individualFuncotationField);
        Assert.assertEquals(guess, gt);

    }

    @SuppressWarnings("unchecked")
    @DataProvider
    public Object[][] provideForRenderSanitizedFuncotationForVcf() {

        return new Object[][]{
                // Test a very basic case where all fields are included
                {OutputRenderer.createFuncotationFromLinkedHashMap(
                        (LinkedHashMap) MapUtils.putAll(new LinkedHashMap<String, String>(),
                        new String[][]{{"FOO", "BAR"},{"BAZ", "HUH?"}}),
                        Allele.create("T"), "FAKEDATA"), Arrays.asList("FOO", "BAZ"),
                        "BAR|HUH?"
                },

                // Test case where only one field is included
                {OutputRenderer.createFuncotationFromLinkedHashMap(
                        (LinkedHashMap) MapUtils.putAll(new LinkedHashMap<String, String>(),
                        new String[][]{{"FOO", "BAR"},{"BAZ", "HUH?"}}),
                        Allele.create("T"), "FAKEDATA"), Arrays.asList("FOO"),
                        "BAR"
                },

                // Make sure that specifying a non-existent included field (NOTHERE) has no effect on the output,
                //  even when another field is excluded.
                {OutputRenderer.createFuncotationFromLinkedHashMap(
                        (LinkedHashMap) MapUtils.putAll(new LinkedHashMap<String, String>(),
                        new String[][]{{"FOO", "BAR"},{"BAZ", "HUH?"}}),
                        Allele.create("T"), "FAKEDATA"), Arrays.asList("FOO", "NOTHERE"),
                        "BAR"
                },

                // Make sure that specifying a non-existent included field (NOTHERE) has no effect on the output,
                //  even when all fields are included..
                {OutputRenderer.createFuncotationFromLinkedHashMap(
                        (LinkedHashMap) MapUtils.putAll(new LinkedHashMap<String, String>(),
                        new String[][]{{"FOO", "BAR"},{"BAZ", "HUH?"}}),
                        Allele.create("T"), "FAKEDATA"), Arrays.asList("FOO", "BAZ", "NOTHERE"),
                        "BAR|HUH?"
                }
        };
    }

    @Test(dataProvider = "provideForRenderSanitizedFuncotationForVcf" )
    public void testRenderSanitizedFuncotationForVcf(final Funcotation funcotation, final List<String> includedFields, final String gt) {
        final String guess = FuncotatorUtils.renderSanitizedFuncotationForVcf(funcotation, includedFields);
        Assert.assertEquals(guess, gt);
    }

    @DataProvider
    public Object[][] provideForTestSanitizeFuncotationFieldForVcf() {
        return new Object[][] {
                { "", "" },
                { ",", "_%2C_" },
                { ";", "_%3B_" },
                { "=", "_%3D_" },
                { "\t", "_%09_" },
                { "|", "_%7C_" },
                { " ", "_%20_" },
                { "\n", "_%0A_" },
                { "#", "_%23_" },
                { ",;=\t| \n#", "_%2C__%3B__%3D__%09__%7C__%20__%0A__%23_" },
                { "ASDFJKL:", "ASDFJKL:" },
                { "ASDF;JKL:", "ASDF_%3B_JKL:" },
        };
    }

    @Test(dataProvider = "provideForTestSanitizeFuncotationFieldForVcf" )
    public void testSanitizeFuncotationFieldForVcf(final String input, final String expected) {
        Assert.assertEquals( FuncotatorUtils.sanitizeFuncotationFieldForVcf(input), expected );
    }

    @DataProvider
    public Object[][] provideCreateLinkedHashMapFromLists() {
        final LinkedHashMap<String,String> gtMap1 = new LinkedHashMap<>();
        gtMap1.put("K1", "V1");
        gtMap1.put("K2", "V2");
        final LinkedHashMap<Allele,Allele> gtMap2 = new LinkedHashMap<>();
        gtMap2.put(Allele.create("C"), Allele.create("G"));
        gtMap2.put(Allele.create("T"), Allele.create("A"));

        return new Object[][] {
                {Arrays.asList("K1", "K2"), Arrays.asList("V1", "V2"), gtMap1},
                {Arrays.asList(Allele.create("C"), Allele.create("T")), Arrays.asList(Allele.create("G"), Allele.create("A")), gtMap2}
        };
    }


    @Test(dataProvider = "provideCreateLinkedHashMapFromLists")
    public <T,U> void testCreateLinkedHashMapFromLists(final List<T> keys, final List<U> values, final LinkedHashMap<T, U> gtMap) {
        final LinkedHashMap<T, U> guess = FuncotatorUtils.createLinkedHashMapFromLists(keys, values);
        Assert.assertEquals(guess, gtMap);
    }

    @DataProvider
    public Object[][] provideCreateLinkedHashMapFromListsWithIllegalArgs() {
        return new Object[][] {
                // Lengths don't match
                {Arrays.asList("K1", "K2"), Collections.singletonList("V1")},
                {Arrays.asList("K1", "K2"), Arrays.asList("V1", "V2", "V3")},
                {Arrays.asList("K1", "K2"), Collections.emptyList()},
                {Collections.singletonList("K1"), Arrays.asList("V1", "V2")},
                {Arrays.asList("K1", "K2", "K3"), Arrays.asList("V1", "V2")},
                {Collections.emptyList(), Arrays.asList("V1", "V2")},
                // Duplicate keys
                {Arrays.asList("K1", "K2", "K1"), Arrays.asList("V1", "V2", "V3")},
                // Duplicate keys, even if the value is the same.  LinkedHashMap implies an ordering, which falls
                //  apart if keys are not unique.
                {Arrays.asList("K1", "K2", "K1"), Arrays.asList("V1", "V2", "V1")}
        };
    }
    @Test(dataProvider = "provideCreateLinkedHashMapFromListsWithIllegalArgs",
            expectedExceptions = IllegalArgumentException.class)
    public void testCreateLinkedHashMapFromListsWithIllegalArgs(final List<String> keys, final List<String> values) {
        FuncotatorUtils.createLinkedHashMapFromLists(keys, values);
    }

}
