package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
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

//    @Test
//    void createRefBaseFile() {
////        printReferenceBases(new File(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()), "chr3", 100000000, 110000000);
////        printReferenceBases(new File("/Users/jonn/Development/references/GRCh37.p13.genome.fasta"), "chr1", 860000,  880000);
////        printReferenceBases();
//    }

    private static Object[] helpCreateDataForTestGetBasesInWindowAroundReferenceAllele(final String refAlelleBases,
                                                                                       final String altAlleleBases,
                                                                                       final String strand,
                                                                                       final int windowSizeInBases,
                                                                                       final int startPos,
                                                                                       final int endPos,
                                                                                       final String expected) {
        return new Object[] {
            Allele.create(refAlelleBases, true),
                Allele.create(altAlleleBases),
                Strand.decode(strand),
                windowSizeInBases,
                new ReferenceContext( refDataSourceHg19Ch3, new SimpleInterval("chr3", startPos, endPos) ),
                expected
        };
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
                        "N"
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
                        "N"
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

        return new Object[][] {
                { new SimpleInterval("chr1", 1, 1),     exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 25, 67),   exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 105, 392), exons_forward, Strand.POSITIVE, -1 },
                { new SimpleInterval("chr1", 10, 10),   exons_forward, Strand.POSITIVE,  1 },
                { new SimpleInterval("chr1", 99, 99),   exons_forward, Strand.POSITIVE, 50 },
                { new SimpleInterval("chr1", 50, 67),   exons_forward, Strand.POSITIVE, 21 },
                { new SimpleInterval("chr1", 67, 75),   exons_forward, Strand.POSITIVE, -1 },

                { new SimpleInterval("chr1", 1, 1),     exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 25, 67),   exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 105, 392), exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 10, 10),   exons_backward, Strand.NEGATIVE, 50 },
                { new SimpleInterval("chr1", 99, 99),   exons_backward, Strand.NEGATIVE,  1 },
                { new SimpleInterval("chr1", 50, 67),   exons_backward, Strand.NEGATIVE, -1 },
                { new SimpleInterval("chr1", 67, 75),   exons_backward, Strand.NEGATIVE, 15 },

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
        };
    }

    @DataProvider
    Object[][] provideDataForGetMitochondrialAminoAcidByCodon() {
        return new Object[][]{
                {null, false, null},
                {"", false, null},
                {"XQZ", false, null},
                {null, true, null},
                {"", true, null},
                {"XQZ", true, null},
                {"ATG", false, AminoAcid.METHIONINE},
                {"CCA", false, AminoAcid.PROLINE},
                {"CCC", false, AminoAcid.PROLINE},
                {"CCG", false, AminoAcid.PROLINE},
                {"CCT", false, AminoAcid.PROLINE},
                {"ATT", false, AminoAcid.ISOLEUCINE},
                {"ATT", true, AminoAcid.METHIONINE},
                {"ATA", false, AminoAcid.METHIONINE},
                {"AGA", false, AminoAcid.STOP_CODON},
                {"AGG", false, AminoAcid.STOP_CODON},
                {"TGA", false, AminoAcid.TRYPTOPHAN},
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
    Object[][] provideDataForGetProteinChangeString() {

        //TODO: Add tests for INDELS, not just ONPs.

        return new Object[][] {
                {"N",   1,  1,  "G", "AAT",    "GGT",    "p.N1G"},
                {"NY",  1,  2, "GN", "AATTAT", "GGTAAT", "p.1_2NY>GN"},
                {"YY",  1,  2, "NY", "TATTAT", "AATTAT", "p.Y1N"},
                {"NY",  1,  2, "NG", "AATTAT", "AATGGT", "p.Y2G"},

                {"N",  71, 71,  "G", "AAT",    "GGT",    "p.N71G"},
                {"NY", 71, 72, "GN", "AATTAT", "GGTAAT", "p.71_72NY>GN"},
                {"YY", 71, 72, "NY", "TATTAT", "AATTAT", "p.Y71N"},
                {"NY", 71, 72, "NG", "AATTAT", "AATGGT", "p.Y72G"},
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
//        expected

//                                                 11111111112222222222333333333344444444445555555555666666
//                                       012345678901234567890123456789012345678901234567890123456789012345
        final String referenceSnippet = "AAATTTGGGCCCATGATATAGGCGCCGTAGCAGTAGATAGCCCCCCAACCGGGGCCCGGGTTTAAA";

        return new Object[][] {

                // alignedRefAlleleStart <= 0
                {referenceSnippet, 10, Allele.create("CC", true), 1, -1, "GCCCAT"},
                {referenceSnippet, 11, Allele.create("CA", true), 1, -1, "CCCATG"},
                {referenceSnippet, 12, Allele.create("AT", true), 1, -1, "CCATGA"},

                // (codingSequenceRefAlleleStart <= 0) && (alignedRefAlleleStart <= 0)
                {referenceSnippet, 10, Allele.create("CC", true), 0, -2, "GCCCAT"},
                {referenceSnippet, 11, Allele.create("CA", true), 0, -2, "CCCATG"},
                {referenceSnippet, 12, Allele.create("AT", true), 0, -2, "CCATGA"},

                {referenceSnippet, 10, Allele.create("CCA", true),  -1, -2, "CCCATG"},
                {referenceSnippet, 11, Allele.create("CAT", true),  -1, -2, "CCATGA"},
                {referenceSnippet, 12, Allele.create("ATG", true),  -1, -2, "CATGAT"},

                // Allele same as reference sequence:
                {referenceSnippet, 0, Allele.create("AAA", true), 1, 1, "AAA"},
                {referenceSnippet, 1, Allele.create("AA", true), 2, 1, "AAA"},
                {referenceSnippet, 2, Allele.create("A", true), 3, 1, "AAA"},

                {referenceSnippet, 6, Allele.create("GGG", true), 1, 1, "GGG"},
                {referenceSnippet, 7, Allele.create("GG", true), 2, 1, "GGG"},
                {referenceSnippet, 8, Allele.create("G", true), 3, 1, "GGG"},

                {referenceSnippet, 6, Allele.create("GGG", true), 60, 60, "GGG"},
                {referenceSnippet, 7, Allele.create("GG", true), 61, 60, "GGG"},
                {referenceSnippet, 8, Allele.create("G", true), 62, 60, "GGG"},

                {referenceSnippet, 17, Allele.create("ATAG", true), 18, 17, "TATAGG"},

                {referenceSnippet, 4, Allele.create(referenceSnippet.substring(4), true), 18, 17, referenceSnippet.substring(3)},

                // Allele diffferent from reference sequence:
                {referenceSnippet, 0, Allele.create("TAA", true), 1, 1, "TAA"},
                {referenceSnippet, 1, Allele.create("AG", true), 2, 1, "AAG"},
                {referenceSnippet, 2, Allele.create("T", true), 3, 1, "AAT"},

                {referenceSnippet, 6, Allele.create("GCC", true), 1, 1, "GCC"},
                {referenceSnippet, 7, Allele.create("AG", true), 2, 1, "GAG"},
                {referenceSnippet, 8, Allele.create("A", true), 3, 1, "GGA"},

                {referenceSnippet, 6, Allele.create("GCC", true), 60, 60, "GCC"},
                {referenceSnippet, 7, Allele.create("AG", true), 61, 60, "GAG"},
                {referenceSnippet, 8, Allele.create("A", true), 62, 60, "GGA"},

                {referenceSnippet, 17, Allele.create("ATAT", true), 18, 17, "TATATG"},

                {referenceSnippet, 4, Allele.create(
                        new String(new char[referenceSnippet.length() - 4]).replace("\0", "A"), true),
                        18, 17,
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
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  1, 500000 + offset, 500000 + offset,          "TAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  2, 500000 + offset, 500000 + offset,         "GTACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  3, 500000 + offset, 500000 + offset,        "AGTACAT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  4, 500000 + offset, 500000 + offset,       "TAGTACATT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  5, 500000 + offset, 500000 + offset,      "TTAGTACATTA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  6, 500000 + offset, 500000 + offset,     "GTTAGTACATTAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  7, 500000 + offset, 500000 + offset,    "AGTTAGTACATTAAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  8, 500000 + offset, 500000 + offset,   "TAGTTAGTACATTAACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+",  9, 500000 + offset, 500000 + offset,  "CTAGTTAGTACATTAACAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "T", "+", 10, 500000 + offset, 500000 + offset, "ACTAGTTAGTACATTAACAAT"),

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  1, 500000 + offset, 500000 + offset,          "TACATT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  2, 500000 + offset, 500000 + offset,         "GTACATTA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  3, 500000 + offset, 500000 + offset,        "AGTACATTAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  4, 500000 + offset, 500000 + offset,       "TAGTACATTAAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  5, 500000 + offset, 500000 + offset,      "TTAGTACATTAACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  6, 500000 + offset, 500000 + offset,     "GTTAGTACATTAACAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  7, 500000 + offset, 500000 + offset,    "AGTTAGTACATTAACAAT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  8, 500000 + offset, 500000 + offset,   "TAGTTAGTACATTAACAATG"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+",  9, 500000 + offset, 500000 + offset,  "CTAGTTAGTACATTAACAATGA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("A", "ATTT", "+", 10, 500000 + offset, 500000 + offset, "ACTAGTTAGTACATTAACAATGAC"),

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  1, 500000 + offset, 500000 + offset,          "TACATT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  2, 500000 + offset, 500000 + offset,         "GTACATTA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  3, 500000 + offset, 500000 + offset,        "AGTACATTAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  4, 500000 + offset, 500000 + offset,       "TAGTACATTAAC"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  5, 500000 + offset, 500000 + offset,      "TTAGTACATTAACA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  6, 500000 + offset, 500000 + offset,     "GTTAGTACATTAACAA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  7, 500000 + offset, 500000 + offset,    "AGTTAGTACATTAACAAT"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  8, 500000 + offset, 500000 + offset,   "TAGTTAGTACATTAACAATG"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+",  9, 500000 + offset, 500000 + offset,  "CTAGTTAGTACATTAACAATGA"),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("ACAT", "A", "+", 10, 500000 + offset, 500000 + offset, "ACTAGTTAGTACATTAACAATGAC"),

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  1, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(         "TAC".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  2, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(        "GTACA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  3, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(       "AGTACAT".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  4, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(      "TAGTACATT".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  5, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(     "TTAGTACATTA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  6, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(    "GTTAGTACATTAA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  7, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(   "AGTTAGTACATTAAC".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  8, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(  "TAGTTAGTACATTAACA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-",  9, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement( "CTAGTTAGTACATTAACAA".getBytes() )) ) ,
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "A", "-", 10, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement("ACTAGTTAGTACATTAACAAT".getBytes() )) ) ,

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  1, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(          "TACATT".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  2, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(         "GTACATTA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  3, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(        "AGTACATTAA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  4, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(       "TAGTACATTAAC".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  5, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(      "TTAGTACATTAACA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  6, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(     "GTTAGTACATTAACAA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  7, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(    "AGTTAGTACATTAACAAT".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  8, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(   "TAGTTAGTACATTAACAATG".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-",  9, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(  "CTAGTTAGTACATTAACAATGA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("T", "TAAA", "-", 10, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement( "ACTAGTTAGTACATTAACAATGAC".getBytes() )) ),

                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  1, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(          "TACATT".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  2, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(         "GTACATTA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  3, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(        "AGTACATTAA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  4, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(       "TAGTACATTAAC".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  5, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(      "TTAGTACATTAACA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  6, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(     "GTTAGTACATTAACAA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  7, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(    "AGTTAGTACATTAACAAT".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  8, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(   "TAGTTAGTACATTAACAATG".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-",  9, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement(  "CTAGTTAGTACATTAACAATGA".getBytes() )) ),
                helpCreateDataForTestGetBasesInWindowAroundReferenceAllele("TGTA", "T", "-", 10, 500000 + offset, 500000 + offset, new String( BaseUtils.simpleReverseComplement( "ACTAGTTAGTACATTAACAATGAC".getBytes() )) ),
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
        Assert.assertEquals(FuncotatorUtils.getAlternateSequence(refCodingSeq, startPos, refAllele, altAllele), expected);
    }

    @Test(dataProvider = "provideDataForGetEukaryoticAminoAcidByCodon")
    void testGetEukaryoticAminoAcidByCodon(final String codon, final AminoAcid expected) {
        Assert.assertEquals(FuncotatorUtils.getEukaryoticAminoAcidByCodon(codon), expected);
    }

    @Test(dataProvider = "provideDataForGetMitochondrialAminoAcidByCodon")
    void testGetMitochondrialAminoAcidByCodon(final String codon, final boolean isFirst, final AminoAcid expected) {
        Assert.assertEquals(FuncotatorUtils.getMitochondrialAminoAcidByCodon(codon, isFirst), expected);
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
                        "Valine"
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

    @Test (dataProvider = "provideDataForGetProteinChangeString")
    void testGetProteinChangeString(final String refAminoAcidSeq,
                                    final int protChangeStartPos,
                                    final int protChangeEndPos,
                                    final String altAminoAcidSeq,
                                    final String refAllele,
                                    final String altAllele,
                                    final String expected) {

        final SequenceComparison seqComp = new SequenceComparison();
        seqComp.setReferenceAminoAcidSequence(refAminoAcidSeq);
        seqComp.setProteinChangeStartPosition(protChangeStartPos);
        seqComp.setProteinChangeEndPosition(protChangeEndPos);
        seqComp.setAlternateAminoAcidSequence(altAminoAcidSeq);
        seqComp.setReferenceAllele(refAllele);
        seqComp.setAlternateAllele(altAllele);

        Assert.assertEquals( FuncotatorUtils.getProteinChangeString(seqComp), expected );
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
                                  final String expected) {
        Assert.assertEquals(
                FuncotatorUtils.getAlignedRefAllele(referenceSnippet,referencePadding,refAllele,codingSequenceRefAlleleStart,alignedRefAlleleStart),
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
    void testGetBasesInWindowAroundReferenceAllele(final Allele refAllele, final Allele altAllele, final Strand strand,
                                                   final int referenceWindow, final ReferenceContext referenceContext,
                                                   final String expected) {

        final String basesInWindow = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(refAllele, altAllele, strand, referenceWindow, referenceContext);
        Assert.assertEquals( basesInWindow, expected );
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
    public Object[][] provideCreateFuncotationsFromVariantContext() {
        final Map<String, String> attributes1 = ImmutableMap.of("FOOFIELD", "FOO", "BAZFIELD", "BAZ");
        final List<VCFInfoHeaderLine> attributes1AsVcfHeaderLine = attributes1.keySet().stream()
                .map(k -> new VCFInfoHeaderLine(k, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Description here"))
                .collect(Collectors.toList());

        final Map<String, String> attributes2 = ImmutableMap.of("FOOFIELD", "FOO,FOOINDEL", "BAZFIELD", "BAZ,BAZINDEL");
        final List<VCFInfoHeaderLine> attributes2AsVcfHeaderLine = attributes2.keySet().stream()
                .map(k -> new VCFInfoHeaderLine(k, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Description here"))
                .collect(Collectors.toList());

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
            }
        };
    }

    @Test(dataProvider = "provideCreateFuncotationsFromVariantContext")
    public void testCreateFuncotationsFromVariantContext(final VariantContext vc, final FuncotationMetadata metadata, final String datasourceName) {
        final List<Funcotation> funcotations = FuncotatorUtils.createFuncotations(vc, metadata, datasourceName);

        Assert.assertTrue(funcotations.stream().allMatch(f -> f.getDataSourceName().equals(datasourceName)));
        Assert.assertEquals(funcotations.stream().map(f -> f.getAltAllele()).collect(Collectors.toSet()), new HashSet<>(vc.getAlternateAlleles()));
        Assert.assertEquals(funcotations.stream().map(f -> f.getMetadata()).collect(Collectors.toSet()), new HashSet<>(Collections.singletonList(metadata)));
    }
}
