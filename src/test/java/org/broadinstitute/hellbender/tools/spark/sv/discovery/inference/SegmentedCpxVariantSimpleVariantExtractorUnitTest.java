package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple3;

import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.b38_reference_chr20_chr21;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.CpxVariantCanonicalRepresentation.UNMAPPED_INSERTION;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.MultiSegmentsCpxVariantExtractor.compactifyMissingSegments;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.MultiSegmentsCpxVariantExtractor.findAllSegments;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class SegmentedCpxVariantSimpleVariantExtractorUnitTest extends GATKBaseTest {

    private static final Random random = new Random(1);

    private static final ZeroAndOneSegmentCpxVariantExtractor zeroAndOneSegmentCpxVariantExtractor = new ZeroAndOneSegmentCpxVariantExtractor();
    private static final MultiSegmentsCpxVariantExtractor multiSegmentsCpxVariantExtractor = new MultiSegmentsCpxVariantExtractor();

    private List<Object[]> caseForZeroAndOneSegmentCalls() {
        final List<Object[]> data = new ArrayList<>(20);

        // NOTE ALL VARIANTS HERE ARE ARTIFICIALLY PUT ON CHR20 AND 21 BECAUSE WE NEED REFERENCE, SO SOME NON-CRITICAL VALUES MAY LOOK NON-SENSE

        // 1. zero segment -> insertion
        VariantContext complex = makeTestComplexVariant(new SimpleInterval("chr20:51740560-51740561"), 549,
                "AT", "AATTAGGCTGGGCACAGTGGCTCACACCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGTGGATCACCTGAGGTCGGGAGTTCAAGACCAGCCTAACCAACATGAGGAAACCCCGTCTCTACTAAAAATACAAAATTAGATGGGCGTGGTGGCGCATGCCTGTAATTCAAACTACTTGGAAGGCTGAGGCAGGAGAATTGCTTGAACCCAGGAGACAGAGGTTGTGGTAAGCCAAGATCATGCCATTGTACTCCAGCATGGGCAACAAGAGTGAGACTCCATCTCAAAAAAAAAAAAAATTAGCCAGGCGTGGTGGTGGGCACCTGTAATCCCAGCTACCCTGGAGACTGAGGCAGAAGAATCGCTTGAACCCAGGAGGCGGAGATTGCAGTGAGCCAAGATTACGCCACTGCACTCCAGCCTGGGCACCAAGAGCAAAACCCTGTCTCAAAAAAATTAACAAATAAAAAGATTTCTGTCTGCCACACGGCTGGTCCATGTGTAAAGACACATTCCTGTTGGTTTTATGTGTCTTGAATTCTAATGGGT",
                Arrays.asList("asm028558:tig00002","asm028558:tig00003"), Collections.emptyList(),
                Arrays.asList("-chr18:11642876-11642927","UINS-496"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                Collections.singletonList(makeInsertion("chr20",51740560, 51740560, 549, Allele.create("A", true))
                        .attribute(CPX_EVENT_KEY, "CPX_chr20:51740560-51740561")
                        .attribute(CONTIG_NAMES, "asm028558:tig00002,asm028558:tig00003")
                        .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                        .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        // 2. one segment -> with inversion
        complex = makeTestComplexVariant(new SimpleInterval("chr21:402806-402905"), 0,
                "GAGTCTTACTCTATTGGGCAGGCTGGAGTACAGCGGTGAAATCATGGCTCACTGCAGCCTCGATGTCCTGGCCTCAAACCATCCCCCTGCTTCAGCCTCC", "GAGGCTGAAGCAGGGGGATGGTTTGAGGCCAGGACATCGAGGCTGCAGTGAGCCATGATTTCACCGCTGTACTCCAGCCTGCCCAATAGAGTAAGACT",
                Collections.singletonList("asm002252:tig00003"), Collections.singletonList(new SimpleInterval("chr21:402807-402904")),
                Collections.singletonList("-1"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                Collections.singletonList(makeInversion(new SimpleInterval("chr21:402807-402904"), Allele.create("N", true)) // THE REF ALLELE N HERE IS BECAUSE OF COORDINATE MESSING on TEST DATA MENTIONED ABOVE
                        .attribute(CPX_EVENT_KEY, "CPX_chr21:402806-402905")
                        .attribute(CONTIG_NAMES, "asm002252:tig00003")
                        .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                        .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        // 3. one segment -> when deletion is not allowed

        // 3.0 -> new material is not long enough (49 bp, boundary cases)
        complex = makeTestComplexVariant(new SimpleInterval("chr21:402806-402905"), 0,
                "GAGTCTTACTCTATTGGGCAGGCTGGAGTACAGCGGTGAAATCATGGCTCACTGCAGCCTCGATGTCCTGGCCTCAAACCATCCCCCTGCTTCAGCCTCC", "AAAAAAAAAAAAAAAAAAAAAAAAGAGTCTTACTCTATTGGGCAGGCTGGAGTACAGCGGTGAAATCATGGCTCACTGCAGCCTCGATGTCCTGGCCTCAAACCATCCCCCTGCTTCAGCCTCCAAAAAAAAAAAAAAAAAAAAAAAAA",
                Collections.singletonList("asm002252:tig00003"), Collections.singletonList(new SimpleInterval("chr21:402807-402904")),
                Arrays.asList("UINS-24", "1", "UINS-25"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                Collections.emptyList()
        });

        // 3.1 -> ...., 1
        complex = makeTestComplexVariant(new SimpleInterval("chr20:18675721-18675877"), 408,
                "TATGTGTATATTTACACACATATATATGTAAATATACCTATGTGTATATTTACACACATATATATGTGTAAATATACCTATGTGTATATTTACACACATATATGTAAATATACCTATGTGTATATTTACACATATATATGTAAATATACCTATGTGT", "TATGTGTATATTTACACACATATATATGTAAATATACCTATGTGTATATTTACACACATATATATGTGTAAATATACCTATGTGTATATTTACACACATATATATGTGTAAATATACCTATGTGTATATTTACACATATATATGTAAATATACCTATGTGTATGTTTACACATATATATGTAAATATACCTATGTGTATGTTTACACATATATATGTGTAAATATACCTATGTGTATGTTTACACATATATATGTGTAAATATACCTATGTGTATGTTTACACATATATGTAAATATACCTATGTGTATGTTTACACATATATGTGTAAATATACCGATGTGTATGTTTACACATATATGTGTAAATATACCTATGTGTATGTTTACACATATATGTGTAAATATACCTATGTGTATGTTTACACATATATATGTGTAAATATACCTATGTGTATGTTTACACATATATATGTGTAAATATACCTATGTGTGTGTTTACACATATATATGTGTAAATATACCTATGTGTGTGTTTACACATATATATGTAAATATACCTATGTGT",
                Collections.singletonList("asm028012:tig00004"), Collections.singletonList(new SimpleInterval("chr20:18675721-18675877")),
                Arrays.asList("1","UINS-28","1","UINS-64","1"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                Collections.singletonList(makeInsertion("chr20", 18675720, 18675720, 408, Allele.create("A", true))
                        .attribute(CPX_EVENT_KEY, "CPX_chr20:18675721-18675877").attribute(CONTIG_NAMES, "asm028012:tig00004")
                        .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                        .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        // 3.2 -> 1, .....
        complex = makeTestComplexVariant(new SimpleInterval("chr20:64096905-64097041"), 318,
                "CCACCATCATCACCATCACCACTATCACCACCACCACCATCATTACCATCATCATCACGACCATCACCACCATCATCACCATCACCACTATCACCACCACCACCATCATCACCATCACCATCATCACAGTCATCACC", "CCACCATCATCACCATCACCACTATCACCACCACCACCATCATCACCATCACCATCATCACAGTCATCACTGTCACCATCATCACCATCCTCACTATCACCACCACCACCATCATCACCATCACCATCATCACAGTCATCACCGCCACCATCATCACCATCCTCACTATCACCACCACCACCATCATCACCATCACCATCATCACAATCATCACCGTCACCATCATCACCATCCTCACTATCACCACCACCACCATCATCACCATCACTATCATCACAGTCATCACCGTCACCATCATCACCATCCTCACTATCACCACCACTACCATCATCATCACATTCATCATCACTATTACCATCATCATCACCACCATCACCATCACTATCACCACCATCATTACATTTGTCACCATCACCACCATTATCACCATCACCGCTATCACCACCACCACCGTC",
                Collections.singletonList("asm028821:tig00001"), Collections.singletonList(new SimpleInterval("chr20:64096905-64097041")),
                Arrays.asList("1","1","UINS-166"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                Collections.singletonList(makeInsertion("chr20", 64097041, 64097041, 318, Allele.create("A", true))
                        .attribute(CPX_EVENT_KEY, "CPX_chr20:64096905-64097041").attribute(CONTIG_NAMES, "asm028821:tig00001")
                        .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                        .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        // 3.3 -> ...., 1, ....
        complex = makeTestComplexVariant(new SimpleInterval("chr20:51740560-51741035"), 599,
                "ATTTTGTGTTGTTGTTTTTGTTTTTTGAGACAAGGTCTCATTCTGTCACCCAGGCAGGACTGTGGTGGCACCATCATGGCTCAGCGCAGCCTCCTTTTCCCCAGGCTCAAGTGATCCTCTTGCCTCAGCCTCCCACGTGGCTGGGACTACAGGTGTGTACCACCACTCCCGGATAATTTTTTTTATTTTTTATTTTTAGTAAAGACAGTCTCACTATGTTGCCCAGGCTGGTCTCCAACTCCTGGTCTCAAGCAATCCTCCCAGTTCAGCCTCTCAAAGTGCTGGGATTACAGATGTGAGCCACAATACCCGGCCCCAATTCTAATGTTTAAAGAGTACAGTCTACACCTTAAAGCCTGCATTTTATCATCCTGTCCTCACTGCTCTGACTTCTTTACAGTTGTGCTGTCCACCTTGGCGGCTTCTACCACATGTGGCTATTTTAAGTTTCAATTAATTAAAATTAAATTTTAATT", "AATTAGGCTGGGCACAGTGGCTCACACCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGTGGATCACCTGAGGTCGGGAGTTCAAGACCAGCCTAACCAACATGAGGAAACCCCGTCTCTACTAAAAATACAAAATTAGATGGGCGTGGTGGCGCATGCCTGTAATTCAAACTACTTGGAAGGCTGAGGCAGGAGAATTGCTTGAACCCAGGAGACAGAGGTTGTGGTAAGCCAAGATCATGCCATTGTACTCCAGCATGGGCAACAAGAGTGAGACTCCATCTCAAAAAAAAAAAAAATTAGCCAGGCGTGGTGGTGGGCACCTGTAATCCCAGCTACCCTGGAGACTGAGGCAGAAGAATCGCTTGAACCCAGGAGGCGGAGATTGCAGTGAGCCAAGATTACGCCACTGCACTCCAGCCTGGGCACCAAGAGCAAAACCCTGTCTCAAAAAAATTAACAAATAAAAAGATTTCTGTCTGCCACACGGCTGGTCCATGTGTAAAGACACATTCCTGTTGGTTTTATGTGTCTTGAATTCTAATGGGTTTTGTGTTGTTGTTTTTGTTTTTTGAGACAAGGTCTCATTCTGTCACCCAGGCAGGACTGTGGTGGCACCATCATGGCTCAGCGCAGCCTCCTTTTCCCCAGGCTCAAGTGATCCTCTTGCCTCAGCCTCCCACGTGGCTGGGACTACAGGTGTGTACCACCACTCCCGGATAATTTTTTTTATTTTTTATTTTTAGTAAAGACAGTCTCACTATGTTGCCCAGGCTGGTCTCCAACTCCTGGTCTCAAGCAATCCTCCCAGTTCAGCCTCTCAAAGTGCTGGGATTACAGATGTGAGCCACAATACCCGGCCCCAATTCTAATGTTTAAAGAGTACAGTCTACACCTTAAAGCCTGCATTTTATCATCCTGTCCTCACTGCTCTGACTTCTTTACAGTTGTGCTGTCCACCTTGGCGGCTTCTACCACATGTGGCTATTTTAAGTTTCAATTAATTAAAATTAAATTTTAATTTAATTAATTAAAAATAAATTTTAATTAATTAATTAAAAATAAATTTTAAT",
                Arrays.asList("asm028558:tig00000", "asm028558:tig00001"), Collections.singletonList(new SimpleInterval("chr20:51740561-51741034")),
                Arrays.asList("-chr18:11642876-11642927","UINS-496","1","UINS-49"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                Arrays.asList(
                        makeInsertion("chr20", 51740560, 51740560, 549, Allele.create("A", true))
                                .attribute(CPX_EVENT_KEY, "CPX_chr20:51740560-51741035").attribute(CONTIG_NAMES, "asm028558:tig00000,asm028558:tig00001")
                                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make(),
                        makeInsertion("chr20", 51741034, 51741034,  50, Allele.create("T", true))
                                .attribute(CPX_EVENT_KEY, "CPX_chr20:51740560-51741035").attribute(CONTIG_NAMES, "asm028558:tig00000,asm028558:tig00001")
                                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        // 4. one segment -> whenNoInvAndNoAsIsAppearance

        // 4.1 -> deletion but no insertion
        complex = makeTestComplexVariant(new SimpleInterval("chr20:20269131-20269199"), -34,
                "ATATATATATATATATACACACACACACACACACATACATATATGTATATACACACACATATATACATA", "ACACACACACACACACACACACACACACACACACA",
                Collections.singletonList("asm028026:tig00000"), Collections.singletonList(new SimpleInterval("chr20:20269131-20269199")),
                Collections.singletonList("-chrX:137700299-137700331"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                              Collections.singletonList(makeDeletion(new SimpleInterval("chr20:20269131-20269198"), Allele.create("A", true))
                                      .attribute(CPX_EVENT_KEY, "CPX_chr20:20269131-20269199").attribute(CONTIG_NAMES, "asm028026:tig00000")
                                      .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                                      .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });
        // 4.2 -> deletion and insertion
        complex = makeTestComplexVariant(new SimpleInterval("chr20:54849491-54849615"), 15,
                "CAAATCTCATGTGAAATGTATCCCCAGTGTGGAGGGGGCAGATCCTCATAATGGCTTGGGCCCTTCCATGGTAATAGTGAGTCTTGCTCTGTAGTTCATAGAGAGCTGATTGTTAAAGGAGTCTG", "CCAAATCTCATGTTGAAATGTAATCCCCAGTGTTGGAGGGGGGCAGATCCCTCATGAATGGCTTGGTGCCCTTCCCATGGTAATGAGTGAGTTCTTGCTCTGTTAGTTCATGAGAGAGCTGATTGTTTAAAGGAGTCTGG",
                Collections.singletonList("asm028586:tig00000"), Collections.singletonList(new SimpleInterval("chr20:54849491-54849615")),
                Arrays.asList("UINS-36","-chr14:58474127-58474172","UINS-54"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                 Arrays.asList(
                         makeDeletion(new SimpleInterval("chr20:54849491-54849614"), Allele.create("C", true))
                                 .attribute(CPX_EVENT_KEY, "CPX_chr20:54849491-54849615").attribute(CONTIG_NAMES, "asm028586:tig00000")
                                 .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                                 .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make(),
                         makeInsertion("chr20", 54849491, 54849491, 140, Allele.create("c", true))
                                 .attribute(CPX_EVENT_KEY, "CPX_chr20:54849491-54849615").attribute(CONTIG_NAMES, "asm028586:tig00000")
                                 .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                                 .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        // 4.3 -> fat insertion
        complex = makeTestComplexVariant(new SimpleInterval("chr20:12558793-12558810"), 133,
                "AAAAAAAAAAAAAAAAAA", "AGACAAAGAAACAAACAAACAAAACAAAACTATATATATATATATATACACACACACACACACACACACATTATTAAAATTCAGATTTAAATAAACTGACTATAAAAAAGTACTTTTGAAACAAAAACTTTAATCATGATTATATATATTA",
                Collections.singletonList("asm027960:tig00003"), Collections.singletonList(new SimpleInterval("chr20:12558793-12558810")),
                Arrays.asList("-chrX:99014092-99014129","UINS-101"));
        data.add(new Object[]{complex, b38_reference_chr20_chr21, zeroAndOneSegmentCpxVariantExtractor,
                 Collections.singletonList(makeInsertion("chr20", 12558793, 12558809, 133, Allele.create("AAAAAAAAAAAAAAAAA", true))
                         .attribute(CPX_EVENT_KEY, "CPX_chr20:12558793-12558810").attribute(CONTIG_NAMES, "asm027960:tig00003")
                         .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                         .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make())
        });

        return data;
    }
    private List<Object[]> caseForMultiSegmentsCalls() {
        final List<Object[]> data = new ArrayList<>(20);

        // case 1: long stretch

        // case 1.1: front insertion only
        SimpleInterval affectedInterval = new SimpleInterval("chr21:21264944-21265096");
        int svLen = 215;
        String refAllele = "T";
        String altSeq = "TATATATATGTGTATGTGTATATATACACATATATATTATATATGTGTATGTGTATATATACACATATATATTATATATGTGTATGTGTATATATACACATATATATTATATATGTGTATGTGTATATATATACACATATATATTATATATGTGTATGTGTATATATATACACATATATATTATATATGTGTATATGTATATATACACATATATATTATATATATATGTGTCTGTATATATATACACATATATATTATATATATGTGTCTGTGTATATATATACACATATATATGTGTCTGTGTATATATGTACACATATATACTATATATGTGTATGTGTATATATATACACACATATATTATATAT";
        List<String> ctgNames = Arrays.asList("asm029034:tig00000","asm029034:tig00001");
        List<SimpleInterval> refSegments = Arrays.asList(new SimpleInterval("chr21:21264944-21264988"), new SimpleInterval("chr21:21264988-21265052"), new SimpleInterval("chr21:21265052-21265096"));
        List<String> altArrangements = Arrays.asList("1","2","3","2","1","2","3");
        VariantContext complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        List<VariantContext> expectedSimple = Collections.singletonList(makeInsertion("chr21", 21264943, 21264943, 221, Allele.create("G", true)).attribute(CPX_EVENT_KEY, "CPX_chr21:21264944-21265096").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 1.2: rear insertion only
        affectedInterval = new SimpleInterval("chr20:61919906-61920109");
        svLen = 541;
        refAllele = "T";
        altSeq = "TCGTGATTATGTGGAAGCGTGGTGTCACGGTGATTGCGTGGAAGCGTGTTGTGATTGTGTGGAAGCGTGGTATCGCGGTGATTGCATGGAAGTGTGGTGTCACAGTGATTGCGTGGAAGCGTGTCGTGATTGTGTGGAAGCATGGTATCGTGATTGTGTGGAAGCGTGGTGTCACGGTGATTGCGTGGAAGCATGTTGTGATTGTGTGGAAGCGTGGTGTGATTGTGTGGAAGCATGGTATCGTGATTGTGGAAGCGTGGTATCGCGGCGATTGTGTGGAAGCGTGGTGTCGCAGTGATTGCGTGGAAGCATGTTGTGATTGTGTGGAAGCGTGGTATCGTGATTGTGTGGAAGCATGGTGTCGTGATTGTGTGGAAGCATGTCGTGATTGTGTGGAAGTGTGATGTCACGGTGATTGCGTGGAAGCGTGTTGTGATTGTGTGGAAGCGTGGTATCGTGATTGGAAGTGTGGTGTCACGCTGATTGCATGGAAGTGTGTTGTGATTGTGTGGAAGCGTGATATCGCAGTGATTGTGTGGAAGCGTGGTGTCACGGTGATTGTGTGGAAGCGTGGTGTCACGGTGATTGCGTGGAAGCGTGTTGTGATTGTGTGGAAGCGTGGTATCGTGATCGGAAGCGTGGTGTTGCGGTGATTGCATGGAAGCATGTTGTGATTGTGTGGAAGCATGGTATCGTGATTGTCTGGAAGCATGGTGTCATGGTGATTGGAAGTGTGTCGTGATTG";
        ctgNames = Arrays.asList("asm028707:tig00000");
        refSegments = Arrays.asList(new SimpleInterval("chr20:61919906-61919908"), new SimpleInterval("chr20:61919908-61920054"), new SimpleInterval("chr20:61920054-61920109"));
        altArrangements = Arrays.asList("1","2","3","UINS-177","1","2","2","3");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Collections.singletonList(makeInsertion("chr20", 61920109, 61920109, 531, Allele.create("G", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr20:61919906-61920109").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 1.3: front and rear insertion
        affectedInterval = new SimpleInterval("chr20:38653054-38653283");
        svLen = 485;
        refAllele = "T";
        altSeq = "TGGTGGTGGTGGTGATGGAAATGATGATGTTAGTTGTGGTGTTGATGATGGTAATGATAATGATTGTGATGGTGGTGTTGGTGGTGCTGGTGATGATAATGATGGTGGTGGTGGTGGTGGTGATGATGGTGGTGGTGGTGATGGAAATGATGATGATGTTAGTTGTGGTGTTGATGATGGTAATGATAATGATTGTGATGGTGGTGTTGGTGGTGGTGATAATGATGGTAGTGGTGGTGGTGATGATGGTGATGATGATTATGATGGTGTGTTGGTGGTGCTGGTGATGATAATGGTGGTGGTGGTGGTGGTGATGGAAATGATGATGATGTTAATTGTGGTGTTGATGATGGTAATGATAATGATTGTGATGGTGGTGTTGGTGGTGCTGGTGATGATAATGGTGGTGGTGGTGGTGATGGTGATGATGATTATGATGGTGGTGGTGGTGGTGGTGGTGCTGGTGATAGTGGTGGTGGTGGTGCTGGTGATGATAATGGTGGTGGTGGTGATGATGGTGATGATGATTATGATGGTGGTGTTGGTGGTGCTGGTGATGATAATCATGCTGGTGGTGGTGGCGTTGATGATGGTGACAGTAGTGGTGATGATGGTGGTGGTGGTGATGGAAATGATGATGATGTTAGTTGTGGTGTTGATGATGGTAATGATAATGATTGTGATGATGGTGGTGGTGGTGGTGATGATGGTGATG";
        ctgNames = Arrays.asList("asm028418:tig00000");
        refSegments = Arrays.asList(new SimpleInterval("chr20:38653054-38653113"), new SimpleInterval("chr20:38653113-38653145"), new SimpleInterval("chr20:38653145-38653179"), new SimpleInterval("chr20:38653179-38653273"), new SimpleInterval("chr20:38653273-38653283"));
        altArrangements = Arrays.asList("1","2","3","4","3","1","2","3","4","5","2","3","4","5");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Arrays.asList(
                makeInsertion("chr20", 38653053, 38653053, 259, Allele.create("A", true))
                        .attribute(CPX_EVENT_KEY, "CPX_chr20:38653054-38653283").attribute(CONTIG_NAMES, ctgNames)
                        .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                        .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make(),
                makeInsertion("chr20", 38653283, 38653283, 175, Allele.create("G", true))
                        .attribute(CPX_EVENT_KEY, "CPX_chr20:38653054-38653283").attribute(CONTIG_NAMES, ctgNames)
                        .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                        .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 2: possibly inversion

        // case 2.1: both as-is and inverted, hence no inversion
        affectedInterval = new SimpleInterval("chr20:23122561-23122996");
        svLen = -293;
        refAllele = "C";
        altSeq = "CTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGC";
        ctgNames = Arrays.asList("asm028059:tig00000","asm028059:tig00001");
        refSegments = Arrays.asList(new SimpleInterval("chr20:23122561-23122596"), new SimpleInterval("chr20:23122596-23122666"), new SimpleInterval("chr20:23122666-23122996"));
        altArrangements = Arrays.asList("1","2","-1");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Collections.singletonList(makeDeletion(new SimpleInterval("chr20:23122666-23122995"), Allele.create("C", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr20:23122561-23122996").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 2.2: no as-is but with inverted, but too short (coordinate and allele massage, original event CPX_chr20:34732145-34733344)
        affectedInterval = new SimpleInterval("chr20:34732145-34733344");
        svLen = -1139;
        refAllele = "G";
        altSeq = "GTTGTCCAGGTGGATCTTGAGTATTTTTAGTAGAGACGGGGGGTTCAATTAACTCTTCCAA";
        ctgNames = Arrays.asList("asm010456:tig00000");
        refSegments = Arrays.asList(new SimpleInterval("chr20:34732145-34733303"), new SimpleInterval("chr20:34733303-34733342"), new SimpleInterval("chr20:34733342-34733344"));
        altArrangements = Arrays.asList("-3","-2","UINS-14","3"); // segment 1 deleted, segment 2 appear inverted but length too short
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Collections.singletonList(makeDeletion(new SimpleInterval("chr20:34732145-34733302"), Allele.create("A", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr20:34732145-34733344").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 2.3: no as-is but with inverted, and inverted sequence long enough
        affectedInterval = new SimpleInterval("chr21:26001843-26002386");
        svLen = -1;
        refAllele = "A";
        altSeq = "ACAATGTCAATACAAATGGGGACTAACTATTGTTTTACTTCCCTATACACGTTGCTCTTTCAGTAGATGCAATAAAGTACTGGTAAAACCAGAGGTGGCTACCATCACGATGATGTCAACAGGAGGGACAGTCAGCACTAAGCCCAGAAGGTGTCAAACACTCCGCAGGAGAAATGCGCCATGCAACGGACATGAAGATGATCTGACACTCTTCACGTGGTTTTCAGATGGAAACGTGGCTACGAAAGCATCAACCTCATTATCATCCATCATTAAGGCCATCTCACTCAGTACTGCTGCTTTCAAAGTCCACGCTCCCAAAGCAAATTGGATTTCTGTACACAATACTCTTACAGGATGAAACCCAACCAACTTGTTGACGTAAGTATCCCATCTACTTCGCAATTTTATTTATCTTCCAAAATTAAAGGACTGGCACCCTGATTTATTAAAAGTGAATTGGTTCTAGGGACCATATCCCCTCTGAGTTACTGACAGAGCAGCTTCTGGCCTGTGAAGCTCAAAGCCATGCCTAGATGTGAG";
        ctgNames = Arrays.asList("asm029075:tig00000");
        refSegments = Arrays.asList(new SimpleInterval("chr21:26001844-26002384"), new SimpleInterval("chr21:26002384-26002386"));
        altArrangements = Arrays.asList("-1");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Collections.singletonList(makeInversion(new SimpleInterval("chr21:26001844-26002384"), Allele.create("T", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr21:26001843-26002386").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 3: possibly deletion

        // case 3.1: deleted range too short

        affectedInterval = new SimpleInterval("chr21:23428920-23429023");
        svLen = 131;
        refAllele = "T";
        altSeq = "TATAAATATATATAATAATATATAATTATATATTATATTATATAATATAATAATATATATTATAATACATTATATAATATATTATAATATATATAATATAATATATTATATTATATAATATATATTACATAATATATTATATTATATAATATATATTACATAATATATTATATTATATAATATATATTACATAATATATTATATTATATAATATATATTACATAATATATTATAT";
        ctgNames = Arrays.asList("asm029052:tig00000","asm029052:tig00001");
        refSegments = Arrays.asList(new SimpleInterval("chr21:23428920-23428968"), new SimpleInterval("chr21:23428968-23428998"), new SimpleInterval("chr21:23428998-23429023"));
        altArrangements = Arrays.asList("UINS-84","2","3","UINS-5","2","2","3");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Collections.singletonList(makeInsertion("chr21", 23428920, 23428920, 85, Allele.create("T", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr21:23428920-23429023").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 3.2: deleted range long enough (tested together with 2.1)

        // case 4: not long stretch but front and rear insertion possible

        // case 4.1: insertion length not big enough (coordinate and allele massage, original event CPX_chr10:13062977-13063278)
        affectedInterval = new SimpleInterval("chr20:13062977-13063278");
        svLen = 72;
        refAllele = "C";
        altSeq = "CTGACCAATCAGCACCCTTGGCTCACTGGCTTACCGATTTCATCTCTGACCAATCAGCACTACTTGCCCACTGCCTAGACAGAGCTGATAAATCAAGACAGGGGAACTGCAATAGAGAAAGAGTAATTCACACAGAGCCGGCTGTGTGGGAGACCAGGGTTTTGTTATTACTCAAATCAGTCTCCCCGAACATTCGGGGAGCAAAGTTTTTAAAAATAACTTGGTGGGTGGGGGGTAAGCCAGTGAGCCAGGAGTGCTGATTGGTCAGAGATGAAATTGGTAAGCCAGTGATCCAGGAGTGCTGATTGGTCAGAGATGAAATCGGTAAGCCAGTGAGCCAGGAGTGCTGATTGGTCAGCACCCTTGGTAACCAC";
        ctgNames = Arrays.asList("asm016524:tig00000");
        refSegments = Arrays.asList(new SimpleInterval("chr20:13062977-13063037"), new SimpleInterval("chr20:13063037-13063272"), new SimpleInterval("chr20:13063272-13063278"));
        altArrangements = Arrays.asList("1","-2","-1","UINS-14");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Arrays.asList(makeInversion(new SimpleInterval("chr20:13063037-13063272"), Allele.create("G", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr20:13062977-13063278").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 4.2: front insertion only (tested together with 3.1)

        // case 4.3: rear insertion only (coordinate and allele massage, original event CPX_chr22:36680290-36680686)
        affectedInterval = new SimpleInterval("chr21:36680290-36680686");
        svLen = 318;
        refAllele = "T";
        altSeq = "TTCATGATCACTGCCACCTGCATCATGGCTATAATTACCACCGCCATCACCATCATCACTACCACCATCATCACCATCATCACCATTACCACCACCATCACTGCCACCATCACTACCACCATCATCACCATCATCACCACCACCACCATCACCACCATCACCATCACCACTTTCATCACCACCATCTTTATCACAGTCATTATTACCACCATCAATCATCACCACCTTCATGATCACTGCCACCTGCATCATGGTTACAATTACTACCACCACCATCAGCACCACCTTCATGATCACCACCACCTGCATCATGGCTATAATTACTACCACCATCACCACCACTAACACCACCATCATTATCACCACCATCTTCATGATCACTGCCACCTGCATCATGGCTATAATTACCACCATCATCACCACTATACTACCACCACCATCACCACAACCATCGCTACCACCACCACCACCACCATCACCATCATCACCATCACTACTACTGCCACCACTACCAAAACCACCACCACCACCATCACCACCACCATTGCCACCACTACCATCACCACCATCACCACCATCACCATCACCACCACTACCACCATTACCACCACCACCACCACCACCACCACCATCACCACCACCACCACCATCATCACCACTATC";
        ctgNames = Arrays.asList("asm029759:tig00000","asm029759:tig00001");
        refSegments = Arrays.asList(new SimpleInterval("chr21:36680290-36680331"), new SimpleInterval("chr21:36680331-36680659"), new SimpleInterval("chr21:36680659-36680686"));
        altArrangements = Arrays.asList("1","2","1","UINS-249");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Arrays.asList(makeInsertion("chr21", 36680686, 36680686, 250, Allele.create("A", true))
                .attribute(CPX_EVENT_KEY, "CPX_chr21:36680290-36680686").attribute(CONTIG_NAMES, ctgNames)
                .attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0))
                .attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        // case 4.4: front and rear insertion (coordinate and allele massage, original event CPX_chr19:8888822-8895655)
        affectedInterval = new SimpleInterval("chr20:8888822-8895655");
        svLen = -6179;
        refAllele = "G";
        altSeq = "GTTTAGAGGGTCAAGGCGGTGAGTGCAGATGGTGTCCACGCCGGTGGCTGCCCCACGTTTTTCAGGCCTGAGGAGGATAAGTGAGGGGAATGACAATAAGTGGATTGCCTAGGAAGAGGTCCTGAGATTCCTGAGGGAAGGACAGGCAGGGGCCAGAGGAAGTGGTACAGAGGTGAGAAGCTACACATAGATTCTGTGAGTGCACCAACTCCCAGCCAGCCTGGGTCTATGAGCAGAATCATGGACACAAAAATATTTGTAACATTATAAAGTGGATCCTCTGTTCTCAGAAGATGGGGACAGGAGAGCTCAGACTTGAGACAGGTGAAGAGCTGGAATTGGGGATCAAAGCAAGTGGATGAGGCATGTTTGTGTTTGAAAGAAATGCTGGTTGTAGAGTTTATTATCTGCAGGAAAGGGAGGTGGTGGGAAAGGAGGGAACTTAGAGGGATGACCAGGGTAGGGGCCCTATCCTGCTGATGAGCTCTGCAGCATCACTGTTGAGTGGAGAGATGTTACTATCAGGGACACTTTGGTCCTGGAGCCACTGCCTCCTGGATTCCACCAGTGCTGAGGGCACCCCTAGGGAGGGCAGGAAGGAAGTTTTCCAAACTTCTTGGGACTGGGAACAATTGGAGAAGGTAGGCTGGCTCCT";
        ctgNames = Arrays.asList("asm026939:tig00007");
        refSegments = Arrays.asList(new SimpleInterval("chr20:8888822-8895288"), new SimpleInterval("chr20:8895288-8895361"), new SimpleInterval("chr20:8895361-8895655"));
        altArrangements = Arrays.asList("UINS-297","2","UINS-280");
        complex = makeTestComplexVariant(affectedInterval, svLen, refAllele, altSeq, ctgNames, refSegments, altArrangements);
        expectedSimple = Arrays.asList(makeDeletion(new SimpleInterval("chr20:8888822-8895287"), Allele.create("G", true)).attribute(CPX_EVENT_KEY, "CPX_chr20:8888822-8895655").attribute(CONTIG_NAMES, ctgNames).attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0)).attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make(),
                                       makeDeletion(new SimpleInterval("chr20:8895361-8895654"), Allele.create("T", true)).attribute(CPX_EVENT_KEY, "CPX_chr20:8888822-8895655").attribute(CONTIG_NAMES, ctgNames).attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0)).attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make(),
                                       makeInsertion("chr20", 8888822, 8888822, 298, Allele.create("G", true)).attribute(CPX_EVENT_KEY, "CPX_chr20:8888822-8895655").attribute(CONTIG_NAMES, ctgNames).attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0)).attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make(),
                                       makeInsertion("chr20", 8895655, 8895655, 281, Allele.create("C", true)).attribute(CPX_EVENT_KEY, "CPX_chr20:8888822-8895655").attribute(CONTIG_NAMES, ctgNames).attribute(MAX_ALIGN_LENGTH, complex.getAttributeAsInt(MAX_ALIGN_LENGTH, 0)).attribute(MAPPING_QUALITIES, complex.getAttributeAsString(MAPPING_QUALITIES, "")).make());
        data.add(new Object[]{complex, b38_reference_chr20_chr21, multiSegmentsCpxVariantExtractor, expectedSimple});

        return data;
    }
    @DataProvider(name = "forTestSegmentedCpxVariantExtractor")
    private Object[][] forTestSegmentedCpxVariantExtractor() {
        final List<Object[]> data = caseForZeroAndOneSegmentCalls();
        data.addAll(caseForMultiSegmentsCalls());
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestSegmentedCpxVariantExtractor")
    public void testSegmentedCpxVariantExtractor(final VariantContext complexVC, final ReferenceMultiSparkSource reference,
                                                 final SegmentedCpxVariantSimpleVariantExtractor worker,
                                                 final List<VariantContext> expected) {
        assertVariantsAreEqual(worker.extract(complexVC, reference), expected, Collections.emptyList(), b38_seqDict_chr20_chr21);
    }

    @DataProvider(name = "forTestMultiSegmentsCpxVariantExtractorFindAllSegments")
    private Object[][] forTestMultiSegmentsCpxVariantExtractorFindAllSegments() {
        final List<Object[]> data = new ArrayList<>(20);

        List<String> altArrangement = Arrays.asList("UINS-58", "1", "2", "2");
        data.add(new Object[]{altArrangement, 2, 1});
        data.add(new Object[]{altArrangement, 3, -1});

        altArrangement = Arrays.asList("UINS-58","1","2","3","1","chrX:10000-10200","1","3","3","UINS-15");
        data.add(new Object[]{altArrangement, 3, 1});
        altArrangement = Arrays.asList("UINS-58","1","2","3","1","chrX:10000-10200","1","2","3","UINS-15");
        data.add(new Object[]{altArrangement, 3, 6});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestMultiSegmentsCpxVariantExtractorFindAllSegments")
    public void testMultiSegmentsCpxVariantExtractorFindAllSegments(final List<String> altArrangement, final int segmentCount,
                                                                    final int expectedIdx) {
        Assert.assertEquals(findAllSegments(altArrangement, segmentCount),
                            expectedIdx);
    }

    @DataProvider(name = "forTestMultiSegmentsCpxVariantExtractorCompactifyMissingSegments")
    private Object[][] forTestMultiSegmentsCpxVariantExtractorCompactifyMissingSegments() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{Sets.newHashSet(new SimpleInterval("chr1:10000-10010"), new SimpleInterval("chr1:10012-10020")),
                              Arrays.asList(new SimpleInterval("chr1:10000-10010"), new SimpleInterval("chr1:10012-10020"))});

        data.add(new Object[]{Sets.newHashSet(new SimpleInterval("chr1:10000-10010"), new SimpleInterval("chr1:10011-10020")),
                              Arrays.asList(new SimpleInterval("chr1:10000-10020"))});

        data.add(new Object[]{Sets.newHashSet(new SimpleInterval("chr1:10000-10010"), new SimpleInterval("chr1:10010-10020")),
                              Arrays.asList(new SimpleInterval("chr1:10000-10020"))});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestMultiSegmentsCpxVariantExtractorCompactifyMissingSegments")
    public void testMultiSegmentsCpxVariantExtractorCompactifyMissingSegments(final Set<SimpleInterval> missingSegments,
                                                                              final List<SimpleInterval> expected) {
        Assert.assertEquals(compactifyMissingSegments(missingSegments),
                            expected);
    }

    //==================================================================================================================

    @DataProvider(name = "forTestIsConsistentWithCPX")
    private Object[][] forTestIsConsistentWithCPX() {
        final List<Object[]> data = new ArrayList<>(20);

        List<SimpleInterval> refSegments = Arrays.asList(new SimpleInterval("chr1:4939507-4939535"), new SimpleInterval("chr1:4939535-4939614"));
        List<String> altArrangements = Arrays.asList("UINS-58", "1", "2", "-2");
        VariantContext complex = makeTestComplexVariant(new SimpleInterval("chr1:4939506-4939614"), 109, "ATATATAT", "TTTTTTTTTTTTTTTTTTTTT", Collections.singletonList("dummy"),
                refSegments, altArrangements);
        RelevantAttributes relevantAttributes = new RelevantAttributes(complex);
        data.add(new Object[]{makeInversion(new SimpleInterval("chr1:4939535-4939614"), Allele.create("A", true)).make(), relevantAttributes,
                              false});
        data.add(new Object[]{makeInsertion("chr1", 4939507, 4939507, 59, Allele.create("A", true)).make(), relevantAttributes, true});
        data.add(new Object[]{makeDeletion(new SimpleInterval("chr1:4939535-4939614"), Allele.create("A", true)).make(), relevantAttributes, false});

        complex = makeTestComplexVariant(new SimpleInterval("chr6:857170-857852"), -477,
                "CTCTCTTCAGAGGAAATAAATTAAAATATACTAATTGTGTTAGAAAAGCCTAAACCTTAAAATTCAATATAATTGTGGTCAAATAATGCAGATTATGAAATGTGCATGTGAGAGTCCTAGCTCAAGGAGAAGTCGTGCCTCAGTCACTGCCACTCGTCCACCCATCCATCCTCCTATCCGCTCATTTATCCATCCATCCACCTATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCTCCTATCCTCTCATTTATCCATCCATCCACCTATCCACCCATCCATCCATCCTCCTATCCTCTCATTTATCCATCCATCCACCTATCCATTCCTCCATCCATCCTCCTCTCATTTATCCATCCATCCACCTATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCTCCTGTCCTCTCATTTATCCATCCATCCACCTATCCACCCATCCATCCATCCTCCTATCCTCTCATTTATCCATCCATCCACCTAGTCACCCATCCATCCATCCTCTTATCCTCTCATTTATCCATCCATCCACCTATCCGTCCATCCATCCTCCTATCCTCTCATTTATCCATCCATCCACCTATCAATCCATCCATCCATCCACCTATCAATCCATCCATCCATCCTCCTATCCTCTCATTTATCCATCGATCCACCTA",
                "TTTTTTTTTTTTTTTTTTTTT",Collections.singletonList("asm009963:tig00000"),
                Arrays.asList(new SimpleInterval("chr6:857170-857564"), new SimpleInterval("chr6:857564-857767"), new SimpleInterval("chr6:857767-857852")), Arrays.asList("2"));
        relevantAttributes = new RelevantAttributes(complex);
        data.add(new Object[]{makeDeletion(new SimpleInterval("chr6:857169-857562"), Allele.create("G", true)).make(), relevantAttributes, true});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestIsConsistentWithCPX")
    public void testIsConsistentWithCPX(final VariantContext simpleVariant, final RelevantAttributes complexVarAttributes,
                                        final boolean expected) {
        Assert.assertEquals(isConsistentWithCPX(simpleVariant, complexVarAttributes),
                            expected);
    }

    @DataProvider(name = "forTestDeletionConsistencyCheck")
    private Object[][] forTestDeletionConsistencyCheck() {
        final List<Object[]> data = new ArrayList<>(20);

        VariantContext del = makeDeletion(new SimpleInterval("chr6:857169-857562"), Allele.create("G", true)).make();
        HashSet<SimpleInterval> missingSegments = Sets.newHashSet(new SimpleInterval("chr6:857170-857564"), new SimpleInterval("chr6:857767-857852"));

        data.add(new Object[]{del, Collections.emptySet(), false});
        data.add(new Object[]{del, missingSegments, true});

        data.add(new Object[]{del, Collections.singleton(new SimpleInterval("chr6:857767-857852")), false});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestDeletionConsistencyCheck")
    public void testDeletionConsistencyCheck(final VariantContext simple, final Set<SimpleInterval> missingSegments,
                                             final boolean expected) {
        Assert.assertEquals(deletionConsistencyCheck(simple, missingSegments),
                            expected);
    }

    @DataProvider(name = "forPostProcessConvertShortDupToIns")
    private Object[][] forPostProcessConvertShortDupToIns() {
        final List<Object[]> data = new ArrayList<>(20);

        // inversion, no effect
        VariantContext var = makeInversion(new SimpleInterval("chr1:10001-10100"), Allele.create("A", true)).make();
        data.add(new Object[]{var, var});

        // deletion, no effect
        var = makeDeletion(new SimpleInterval("chr1:10001-10100"), Allele.create("A", true)).make();
        data.add(new Object[]{var, var});

        // insertion, no effect
        var = makeInsertion("chr1", 11111, 11111, 100, Allele.create("A", true)).make();
        data.add(new Object[]{var, var});

        // dup large enough
        var = new VariantContextBuilder().chr("chr2").start(241987322).stop(241987322)
                .id("INS-DUPLICATION-TANDEM-EXPANSION_chr2_241987322_241987322_CPX_DERIVED")
                .alleles(Arrays.asList(Allele.create("T", true), Allele.create("<DUP>")))
                .attribute(VCFConstants.END_KEY, 241987322)
                .attribute(CONTIG_NAMES, "asm004634:tig00000")
                .attribute(ALIGN_LENGTHS, 1125)
                .attribute(HQ_MAPPINGS, 1)
                .attribute(MAPPING_QUALITIES, 60)
                .attribute(TOTAL_MAPPINGS,1)
                .attribute(MAX_ALIGN_LENGTH, 1125)
                .attribute(INSERTED_SEQUENCE, "AGTGGATGGCGTTGGGTGTACTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCAGAGATCCAGCTGATGGCATTCAGCGTACTCGGAGATCCAGTTGATGGTGTTGGGTGTTCTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCAGAGATCCAGTTGATGGCATTCAGTGTACTCGGAGATCTAGTCGATGGC")
                .attribute(INSERTED_SEQUENCE_LENGTH, 196)
                .attribute(INSERTED_SEQUENCE_MAPPINGS, "1102_1242_chr2:241987361-241987501_-_1101H141M1284H_60_8_101_O,1176_1373_chr2:241987323-241987520_-_1175H198M1153H_60_12_138_O")
                .attribute(SEQ_ALT_HAPLOTYPE, "GTTGGGTGTACTCGGAGATCCGGTCGATGGCGTTGGGTGTACTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCGGAGATCCAGTGGATGGCGTTGGGTGTACTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCGGAGATCCAGTGGATGGCGTTGGGTGTACTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCAGAGATCCAGCTGATGGCATTCAGCGTACTCGGAGATCCAGTTGATGGTGTTGGGTGTTCTCGGAGATCCAGTCGGTGGCGTTGGGTGTACTCAGAGATCCAGTTGATGGCATTCAGTGTACTCGGAGATCTAGTCGATGGCGTTGGGTGTACTCGGAGATCCAGTTGATGGCATTCAGCGTACTCGGAGATCCAGTTGATGGTGTTGGGTGTTCTCGGAGATCCAGTCGGTGGTGTTGGGTGTACTCAGAGATCCAGTTGATGGCATTCATTGTACTCGGAGATCCAGTCGATGGCGTTGGGTGTACTTGGAGATCCAGTCGGTGGCGTTGGGTGTACTTGGAGATCC")
                .attribute(SVLEN, 403)
                .attribute(SVTYPE, "DUP")
                .attribute(CPX_EVENT_KEY, "CPX_chr2:241987323-241987529")
                .attribute(DUPLICATION_NUMBERS, "1,2")
                .attribute(DUP_ORIENTATIONS, "++")
                .attribute(DUP_REPEAT_UNIT_REF_SPAN, "chr2:241987323-241987529")
                .attribute(DUP_SEQ_CIGARS,"207M,207M")
                .attribute(DUP_TAN_EXPANSION_STRING, "")
                .make();
        data.add(new Object[]{var, var});

        // dup not large enough
        var = new VariantContextBuilder().chr("chr2").start(83340906).stop(83340906)
                .id("INS-DUPLICATION-TANDEM-EXPANSION_chr2_83340906_83340906_CPX_DERIVED")
                .alleles(Arrays.asList(Allele.create("T", true), Allele.create("<DUP>")))
                .attribute(VCFConstants.END_KEY, 83340906)
                .attribute(CONTIG_NAMES, "asm003204:tig00000")
                .attribute(ALIGN_LENGTHS, 58)
                .attribute(HQ_MAPPINGS, 1)
                .attribute(MAPPING_QUALITIES, 60)
                .attribute(TOTAL_MAPPINGS,1)
                .attribute(MAX_ALIGN_LENGTH, 58)
                .attribute(INSERTED_SEQUENCE, "TTAATATTATATTATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAA")
                .attribute(INSERTED_SEQUENCE_LENGTH, 60)
                .attribute(INSERTED_SEQUENCE_MAPPINGS, "1280_1328_chr2:83340902-83340950_+_1279H49M62H_60_0_49_O")
                .attribute(SEQ_ALT_HAPLOTYPE, "TATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAATATATTATAATATATTTTAATATTATATTATATTATAATATATT")
                .attribute(SVLEN, 104)
                .attribute(SVTYPE, "DUP")
                .attribute(CPX_EVENT_KEY, "CPX_chr2:83340902-83340950")
                .attribute(DUPLICATION_NUMBERS, "1,2")
                .attribute(DUP_ORIENTATIONS, "+-")
                .attribute(DUP_REPEAT_UNIT_REF_SPAN, "chr2:83340907-83340950")
                .attribute(DUP_SEQ_CIGARS,"44M,44M")
                .attribute(DUP_TAN_EXPANSION_STRING, "")
                .make();
        final VariantContext insertion = new VariantContextBuilder().chr("chr2").start(83340906).stop(83340906)
                .id("INS-DUPLICATION-TANDEM-EXPANSION_chr2_83340906_83340906_CPX_DERIVED")
                .alleles(Arrays.asList(Allele.create("T", true), Allele.create("<INS>")))
                .attribute(VCFConstants.END_KEY, 83340906)
                .attribute(CONTIG_NAMES, "asm003204:tig00000")
                .attribute(ALIGN_LENGTHS, 58)
                .attribute(HQ_MAPPINGS, 1)
                .attribute(MAPPING_QUALITIES, 60)
                .attribute(TOTAL_MAPPINGS,1)
                .attribute(MAX_ALIGN_LENGTH, 58)
                .attribute(INSERTED_SEQUENCE, "TTAATATTATATTATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAA")
                .attribute(INSERTED_SEQUENCE_LENGTH, 60)
                .attribute(INSERTED_SEQUENCE_MAPPINGS, "1280_1328_chr2:83340902-83340950_+_1279H49M62H_60_0_49_O")
                .attribute(SEQ_ALT_HAPLOTYPE, "TATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAATATTATATTATATTATAATATATTTTAATATATTATAATATATTTTAATATTATATTATATTATAATATATT")
                .attribute(SVLEN, 104)
                .attribute(SVTYPE, "INS")
                .attribute(CPX_EVENT_KEY, "CPX_chr2:83340902-83340950")
                .attribute(DUPLICATION_NUMBERS, "1,2")
                .attribute(DUP_ORIENTATIONS, "+-")
                .attribute(DUP_REPEAT_UNIT_REF_SPAN, "chr2:83340907-83340950")
                .attribute(DUP_SEQ_CIGARS,"44M,44M")
                .attribute(DUP_TAN_EXPANSION_STRING, "")
                .make();
        data.add(new Object[]{var, insertion});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forPostProcessConvertShortDupToIns")
    public void testPostProcessConvertShortDupToIns(final VariantContext simpleVariant,
                                                    final VariantContext expected) {
        VariantContextTestUtils.assertVariantContextsAreEqual(postProcessConvertShortDupToIns(simpleVariant), expected, Collections.emptyList());
    }

    @DataProvider(name = "forPostProcessConvertReplacementToFatInsOrInsAndDel")
    private Object[][] forPostProcessConvertReplacementToFatInsOrInsAndDel() {
        final List<Object[]> data = new ArrayList<>(20);

        // rare but possible: no variant would be emitted
        data.add(new Object[]{new VariantContextBuilder().chr("chr9").start(108455205).stop(108455252).alleles(Arrays.asList(Allele.create("T", true), Allele.create("<DEL>"))).attribute(VCFConstants.END_KEY, 108455252).attribute(SVTYPE, "DEL").attribute(SVLEN, -47).attribute(INSERTED_SEQUENCE_LENGTH, 7).attribute(INSERTED_SEQUENCE, "ATCTTAT").make(),
                              Collections.emptyList()
        });

        // fat insertion
        VariantContext deletion = makeDeletion(new SimpleInterval("chr21:23428920-23428967"), Allele.create("T", true)).attribute(ALIGN_LENGTHS, "56,56").attribute(CPX_EVENT_KEY, "CPX_chr21:23428920-23429023").attribute(CONTIG_NAMES, "asm029052:tig00000,asm029052:tig00001").attribute(HQ_MAPPINGS, 2).attribute(INSERTED_SEQUENCE, "ATAAATATATATAATAATATATAATTATATATTATATTATATAATATAATAATATATATTATAATACATTATATAATATATTATA").attribute(INSERTED_SEQUENCE_LENGTH, 85).attribute(INSERTED_SEQUENCE_MAPPINGS, "1330_1385_chr21:23428968-23429023_+_1329H56M1200H_49_4_36_O,1330_1385_chr21:23428968-23429023_+_1329H56M1200H_49_4_36_O").attribute(MAPPING_QUALITIES, "60,60").attribute(MAX_ALIGN_LENGTH, 56).attribute(SEQ_ALT_HAPLOTYPE, "ATAAATATATATAATAATATATAATTATATATTATATTATATAATATAATAATATATATTATAATACATTATATAATATATTATA").attribute(TOTAL_MAPPINGS, 2).make();
        VariantContext fatInsertion = makeInsertion("chr21", 23428920, 23428967, 85, Allele.create("TTTATATAAATATATATAAATATATAATATATAATAATATAATATAAT", true)).attribute(ALIGN_LENGTHS, "56,56").attribute(CPX_EVENT_KEY, "CPX_chr21:23428920-23429023").attribute(CONTIG_NAMES, "asm029052:tig00000,asm029052:tig00001").attribute(HQ_MAPPINGS, 2).attribute(INSERTED_SEQUENCE, "ATAAATATATATAATAATATATAATTATATATTATATTATATAATATAATAATATATATTATAATACATTATATAATATATTATA").attribute(INSERTED_SEQUENCE_LENGTH, 85).attribute(MAPPING_QUALITIES, "60,60").attribute(MAX_ALIGN_LENGTH, 56).attribute(SEQ_ALT_HAPLOTYPE, "ATAAATATATATAATAATATATAATTATATATTATATTATATAATATAATAATATATATTATAATACATTATATAATATATTATA").attribute(TOTAL_MAPPINGS, 2).make();
        data.add(new Object[]{deletion,
                              Collections.singletonList(fatInsertion)
        });

        // deletion with small insertion, i.e. no modification
        VariantContext deletionWithMicroInsertion = new VariantContextBuilder().chr("chr20").start(63093346).stop(63094245).alleles(Arrays.asList(Allele.create("G", true), Allele.create("<DEL>"))).attribute(VCFConstants.END_KEY, 63094245).attribute(SVTYPE, "DEL").attribute(SVLEN, -899).attribute(INSERTED_SEQUENCE_LENGTH, 1).attribute(INSERTED_SEQUENCE, "T").attribute(ALIGN_LENGTHS,942).attribute(CPX_EVENT_KEY, "CPX_chr20:63092255-63094246").attribute(CONTIG_NAMES, "asm028762:tig00002").attribute(HQ_MAPPINGS, 1).attribute(MAPPING_QUALITIES, "60").attribute(MAX_ALIGN_LENGTH,942).attribute(SEQ_ALT_HAPLOTYPE, "T").attribute(TOTAL_MAPPINGS, 1).make();
        data.add(new Object[]{deletionWithMicroInsertion,
                              Collections.singletonList(deletionWithMicroInsertion)
        });

        // deletion and insertion at the same time (location massaging)
        VariantContext sourceDeletion = makeDeletion(new SimpleInterval("chr20:440444-440697"), Allele.create("A", true)).attribute(CPX_EVENT_KEY, "CPX_chrX:439692-440698").attribute(CONTIG_NAMES, "asm030101:tig00001").attribute(TOTAL_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(HQ_MAPPINGS, 1).attribute(ALIGN_LENGTHS, 170).attribute(MAX_ALIGN_LENGTH, 170).attribute(INSERTED_SEQUENCE_LENGTH, 60).attribute(INSERTED_SEQUENCE, "TTCATACACACACAGATACACACCCGCGCACACACAGATGCACACACACACCCGTACACT").attribute(SEQ_ALT_HAPLOTYPE, "TTCATACACACACAGATACACACCCGCGCACACACAGATGCACACACACACCCGTACACT").make();
        VariantContext linkedDel = makeDeletion(new SimpleInterval("chr20:440444-440697"), Allele.create("A", true)).attribute(LINK, "INS_chr20_440444_440444_CPX_DERIVED").attribute(CPX_EVENT_KEY, "CPX_chrX:439692-440698").attribute(CONTIG_NAMES, "asm030101:tig00001").attribute(TOTAL_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(HQ_MAPPINGS, 1).attribute(ALIGN_LENGTHS, 170).attribute(MAX_ALIGN_LENGTH, 170).make();
        VariantContext linkedIns = makeInsertion("chr20", 440444, 440444, 60, Allele.create("A", true)).attribute(LINK, "DEL_chr20_440444_440697_CPX_DERIVED").attribute(CPX_EVENT_KEY, "CPX_chrX:439692-440698").attribute(CONTIG_NAMES, "asm030101:tig00001").attribute(TOTAL_MAPPINGS, 1).attribute(MAPPING_QUALITIES, 60).attribute(HQ_MAPPINGS, 1).attribute(ALIGN_LENGTHS, 170).attribute(MAX_ALIGN_LENGTH, 170).attribute(INSERTED_SEQUENCE_LENGTH, 60).attribute(INSERTED_SEQUENCE_LENGTH, 60).attribute(INSERTED_SEQUENCE, "TTCATACACACACAGATACACACCCGCGCACACACAGATGCACACACACACCCGTACACT").attribute(SEQ_ALT_HAPLOTYPE, "TTCATACACACACAGATACACACCCGCGCACACACAGATGCACACACACACCCGTACACT").make();
        data.add(new Object[]{sourceDeletion,
                              Arrays.asList(linkedIns, linkedDel)
        });

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forPostProcessConvertReplacementToFatInsOrInsAndDel")
    public void testPostProcessConvertReplacementToFatInsOrInsAndDel(final VariantContext simpleVariant,
                                                                     final List<VariantContext> expected) {
        assertVariantsAreEqual(postProcessConvertReplacementToFatInsOrInsAndDel(simpleVariant, b38_reference_chr20_chr21).collect(Collectors.toList()),
                               expected, Collections.emptyList(), b38_seqDict_chr20_chr21);
    }

    @DataProvider(name = "forTestRemoveDuplicates")
    private Object[][] forTestRemoveDuplicates() {
        final List<Object[]> data = new ArrayList<>(20);

        final List<VariantContext> sourceWithLessAnnotations = new ArrayList<>();
        final List<VariantContext> sourceWithMoreAnnotations = new ArrayList<>();
        final List<VariantContext> expected = new ArrayList<>();

        final VariantContext firstInsertion = makeInsertion("chr21", 46069065, 46069065, 60, Allele.create("C", true)).attribute(CPX_EVENT_KEY, "CPX_chr21:46069065-46069209").attribute(CONTIG_NAMES, "asm029362:tig00001,asm029362:tig00002").make();
        final VariantContext firstInsertionWithMoreAnnotations = makeInsertion("chr21", 46069065, 46069065, 60, Allele.create("C", true)).attribute(CPX_EVENT_KEY, "CPX_chr21:46069065-46069209").attribute(CONTIG_NAMES, "asm029362:tig00001,asm029362:tig00002")
                .attribute(TOTAL_MAPPINGS, 2).attribute(MAPPING_QUALITIES, "60,60").attribute(HQ_MAPPINGS, "2").attribute(ALIGN_LENGTHS, "91,91").attribute(MAX_ALIGN_LENGTH, "91").attribute(INSERTED_SEQUENCE, "CTAGGTGTGTGCATGTGTGCACACGTGTGTGCATGTGTGTGCATGTGTGCACACGTGTGT").attribute(INSERTED_SEQUENCE_LENGTH, 60).attribute(SEQ_ALT_HAPLOTYPE, "CTAGGTGTGTGCATGTGTGCACACGTGTGTGCATGTGTGTGCATGTGTGCACACGTGTGT").make();
        sourceWithLessAnnotations.add(firstInsertion);
        sourceWithMoreAnnotations.add(firstInsertionWithMoreAnnotations);
        expected.add(firstInsertionWithMoreAnnotations);

        final VariantContext firstDeletion = makeDeletion(new SimpleInterval("chr21:46069156-46069208"), Allele.create("A", true)).attribute(CPX_EVENT_KEY, "CPX_chr21:46069065-46069209").attribute(CONTIG_NAMES, "asm029362:tig00001,asm029362:tig00002").make();
        final VariantContext firstDeletionWitMoreAnnotations = makeDeletion(new SimpleInterval("chr21:46069156-46069208"), Allele.create("A", true)).attribute(CPX_EVENT_KEY, "CPX_chr21:46069065-46069209").attribute(CONTIG_NAMES, "asm029362:tig00001,asm029362:tig00002")
                .attribute(TOTAL_MAPPINGS, 2).attribute(MAPPING_QUALITIES, "60,60").attribute(HQ_MAPPINGS, "2").attribute(ALIGN_LENGTHS, "91,91").attribute(MAX_ALIGN_LENGTH, "91").make();
        sourceWithLessAnnotations.add(firstDeletion);
        sourceWithMoreAnnotations.add(firstDeletionWitMoreAnnotations);
        expected.add(firstDeletionWitMoreAnnotations);

        // locations below seems inconsistent with the annotations, but that's because we artificially put data on chr20 and 21
        final VariantContext insertionFromStringParsing = makeInsertion("chr20", 439692, 439692, 130, Allele.create("A", true)).attribute(CPX_EVENT_KEY, "CPX_chrX:439692-440698").attribute(CONTIG_NAMES, "asm030101:tig00001").make();
        final VariantContext deletionFromStringParsing = makeDeletion(new SimpleInterval("chr20:439692-440161"), Allele.create("A", true)).attribute(CPX_EVENT_KEY, "CPX_chrX:439692-440698").attribute(CONTIG_NAMES, "asm030101:tig00001").make();
        sourceWithLessAnnotations.add(insertionFromStringParsing);
        expected.add(insertionFromStringParsing);
        sourceWithLessAnnotations.add(deletionFromStringParsing);
        expected.add(deletionFromStringParsing);

        final VariantContext inversionFromStringParsing = makeInversion(new SimpleInterval("chr21:187497346-187497595"), Allele.create("A", true)).attribute(CPX_EVENT_KEY, "CPX_chr1:187495696-187497598").attribute(CONTIG_NAMES, "asm001762:tig00000,asm001762:tig00001,asm001763:tig00000").make();
        sourceWithLessAnnotations.add(inversionFromStringParsing);
        expected.add(inversionFromStringParsing);

        final VariantContext tandemDuplicationFromPairIteration = new VariantContextBuilder()
                .chr("chr20").start(56839685).stop(56839685)
                .id("INS-DUPLICATION-TANDEM-EXPANSION_chrY_56839685_56839685_CPX_DERIVED")
                .alleles(Arrays.asList(Allele.create("A", true), Allele.create("<DUP>")))
                .attribute(VCFConstants.END_KEY, 56839685)
                .attribute(TOTAL_MAPPINGS, 1)
                .attribute(MAPPING_QUALITIES, "60")
                .attribute(HQ_MAPPINGS, 1)
                .attribute(ALIGN_LENGTHS, 770)
                .attribute(MAX_ALIGN_LENGTH, 770)
                .attribute(INSERTED_SEQUENCE, "CTGTTGACTAGTCTTTGCCTACAGAGGGCGTTGTGACATATCTCTGCACTGATCTCTCAGGTGAGGTAACTTCTCTAGTCTCTGCCTACAGAGGG")
                .attribute(INSERTED_SEQUENCE_LENGTH, 95)
                .attribute(DUP_REPEAT_UNIT_REF_SPAN, "chrY:56839686-56839794")
                .attribute(DUPLICATION_NUMBERS, "1,2")
                .attribute(DUP_ORIENTATIONS, "++")
                .attribute(DUP_SEQ_CIGARS, "109M,109M")
                .attribute(DUP_TAN_EXPANSION_STRING, "")
                .attribute(SEQ_ALT_HAPLOTYPE, "CATTGTGACTTATCTCTGCACTGATCACCCAGGTGATGTAACTCTTGTCTAGGCTCTGGCCACAGGGACATAGTGACATATATCTGCACTGATCACACAGGTAATGTAACTGTTGACTAGTCTTTGCCTACAGAGGGCGTTGTGACATATCTCTGCACTGATCTCTCAGGTGAGGTAACTTCTCTAGTCTCTGCCTACAGAGGGCATTGTGACATCACTCTGCAATGATCACCCAGGTGATGTAACTCTTGTCTAGGCTCTGCCTACATGGACATTGTGACATGTCTCTGCACTGATCACCCAGGTGATGTAA")
                .attribute(CPX_EVENT_KEY, "CPX_chrY:56838784-56839794")
                .attribute(CONTIG_NAMES, "asm031346:tig00004")
                .make();
        sourceWithMoreAnnotations.add(tandemDuplicationFromPairIteration);
        expected.add(tandemDuplicationFromPairIteration);

        data.add(new Object[]{sourceWithLessAnnotations, sourceWithMoreAnnotations, expected});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestRemoveDuplicates")
    public void testRemoveDuplicates(final List<VariantContext> sourceWithLessAnnotations, final List<VariantContext> sourceWithMoreAnnotations,
                                     final List<VariantContext> expected) {
        assertVariantsAreEqual(removeDuplicates(sourceWithLessAnnotations, sourceWithMoreAnnotations), expected,
                Collections.emptyList(), bareBoneHg38SAMSeqDict);
    }

    //==================================================================================================================

    @DataProvider(name = "forTestGetInsFromOneEnd")
    private Object[][] forTestGetInsFromOneEnd() {
        final List<Object[]> data = new ArrayList<>(20);

        final SimpleInterval dummyInsertionPos = new SimpleInterval("chr1", 100, 100);
        final Allele dummyRefAllele = Allele.create("A", true);
        final List<Integer> refSegmentLengths = Arrays.asList(30, 40, 50, 10);

        data.add(new Object[]{true, 4, dummyInsertionPos, dummyRefAllele, refSegmentLengths, Arrays.asList("1","2","3","4","1","2","3","4"), true,
                              makeInsertion(dummyInsertionPos.getContig(), dummyInsertionPos.getStart(), dummyInsertionPos.getEnd(), 131, dummyRefAllele).make()});
        data.add(new Object[]{false, 7, dummyInsertionPos, dummyRefAllele, refSegmentLengths, Arrays.asList("1","2","3","4","1","2","3","4"), true,
                              null});

        data.add(new Object[]{true, 2, dummyInsertionPos, dummyRefAllele, refSegmentLengths, Arrays.asList(UNMAPPED_INSERTION + "-50", "-chrX:10001-10009", "1","2","3","3","1","4"), true,
                              makeInsertion(dummyInsertionPos.getContig(), dummyInsertionPos.getStart(), dummyInsertionPos.getEnd(), 60, dummyRefAllele).make()});
        data.add(new Object[]{false, 5, dummyInsertionPos, dummyRefAllele, refSegmentLengths, Arrays.asList("1","2","3","3","1","4", UNMAPPED_INSERTION + "-49"), true,
                              makeInsertion(dummyInsertionPos.getContig(), dummyInsertionPos.getStart(), dummyInsertionPos.getEnd(), 50, dummyRefAllele).make()});
        data.add(new Object[]{false, 5, dummyInsertionPos, dummyRefAllele, refSegmentLengths, Arrays.asList("1","2","3","3","1","4", UNMAPPED_INSERTION + "-49"), false,
                              null});

        data.add(new Object[]{true, 1, dummyInsertionPos, dummyRefAllele, Collections.singletonList(20), Arrays.asList(UNMAPPED_INSERTION + "-100", "1"), true,
                              makeInsertion(dummyInsertionPos.getContig(), dummyInsertionPos.getStart(), dummyInsertionPos.getEnd(), 101, dummyRefAllele).make()});
        data.add(new Object[]{false, 0, dummyInsertionPos, dummyRefAllele, Collections.singletonList(20), Arrays.asList("1", UNMAPPED_INSERTION + "-100"), true,
                              makeInsertion(dummyInsertionPos.getContig(), dummyInsertionPos.getStart(), dummyInsertionPos.getEnd(), 101, dummyRefAllele).make()});
        data.add(new Object[]{true, 2, dummyInsertionPos, dummyRefAllele, Collections.singletonList(20), Arrays.asList("1", "-chrX:10001-10039", "1"), true,
                              makeInsertion(dummyInsertionPos.getContig(), dummyInsertionPos.getStart(), dummyInsertionPos.getEnd(), 60, dummyRefAllele).make()});
        data.add(new Object[]{true, 2, dummyInsertionPos, dummyRefAllele, Collections.singletonList(20), Arrays.asList("1", "-chrX:10001-10009", "1"), true,
                              null});
        data.add(new Object[]{false, 2, dummyInsertionPos, dummyRefAllele, Collections.singletonList(20), Arrays.asList("1", "-chrX:10001-10039", "1"), true,
                              null});
        data.add(new Object[]{false, 2, dummyInsertionPos, dummyRefAllele, Collections.singletonList(20), Arrays.asList("1", "-chrX:10001-10009", "1"), true,
                              null});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestGetInsFromOneEnd")
    public void testGetInsFromOneEnd(final boolean fromFront, final int idxFirstMatch,
                                     final SimpleInterval insertionStartAndStop, final Allele anchorBaseRefAllele,
                                     final List<Integer> refSegmentLengths, final List<String> altArrangement,
                                     final boolean shouldIncreaseInsLenByOne,
                                     final VariantContext expected) {
        if (expected != null) {
            VariantContextTestUtils.assertVariantContextsAreEqual(getInsFromOneEnd(fromFront, idxFirstMatch, insertionStartAndStop, anchorBaseRefAllele, refSegmentLengths, altArrangement, shouldIncreaseInsLenByOne).make(),
                                                                    expected, Collections.emptyList());
        } else {
            Assert.assertNull(getInsFromOneEnd(fromFront, idxFirstMatch, insertionStartAndStop, anchorBaseRefAllele, refSegmentLengths, altArrangement, shouldIncreaseInsLenByOne));
        }
    }

    @DataProvider(name = "forTestGetInsLen")
    private Object[][] forTestGetInsLen() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{UNMAPPED_INSERTION + "-12", Collections.emptyList(), 12});

        data.add(new Object[]{UNMAPPED_INSERTION + "-12", Collections.singletonList(89), 12});

        data.add(new Object[]{"3", Arrays.asList(30, 40, 50, 10), 50});

        data.add(new Object[]{"-4", Arrays.asList(30, 40, 50, 10), 10});

        data.add(new Object[]{"chr1:1000-1000", Arrays.asList(30, 40, 50, 10), 1});

        data.add(new Object[]{"chr1:1000-1001", Arrays.asList(30, 40, 50, 10), 2});

        data.add(new Object[]{"-chr1:1001-1050", Arrays.asList(30, 40, 50, 10), 50});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestGetInsLen")
    public void testGetInsLen(final String description, final List<Integer> refSegmentLengths,
                              final int expected) {
        Assert.assertEquals(getInsLen(description, refSegmentLengths),
                            expected);
    }

    @DataProvider(name = "forGetMissingAndPresentAndInvertedSegments")
    private Object[][] forGetMissingAndPresentAndInvertedSegments() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{Arrays.asList(new SimpleInterval("chr1:1001-1100")), Arrays.asList("-1"),
                              new Tuple3<>(Collections.emptySet(), Collections.emptySet(), Arrays.asList(1))});

        final List<SimpleInterval> refSegments = Arrays.asList(new SimpleInterval("chr1:1001-1100"), new SimpleInterval("chr1:1100-1200"), new SimpleInterval("chr1:1200-1300"), new SimpleInterval("chr1:1300-1400"));
        data.add(new Object[]{refSegments, Arrays.asList("1","3","4"),
                              new Tuple3<>(Sets.newHashSet(new SimpleInterval("chr1:1100-1200")), new TreeSet<>(Sets.newHashSet(1,3,4)), Collections.emptyList())});

        data.add(new Object[]{refSegments, Arrays.asList("1","4","2","3"),
                              new Tuple3<>(Collections.emptySet(), new TreeSet<>(Sets.newHashSet(1,2,3,4)), Collections.emptyList())});

        data.add(new Object[]{refSegments, Arrays.asList("-4","1","3","2"),
                              new Tuple3<>(Collections.emptySet(), new TreeSet<>(Sets.newHashSet(1,2,3)), Arrays.asList(4))});

        data.add(new Object[]{refSegments, Arrays.asList("1","3","4","-1"),
                              new Tuple3<>(Sets.newHashSet(new SimpleInterval("chr1:1100-1200")), new TreeSet<>(Sets.newHashSet(1,3,4)), Arrays.asList(1))});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forGetMissingAndPresentAndInvertedSegments")
    public void testGetMissingAndPresentAndInvertedSegments(final List<SimpleInterval> refSegments,
                                                            final List<String> altArrangements,
                                                            Tuple3<Set<SimpleInterval>, Set<Integer>, List<Integer>> expected) {

        final Tuple3<Set<SimpleInterval>, Set<Integer>, List<Integer>> actual =
                getMissingAndPresentAndInvertedSegments(refSegments, altArrangements);
        Assert.assertEquals(actual._1(), expected._1());
        Assert.assertEquals(actual._2(), expected._2());
        Assert.assertEquals(actual._3(), expected._3());
    }

    //==================================================================================================================

    private VariantContext makeTestComplexVariant(final SimpleInterval affectedRefRegion, final int svLen,
                                                  final String referenceBases, final String altSeqBases,
                                                  final List<String> contigNames,
                                                  final List<SimpleInterval> referenceSegments, final List<String> altArrangement) {
        final int maxAlignmentLength = random.nextInt(4000) + 1;
        int evidenceContigCount = contigNames.size();
        StringBuilder mqs = new StringBuilder();
        for (int i = 0; i <= evidenceContigCount; ++i) mqs.append(random.nextInt(61)).append(",");
        String mq = mqs.toString();

        final VariantContextBuilder builder = new VariantContextBuilder()
                .chr(affectedRefRegion.getContig()).start(affectedRefRegion.getStart()).stop(affectedRefRegion.getEnd())
                .alleles(Arrays.asList(Allele.create(referenceBases, true),
                        Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR))))
                .id(CPX_SV_SYB_ALT_ALLELE_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + affectedRefRegion.toString())
                .attribute(VCFConstants.END_KEY, affectedRefRegion.getEnd())
                .attribute(SVLEN, svLen)
                .attribute(SVTYPE, CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(CPX_EVENT_ALT_ARRANGEMENTS, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, altArrangement))
                .attribute(CONTIG_NAMES, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, contigNames))
                .attribute(SEQ_ALT_HAPLOTYPE, altSeqBases)
                .attribute(MAX_ALIGN_LENGTH, maxAlignmentLength)
                .attribute(MAPPING_QUALITIES, mq.substring(0, mq.length() - 1)); // drop last coma
        if (referenceSegments.isEmpty())
            return builder.make();
        else
            return builder
                    .attribute(CPX_SV_REF_SEGMENTS, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR,
                            referenceSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())))
                    .make();
    }

    private static void assertVariantsAreEqual(final Iterable<VariantContext> actual, final Iterable<VariantContext> expected,
                                               final List<String> attributesToIgnore, final SAMSequenceDictionary refSeqDict) {

        final List<VariantContext> actualList = SVVCFWriter.sortVariantsByCoordinate(Utils.stream(actual).collect(Collectors.toList()), refSeqDict);
        final List<VariantContext> expectedList = SVVCFWriter.sortVariantsByCoordinate(Utils.stream(expected).collect(Collectors.toList()), refSeqDict);
        if (actualList.size() != expectedList.size()) {
            throw new AssertionError("Two sources of variants are not of the same size. expected size: " + expectedList.size() + "actual size: " + actualList.size());
        }
        for (int i = 0; i < actualList.size(); ++i) {
            VariantContextTestUtils.assertVariantContextsAreEqual(actualList.get(i), expectedList.get(i), attributesToIgnore);
        }
    }

}
