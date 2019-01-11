package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import java.util.Arrays;
import java.util.List;

/**
 * A simple holder class for MNP data providers for the PIK3CA gene.
 * These cases were pulled directly from Oncotator (MUC16ChangeTestdata.py).
 * Created by jonn on 9/21/17.
 */
public class DataProviderForPik3caTestData {
    public static List<Object[]> providePik3caMnpData() {
        return Arrays.asList(
                new Object[]{"PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.SILENT,   GencodeFuncotation.VariantType.SNP, "T", "A", "g.chr3:178916619T>A", "+", "c.6T>A",    "c.(4-6)ccT>ccA",       "p.P2P"},
                new Object[]{"PIK3CA", 3, 178916620, 178916620, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178916620C>T", "+", "c.7C>T",    "c.(7-9)Cca>Tca",       "p.P3S"},
                new Object[]{"PIK3CA", 3, 178916617, 178916617, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178916617C>T", "+", "c.4C>T",    "c.(4-6)Cct>Tct",       "p.P2S"},
                new Object[]{"PIK3CA", 3, 178919220, 178919220, GencodeFuncotation.VariantClassification.SILENT,   GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178919220C>T", "+", "c.705C>T",  "c.(703-705)tcC>tcT",   "p.S235S"},
                new Object[]{"PIK3CA", 3, 178921433, 178921433, GencodeFuncotation.VariantClassification.SILENT,   GencodeFuncotation.VariantType.SNP, "A", "T", "g.chr3:178921433A>T", "+", "c.915A>T",  "c.(913-915)ccA>ccT",   "p.P305P"},
                new Object[]{"PIK3CA", 3, 178922366, 178922366, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "T", "A", "g.chr3:178922366T>A", "+", "c.1135T>A", "c.(1135-1137)Tcc>Acc", "p.S379T"},
                new Object[]{"PIK3CA", 3, 178928317, 178928317, GencodeFuncotation.VariantClassification.SILENT,   GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178928317C>T", "+", "c.1503C>T", "c.(1501-1503)tcC>tcT", "p.S501S"},
                new Object[]{"PIK3CA", 3, 178936091, 178936091, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "G", "A", "g.chr3:178936091G>A", "+", "c.1633G>A", "c.(1633-1635)Gag>Aag", "p.E545K"},
                new Object[]{"PIK3CA", 3, 178937063, 178937063, GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178937063C>T", "+", "c.1744C>T", "c.(1744-1746)Cag>Tag", "p.Q582*"},
                new Object[]{"PIK3CA", 3, 178941890, 178941890, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "G", "A", "g.chr3:178941890G>A", "+", "c.2209G>A", "c.(2209-2211)Gag>Aag", "p.E737K"},
                new Object[]{"PIK3CA", 3, 178942511, 178942511, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178942511C>T", "+", "c.2318C>T", "c.(2317-2319)tCc>tTc", "p.S773F"},
                new Object[]{"PIK3CA", 3, 178942523, 178942523, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "G", "A", "g.chr3:178942523G>A", "+", "c.2330G>A", "c.(2329-2331)aGg>aAg", "p.R777K"},
                new Object[]{"PIK3CA", 3, 178943785, 178943785, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "C", "T", "g.chr3:178943785C>T", "+", "c.2452C>T", "c.(2452-2454)Cgt>Tgt", "p.R818C"},
                new Object[]{"PIK3CA", 3, 178947158, 178947158, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "G", "A", "g.chr3:178947158G>A", "+", "c.2594G>A", "c.(2593-2595)gGc>gAc", "p.G865D"},
                new Object[]{"PIK3CA", 3, 178952085, 178952085, GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantType.SNP, "A", "T", "g.chr3:178952085A>T", "+", "c.3140A>T", "c.(3139-3141)cAt>cTt", "p.H1047L"}
        );
    }

    //TODO: Get indel data from a negative strand gene.

    /**
     * @return Test data for PIK3CA INDELs as taken from Oncotator:VariantClassifierTest.py:test_pik3ca_change_transcript:203-222
     */
    public static List<Object[]> providePik3caInDelData() {
        return Arrays.asList(

                // TODO:
                // 1 - must make sure that for insertions you keep the codons that are spanned by the
                // 2 - shift insertions 1 codon (3 bases) right to fix the leading base bug
                // 3 - properly align the alleles for capitalization and correct insertion location

                // For insertions, the inserted bases occur just AFTER the given start/end position (because start and end are the same).
                // For deletions, the deleted bases occur just AFTER the given start position.

                // Insertions:
                new Object[] { "PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "T",          "TCGA",   "g.chr3:178916619_178916620insCGA",      "+", "c.6_7insCGA",             "c.(7-9)cca>CGAcca",       "p.2_3insR" },
                new Object[] { "PIK3CA", 3, 178948159, 178948159, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "T",          "TGAG",   "g.chr3:178948159_178948160insGAG",      "+", "c.2931_2932insGAG",       "c.(2932-2934)gag>GAGgag", "p.977_978insE" },
                new Object[] { "PIK3CA", 3, 178948165, 178948165, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.INS, "G",          "GT",     "g.chr3:178948165_178948166insT",        "+", "c.e20+1G>GT",             "c.e20+1",                 null },
                new Object[] { "PIK3CA", 3, 178948166, 178948166, GencodeFuncotation.VariantClassification.INTRON,          GencodeFuncotation.VariantType.INS, "T",          "TT",     "g.chr3:178948166_178948167insT",        "+", "c.e20+2T>TT",             null,                      null },

                new Object[] { "PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",          "TCG",    "g.chr3:178916619_178916620insCG",       "+", "c.6_7insCG",              "c.(7-9)ccafs",            "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178948159, 178948159, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",          "TGA",    "g.chr3:178948159_178948160insGA",       "+", "c.2931_2932insGA",        "c.(2932-2934)gagfs",      "p.F980fs" },
                new Object[] { "PIK3CA", 3, 178948163, 178948163, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.INS, "A",          "AT",     "g.chr3:178948163_178948164insT",        "+", "c.2935_2936insT",         "c.(2935-2937)aggfs",      "p.R979fs" },
                new Object[] { "PIK3CA", 3, 178948154, 178948154, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "G",          "GGAATT", "g.chr3:178948154_178948155insGAATT",    "+", "c.2926_2927insGAATT",     "c.(2926-2928)gaafs",      "p.E976fs" },
                new Object[] { "PIK3CA", 3, 178948155, 178948155, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "A",          "AGAATT", "g.chr3:178948155_178948156insGAATT",    "+", "c.2927_2928insGAATT",     "c.(2926-2928)gaafs",      "p.F977fs" },
                new Object[] { "PIK3CA", 3, 178948156, 178948156, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "A",          "AGAATT", "g.chr3:178948156_178948157insGAATT",    "+", "c.2928_2929insGAATT",     "c.(2929-2931)tttfs",      "p.F977fs" },
                new Object[] { "PIK3CA", 3, 178948157, 178948157, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",          "TGAATT", "g.chr3:178948157_178948158insGAATT",    "+", "c.2929_2930insGAATT",     "c.(2929-2931)tttfs",      "p.F977fs" },

                // For this, had to modify protein change and codon change to be different from Oncotator.
                // Oncotator wanted to include 2 codons (TTT, GAG) and 2 Amino Acids (F977 and E978).
                // Because the first affected amino acid / codon is codon 977, I chose to only include that in the resulting
                // annotations to keep Funcotator consistent.
                new Object[] { "PIK3CA", 3, 178948158, 178948158, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",          "TGAATT", "g.chr3:178948158_178948159insGAATT",    "+", "c.2930_2931insGAATT",     "c.(2929-2931)tttfs",      "p.F977fs" },

                // Problemmatic case from Oncotator, but Funcotator provides reasonable output, so ground-truth is changed:
                // TODO: Check this case against issue #3482 - https://github.com/broadinstitute/gatk/issues/3842
//                new Object[] { "PIK3CA", 3, 178948163, 178948163, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.INS, "A",          "ATGA",   "g.chr3:178948163_178948164insTGA",      "+", "c.2935_2936insTGA",       "c.(2935-2937)agg>aTGAgg", "p.978_979insM" },
                new Object[] { "PIK3CA", 3, 178948163, 178948163, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.INS, "A",          "ATGA",   "g.chr3:178948163_178948164insTGA",      "+", "c.2935_2936insTGA",       "c.(2935-2937)agg>aTGAgg", "p.978_979insM" },

                // Deletions:
                new Object[] { "PIK3CA", 3, 178916937, 178916940, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "TGAA",       "T",     "g.chr3:178916938_178916940delGAA",       "+", "c.325_327delGAA",         "c.(325-327)gaadel",       "p.E109del" },
                new Object[] { "PIK3CA", 3, 178948159, 178948162, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "TGAG",       "T",     "g.chr3:178948160_178948162delGAG",       "+", "c.2932_2934delGAG",       "c.(2932-2934)gagdel",     "p.E978del" },
                new Object[] { "PIK3CA", 3, 178948159, 178948161, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "TGA",        "T",     "g.chr3:178948160_178948161delGA",        "+", "c.2932_2933delGA",        "c.(2932-2934)gagfs",      "p.R979fs" },
                new Object[] { "PIK3CA", 3, 178948153, 178948158, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "AGAATT",     "A",     "g.chr3:178948154_178948158delGAATT",     "+", "c.2926_2930delGAATT",     "c.(2926-2931)gaatttfs",   "p.E976fs" },
                new Object[] { "PIK3CA", 3, 178948153, 178948157, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "AGAAT",      "A",     "g.chr3:178948154_178948157delGAAT",      "+", "c.2926_2929delGAAT",      "c.(2926-2931)gaatttfs",   "p.E976fs" },
                new Object[] { "PIK3CA", 3, 178948159, 178948164, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.DEL, "TGAGAG",     "T",     "g.chr3:178948160_178948164delGAGAG",     "+", "c.2932_2936delGAGAG",     "c.(2932-2937)gagaggfs",   "p.E978fs" },
                new Object[] { "PIK3CA", 3, 178948165, 178948168, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.DEL, "GTGA",       "G",     "g.chr3:178948166_178948168delTGA",       "+", "c.e20+1GTGA>G",           "c.e20+2",                 null }

                // Known issue #3749 - https://github.com/broadinstitute/gatk/issues/3749.
                // When fixed these tests should be uncommented.
                //
                // Pathological cases - deletions run off the end of an exon.
                // Currently there's a bug that makes the reference only look at exonal bases, so Funcotator believes
                // that these are new alternate ref alleles and doesn't handle them correctly.
//                new Object[] { "PIK3CA", 3, 178948159, 178948167, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.DEL, "TGAGAGGTG",  "T",     "g.chr3:178948160_178948167delGAGAGGTG",  "+", "c.2932_2936delGAGAGGTG",  "c.(2932-2937)gagaggfs",   "p.ER978fs" },
//                new Object[] { "PIK3CA", 3, 178948159, 178948168, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.DEL, "TGAGAGGTGA", "T",     "g.chr3:178948160_178948168delGAGAGGTGA", "+", "c.2932_2936delGAGAGGTGA", "c.(2932-2937)gagagg>g",   "p.ER978del" }
        );
    }

    /**
     * @return Test data for PIK3CA INDELs as taken from Oncotator:VariantClassifierTest.py:test_reference_sequence_codon_construction_positive_strand:510-554
     */
    public static List<Object[]> providePik3caInDelData2() {

        //TODO: Uncomment failing tests here and fix the bugs they imply!

        return Arrays.asList(
                new Object[] { "PIK3CA", 3, 178916621, 178916621, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "C",       "CTAT",       "g.chr3:178916621_178916622insTAT",       "+", "c.8_9insTAT",        "c.(7-9)cca>ccTATa",  "p.3_4insI" },
                new Object[] { "PIK3CA", 3, 178916622, 178916622, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "A",       "ATAT",       "g.chr3:178916622_178916623insTAT",       "+", "c.9_10insTAT",       "c.(10-12)cga>TATcga",       "p.3_4insY" },
                new Object[] { "PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "T",       "TTAT",       "g.chr3:178916619_178916620insTAT",       "+", "c.6_7insTAT",        "c.(7-9)cca>TATcca",         "p.2_3insY" },
                new Object[] { "PIK3CA", 3, 178916620, 178916620, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "C",       "CTAT",       "g.chr3:178916620_178916621insTAT",       "+", "c.7_8insTAT",        "c.(7-9)cca>cTATca",         "p.P3LS" },
                new Object[] { "PIK3CA", 3, 178916622, 178916622, GencodeFuncotation.VariantClassification.IN_FRAME_INS,    GencodeFuncotation.VariantType.INS, "A",       "ACTTGAAGAA", "g.chr3:178916622_178916623insCTTGAAGAA", "+", "c.9_10insCTTGAAGAA", "c.(10-12)cga>CTTGAAGAAcga", "p.3_4insLEE" },

                new Object[] { "PIK3CA", 3, 178916621, 178916621, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "C",       "CTA",        "g.chr3:178916621_178916622insTA",        "+", "c.8_9insTA",         "c.(7-9)ccafs",          "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916622, 178916622, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "A",       "ATA",        "g.chr3:178916622_178916623insTA",        "+", "c.9_10insTA",        "c.(10-12)cgafs",            "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",       "TTA",        "g.chr3:178916619_178916620insTA",        "+", "c.6_7insTA",         "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916620, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "C",       "CTA",        "g.chr3:178916620_178916621insTA",        "+", "c.7_8insTA",         "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916621, 178916621, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "C",       "CT",         "g.chr3:178916621_178916622insT",         "+", "c.8_9insT",          "c.(7-9)ccafs",          "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916622, 178916622, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "A",       "AT",         "g.chr3:178916622_178916623insT",         "+", "c.9_10insT",         "c.(10-12)cgafs",            "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",       "TT",         "g.chr3:178916619_178916620insT",         "+", "c.6_7insT",          "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916620, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "C",       "CT",         "g.chr3:178916620_178916621insT",         "+", "c.7_8insT",          "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916621, 178916621, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "C",       "CTATT",      "g.chr3:178916621_178916622insTATT",      "+", "c.8_9insTATT",       "c.(7-9)ccafs",          "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916622, 178916622, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "A",       "ATATT",      "g.chr3:178916622_178916623insTATT",      "+", "c.9_10insTATT",      "c.(10-12)cgafs",            "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916619, 178916619, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "T",       "TTATT",      "g.chr3:178916619_178916620insTATT",      "+", "c.6_7insTATT",       "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916620, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "C",       "CTATT",      "g.chr3:178916620_178916621insTATT",      "+", "c.7_8insTATT",       "c.(7-9)ccafs",              "p.P3fs" },

                // Problemmatic case from Oncotator, but Funcotator provides reasonable output, so ground-truth is changed:
                // TODO: Check this case against issue #3482 - https://github.com/broadinstitute/gatk/issues/3842
//                new Object[] { "PIK3CA", 3, 178916618, 178916621, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "CTCC",     "C",         "g.chr3:178916619_178916621delTCC",       "+", "c.6_8delTCC",        "c.(4-9)cctcca>cca",         "p.2_3PP>P" },
                new Object[] { "PIK3CA", 3, 178916618, 178916621, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "CTCC",     "C",         "g.chr3:178916619_178916621delTCC",       "+", "c.6_8delTCC",        "c.(4-9)cctcca>cca",         "p.P3del" },


                new Object[] { "PIK3CA", 3, 178916619, 178916622, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "TCCA",     "T",         "g.chr3:178916620_178916622delCCA",       "+", "c.7_9delCCA",        "c.(7-9)ccadel",             "p.P3del" },
                new Object[] { "PIK3CA", 3, 178916620, 178916623, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "CCAC",     "C",         "g.chr3:178916621_178916623delCAC",       "+", "c.8_10delCAC",       "c.(7-12)ccacga>cga",        "p.P3del" },
                new Object[] { "PIK3CA", 3, 178916621, 178916624, GencodeFuncotation.VariantClassification.IN_FRAME_DEL,    GencodeFuncotation.VariantType.DEL, "CACG",     "C",         "g.chr3:178916622_178916624delACG",       "+", "c.9_11delACG",       "c.(7-12)ccacga>cca",        "p.R4del" },

                // Problemmatic case from Oncotator, but Funcotator provides reasonable output, so ground-truth is changed:
                // TODO: Check this case against issue #3482 - https://github.com/broadinstitute/gatk/issues/3842
//                new Object[] { "PIK3CA", 3, 178916618, 178916620, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CTC",      "C",         "g.chr3:178916619_178916620delTC",        "+", "c.6_7delTC",         "c.(4-9)cctccafs",           "p.PP2fs" },
                new Object[] { "PIK3CA", 3, 178916618, 178916620, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CTC",      "C",         "g.chr3:178916619_178916620delTC",        "+", "c.6_7delTC",         "c.(4-9)cctccafs",           "p.P3fs" },

                new Object[] { "PIK3CA", 3, 178916619, 178916621, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "TCC",      "T",         "g.chr3:178916620_178916621delCC",        "+", "c.7_8delCC",         "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916622, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CCA",      "C",         "g.chr3:178916621_178916622delCA",        "+", "c.8_9delCA",         "c.(7-9)ccafs",              "p.R4fs" },

                // Need to qualify these:
                new Object[] { "PIK3CA", 3, 178916621, 178916623, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CAC",      "C",         "g.chr3:178916622_178916623delAC",        "+", "c.9_10delAC",        "c.(7-12)ccacgafs",          "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916618, 178916619, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CT",       "C",         "g.chr3:178916619delT",                   "+", "c.6delT",            "c.(4-6)cctfs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916619, 178916620, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "TC",       "T",         "g.chr3:178916620delC",                   "+", "c.7delC",            "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916621, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CC",       "C",         "g.chr3:178916621delC",                   "+", "c.8delC",            "c.(7-9)ccafs",              "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916621, 178916622, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CA",       "C",         "g.chr3:178916622delA",                   "+", "c.9delA",            "c.(7-9)ccafs",              "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916618, 178916622, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CTCCA",    "C",         "g.chr3:178916619_178916622delTCCA",      "+", "c.6_9delTCCA",       "c.(4-9)cctccafs",           "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916619, 178916623, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "TCCAC",    "T",         "g.chr3:178916620_178916623delCCAC",      "+", "c.7_10delCCAC",      "c.(7-12)ccacgafs",          "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916624, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CCACG",    "C",         "g.chr3:178916621_178916624delCACG",      "+", "c.8_11delCACG",      "c.(7-12)ccacgafs",          "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916621, 178916625, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CACGA",    "C",         "g.chr3:178916622_178916625delACGA",      "+", "c.9_12delACGA",      "c.(7-12)ccacgafs",          "p.R4fs" },
                new Object[] { "PIK3CA", 3, 178916618, 178916625, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CTCCACGA", "C",         "g.chr3:178916619_178916625delTCCACGA",   "+", "c.6_12delTCCACGA",   "c.(4-12)cctccacgafs",       "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916619, 178916626, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "TCCACGAC", "T",         "g.chr3:178916620_178916626delCCACGAC",   "+", "c.7_13delCCACGAC",   "c.(7-15)ccacgaccafs",       "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916620, 178916627, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CCACGACC", "C",         "g.chr3:178916621_178916627delCACGACC",   "+", "c.8_14delCACGACC",   "c.(7-15)ccacgaccafs",       "p.P3fs" },
                new Object[] { "PIK3CA", 3, 178916621, 178916628, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "CACGACCA", "C",         "g.chr3:178916622_178916628delACGACCA",   "+", "c.9_15delACGACCA",   "c.(7-15)ccacgaccafs",       "p.R4fs" }
            );
    }
}
