package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import java.util.Arrays;
import java.util.List;

/**
 * A class to hold Indel data for MUC16 tests.
 * Created by jonn on 5/15/18.
 */
public class DataProviderForMuc16IndelData {

    static List<Object[]> provideIndelDataForMuc16() {
        return Arrays.asList(
                // ==============================
                // These exercise a test case for issue 4712 (https://github.com/broadinstitute/gatk/issues/4712):
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantType.INS, "T", "ATTCG", "g.chr19:9091890_9091891insATTCG", "-", null, null, null },
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantType.INS, "T", "CGCAT", "g.chr19:9091890_9091891insCGCA", "-", null, null, null },
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantType.INS, "T", "GCATC", "g.chr19:9091890_9091891insGCATC", "-", null, null, null },
                // ==============================

                // ==============================
                // Test cases for issue 5050:
                new Object[]{ "MUC16", 19, 8966646, 8966648, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.DEL, "TTA", "T", "g.chr19:8966647_8966648delTA", "-", "c.e81-4TAA>A", "c.e81+5", null },
                new Object[]{ "MUC16", 19, 8966818, 8966820, GencodeFuncotation.VariantClassification.INTRON,          GencodeFuncotation.VariantType.DEL, "TGG", "T", "g.chr19:8966819_8966820delGG", "-", "c.e81+2CCA>A", null, null },
                new Object[]{ "MUC16", 19, 8966813, 8966815, GencodeFuncotation.VariantClassification.SPLICE_SITE,     GencodeFuncotation.VariantType.DEL, "AGA", "A", "g.chr19:8966814_8966815delGA", "-", "c.43138_43139delTC", "c.(43138-43140)tctfs", "p.S14380fs" },
                new Object[]{ "MUC16", 19, 8966651, 8966653, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "ATC", "A", "g.chr19:8966652_8966653delTC", "-", "c.43300_43301delGA", "c.(43300-43302)gatfs", "p.D14434fs" },

                new Object[]{ "MUC16", 19, 8966647, 8966647, GencodeFuncotation.VariantClassification.INTRON,          GencodeFuncotation.VariantType.INS, "T", "TAG", "g.chr19:8966647_8966648insAG", "-", "c.e81-3A>CTA", null, null },
                new Object[]{ "MUC16", 19, 8966818, 8966818, GencodeFuncotation.VariantClassification.INTRON,          GencodeFuncotation.VariantType.INS, "T", "TAG", "g.chr19:8966818_8966819insAG", "-", "c.e81+2A>CTA", null, null },
                new Object[]{ "MUC16", 19, 8966814, 8966814, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "G", "GTC", "g.chr19:8966814_8966815insTC", "-", "c.43138_43139insGA", "c.(43138-43140)tctfs", "p.S14380fs" },
                new Object[]{ "MUC16", 19, 8966651, 8966651, GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantType.INS, "A", "ATG", "g.chr19:8966651_8966652insTG", "-", "c.43301_43302insCA", "c.(43300-43302)gatfs", "p.G14435fs" },

                new Object[]{ "MUC16", 19, 8966648, 8966648, GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantType.INS, "A", "AC", "g.chr19:8966648_8966649insC", "-", "c.e81-2T>GT", "c.e81+2", null },
                new Object[]{ "MUC16", 19, 8966817, 8966817, GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantType.INS, "C", "CT", "g.chr19:8966817_8966818insT", "-", "c.e81+1G>AG", "c.e81-1", null },
                new Object[]{ "MUC16", 19, 8966815, 8966815, GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantType.INS, "A", "AT", "g.chr19:8966815_8966816insT", "-", "c.43137_43138insA", "c.(43138-43140)tctfs", "p.S14380fs" },
                new Object[]{ "MUC16", 19, 8966650, 8966650, GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantType.INS, "C", "CG", "g.chr19:8966650_8966651insG", "-", "c.43302_43303insC", "c.(43303-43305)gggfs", "p.G14435fs" },
                // ==============================

                // ==============================
                // Test cases for FS INS issue:
                new Object[]{ "MUC16", 19, 8966813, 8966813, GencodeFuncotation.VariantClassification.IN_FRAME_INS, GencodeFuncotation.VariantType.INS, "A", "AGTGAGA", "g.chr19:8966813_8966814insGTGAGA", "-", "c.43139_43140insTCTCAC", "c.(43138-43140)tct>tcTCTCACt", "p.14380_14381insLT" }
        );
    }
}
