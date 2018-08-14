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
                // These exercise a test case for issue 4712 (https://github.com/broadinstitute/gatk/issues/4712):
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, GencodeFuncotation.VariantType.INS, "T", "ATTCG", "g.chr19:9091890_9091891insATTCG", "-", null, null, null },
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantType.INS, "T", "CGCAT", "g.chr19:9091890_9091891insCGCA", "-", null, null, null },
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantType.INS, "T", "GCATC", "g.chr19:9091890_9091891insGCATC", "-", null, null, null },

                // Test cases for issue 5050:
                new Object[]{ "MUC16", 19, 8966646, 8966648, GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantType.DEL, "TTA", "T", "g.chr19:8966647_8966648delTA", "-", null, "c.e81+5", null },
                new Object[]{ "MUC16", 19, 8966818, 8966820, GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantType.DEL, "TGG", "T", "g.chr19:8966819_8966820delGG", "-", null, null, null },
                new Object[]{ "MUC16", 19, 8966813, 8966815, GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantType.DEL, "AGA", "A", "g.chr19:8966814_8966815delGA", "-", "c.43139_43140delCT", "c.(43138-43140)tctfs", "p.S14380fs" },
                new Object[]{ "MUC16", 19, 8966651, 8966653, GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantType.DEL, "ATC", "A", "g.chr19:8966652_8966653delTC", "-", "c.43301_43302delAT", "c.(43300-43302)gatfs", "p.D14434fs" }
        );
    }
}
