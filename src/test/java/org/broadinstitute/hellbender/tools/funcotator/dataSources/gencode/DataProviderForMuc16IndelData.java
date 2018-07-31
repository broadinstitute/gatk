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
                new Object[]{ "MUC16", 19, 9091890, 9091890, GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantType.INS, "T", "GCATC", "g.chr19:9091890_9091891insGCATC", "-", null, null, null }
        );
    }
}
