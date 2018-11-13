package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import java.util.Arrays;
import java.util.List;

/**
 * Hold a data provider that creates trouble variants.
 * Created by jonn on 11/9/18.
 */
public class DataProviderForTroubleVariants {

    public static List<Object[]> provideSymbolicAllelesAndMaskedBasesForHg38() {
        return Arrays.asList(

                // Symbolic alt allele:
        new Object[] {"NPB",
                        17, 81902450, 81902450,
                        GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE,
                        GencodeFuncotation.VariantType.NA,
                        "C",
                        "<DEL>",
                        null,
                        "+",
                        null,
                        null,
                        null},
                // Symbolic alt allele:
                new Object[] {"NPB",
                        17, 81902450, 81902450,
                        GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE,
                        GencodeFuncotation.VariantType.NA,
                        "C",
                        "<DUP>",
                        null,
                        "+",
                        null,
                        null,
                        null},
                // Masked Ref allele + Symbolic alt allele:
                new Object[] {"NPB",
                        17, 81902450, 81902450,
                        GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE,
                        GencodeFuncotation.VariantType.NA,
                        "N",
                        "<DEL>",
                        null,
                        "+",
                        null,
                        null,
                        null},
                // Masked Ref allele + Symbolic alt allele:
                new Object[] {"NPB",
                        17, 81902450, 81902450,
                        GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE,
                        GencodeFuncotation.VariantType.NA,
                        "N",
                        "<DUP>",
                        null,
                        "+",
                        null,
                        null,
                        null},
                // Masked Ref allele:
                new Object[] {"NPB",
                        17, 81902450, 81902450,
                        GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE,
                        GencodeFuncotation.VariantType.SNP,
                        "N",
                        "A",
                        null,
                        "+",
                        null,
                        null,
                        null},
                // Masked Alt allele:
                new Object[] {"NPB",
                        17, 81902450, 81902450,
                        GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE,
                        GencodeFuncotation.VariantType.SNP,
                        "C",
                        "N",
                        null,
                        "+",
                        null,
                        null,
                        null},
                // Masked Ref allele + Symbollic alt allele in IGR:
                new Object[] {null,
                        17, 81971500, 81971500,
                        GencodeFuncotation.VariantClassification.IGR,
                        GencodeFuncotation.VariantType.NA,
                        "N",
                        "<DUP>",
                        null,
                        null,
                        null,
                        null,
                        null},
                // Masked Ref allele in IGR:
                new Object[] {null,
                        17, 81971500, 81971500,
                        GencodeFuncotation.VariantClassification.IGR,
                        GencodeFuncotation.VariantType.SNP,
                        "N",
                        "A",
                        "g.chr17:81971500N>A",
                        null,
                        null,
                        null,
                        null},
                // Masked Alt allele in IGR:
                new Object[] {null,
                        17, 81971500, 81971500,
                        GencodeFuncotation.VariantClassification.IGR,
                        GencodeFuncotation.VariantType.SNP,
                        "C",
                        "N",
                        "g.chr17:81971500C>N",
                        null,
                        null,
                        null,
                        null}
        );
    }

}
