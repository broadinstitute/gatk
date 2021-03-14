package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import org.broadinstitute.hellbender.tools.funcotator.FuncotatorConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Hold a data provider that creates trouble variants.
 * Created by jonn on 11/9/18.
 */
public class DataProviderForTroubleVariants {

    public static List<Object[]> provideTP53AllelesForOrderTestsHg19() {
        return Collections.singletonList(
                new Object[]{ "TP53",
                        17, 7578492, 7578492,
                        GencodeFuncotation.VariantClassification.NONSENSE,
                        GencodeFuncotation.VariantType.SNP,
                        "C",
                        "T",
                        "g.chr17:7578492C>T",
                        "-",
                        "c.438G>A",
                        "c.(436-438)tgG>tgA",
                        "p.W146*" }
        );
    }

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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
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
                        FuncotatorConstants.MASKED_ANY_BASE_STRING,
                        "g.chr17:81971500C>N",
                        null,
                        null,
                        null,
                        null}
        );
    }

}
