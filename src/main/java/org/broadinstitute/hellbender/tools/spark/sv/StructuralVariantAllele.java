package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.base.Functions;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 4/27/17.
 */
public enum StructuralVariantAllele {
    INS, DEL, DUP, INV;

    private final Allele allele;

    @SuppressWarnings("RedundantTypeArguments")
    private static final Map<String, StructuralVariantAllele> instanceByName =
            Stream.of(values()).collect(Collectors.toMap(Enum<StructuralVariantAllele>::name, v -> v));

    /**
     * Creates a new {@link StructuralVariantAllele} instance.
     */
    StructuralVariantAllele() {
        allele = Allele.create("<" + name() + ">", false);
    }

    /**
     * Returns the {@link StructuralVariantAllele} instance that corresponds to a given allele.
     * @param allele the query allele.
     * @throws NullPointerException if the input allele is {@code null}.
     * @throws IllegalArgumentException if the input allele does not represent a known variant allele.
     * @return never {@code null}.
     */
    public static StructuralVariantAllele valueOf(final Allele allele) {
        if (allele == null) {
            throw new IllegalArgumentException();
        } else if (!allele.isSymbolic()) {
            throw new IllegalArgumentException();
        } else {
            final String text = allele.getDisplayString();
            if (!text.startsWith("<") || !text.endsWith(">")) {
                throw new IllegalArgumentException("unexpected symbolic allele name: " + text);
            } else {
                final StructuralVariantAllele result = instanceByName.get(text.substring(1, text.length() - 1));
                if (result == null) {
                    throw new IllegalArgumentException();
                } else {
                    return result;
                }
            }
        }
    }

    /**
     * Checks whether an allele seems to be a structural one.
     * @param allele the allele.
     * @return never {@code null}.
     */
    public static boolean isStructural(final Allele allele) {
        if (allele == null || !allele.isSymbolic()) {
            return false;
        } else {
            final String text = allele.getDisplayString();
            if (!text.startsWith("<") || !text.endsWith(">")) {
                return false;
            } else {
                return instanceByName.containsKey(text.substring(1, text.length() - 1));
            }
        }
    }

    /**
     * Returns a allele instance that represents this
     *
     * @return never {@code null}.
     */
    public Allele allele() {
        return allele;
    }
}

