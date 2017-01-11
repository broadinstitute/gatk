package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.interfaces.AlleleMetadataProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;

import javax.annotation.Nonnull;
import java.util.*;

/**
 * Three-state copy number hidden states enum.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public enum CopyNumberTriState implements AlleleMetadataProducer, CallStringProducer, ScalarProducer {

    /**
     * Represents a loss of one or more copies in a region.
     */
    DELETION("-", Allele.create("<DEL>", false), 0.5),

    /**
     * Represents not a loss nor a gain of any copy in a region.
     */
    NEUTRAL("0", Allele.create("N", true), 1.0),

    /**
     * Represents a gain of one or more copies in a region.
     */
    DUPLICATION("+", Allele.create("<DUP>", false), 1.5);

    /**
     * For making VCF header
     */
    public static final Map<String, String> ALT_CALL_STRING_TO_VCF_DESCRIPTION = ImmutableMap.of(
            "-", "Represents a deletion with respect to reference copy number",
            "+", "Represents a duplication with respect to reference copy number");
    public static final String ALT_KEY = "ALT";

    /**
     * Unmodifiable list alleles where the first allele is the reference allele
     * (corresponding {@link #NEUTRAL}). Other states alternative alleles follow
     * in ascending order based on their ordinal.
     */
    public static List<Allele> ALL_ALLELES =
            Collections.unmodifiableList(Arrays.asList(NEUTRAL.allele, DELETION.allele, DUPLICATION.allele));

    /**
     * Unmodifiable list of alternative alleles order by their corresponding state ordinal.
     */
    public static List<Allele> ALTERNATIVE_ALLELES = ALL_ALLELES.subList(1, ALL_ALLELES.size());

    public final String callString;

    public final Allele allele;

    public final double copyRatio;

    CopyNumberTriState(final String callString, final Allele allele, final double copyRatio) {
        this.callString = callString;
        this.allele = allele;
        this.copyRatio = copyRatio;
    }

    /**
     * Returns the state given a call string.
     *
     * @param callString the query call string.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code callString} is either a {@code null} or an invalid
     *  call string.
     */
    public static CopyNumberTriState fromCallString(final String callString) {
        Utils.nonNull(callString);
        if (callString.equals(DELETION.callString)) {
            return DELETION;
        } else if (callString.equals(NEUTRAL.callString)) {
            return NEUTRAL;
        } else if (callString.equals(DUPLICATION.callString)) {
            return DUPLICATION;
        } else {
            throw new IllegalArgumentException(String.format("'%s' is not a valid call name", callString));
        }
    }

    /**
     * Returns a call string given the state.
     * @return non-null call string
     */
    @Override
    public String getCallString() {
        return callString;
    }

    @Override
    public double getScalar() {
        return copyRatio;
    }

    @Override
    public Allele toAllele() {
        return allele;
    }

    @Override
    public void addHeaderLineTo(@Nonnull VCFHeader header) {
        if (ALT_CALL_STRING_TO_VCF_DESCRIPTION.containsKey(callString)) {
            header.addMetaDataLine(new VCFSimpleHeaderLine(ALT_KEY, allele.getBaseString(),
                    ALT_CALL_STRING_TO_VCF_DESCRIPTION.get(callString)));
        }
    }
}
