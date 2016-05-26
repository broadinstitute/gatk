package org.broadinstitute.hellbender.tools.exome;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Class whose instance represent copy-number-tri-state alleles.
 */
@SuppressWarnings("serial")
public final class CopyNumberTriStateAllele extends Allele {

    public static final CopyNumberTriStateAllele REF = new CopyNumberTriStateAllele(CopyNumberTriState.NEUTRAL, "N", true);
    public static final CopyNumberTriStateAllele DEL = new CopyNumberTriStateAllele(CopyNumberTriState.DELETION, "<DEL>", false);
    public static final CopyNumberTriStateAllele DUP = new CopyNumberTriStateAllele(CopyNumberTriState.DUPLICATION, "<DUP>", false);

    public static final String ALT_KEY = "ALT";
    public static final String DEL_VCF_DESCRIPTION = "Represents a deletion with respect to reference copy number";
    public static final String DUP_VCF_DESCRIPTION = "Represents a duplication with respect to reference copy number";

    public final CopyNumberTriState state;

    /**
     * List of all alleles as they should be included in VCF output files.
     * <p>
     *     This collection is unmodifiable.
     * </p>
     */
    public static final List<CopyNumberTriStateAllele> ALL_ALLELES =
            Collections.unmodifiableList(Arrays.asList(REF, DEL, DUP));

    /**
     * All allele list typed as a {@code List<Allele>} as required by some elements in the VCF framework.
     */
    @SuppressWarnings("unchecked")
    public static final List<Allele> PLAIN_ALL_ALLELES =
            (List<Allele>) (List) ALL_ALLELES;

    /**
     * List of all alleles as they should be included in VCF output files.
     * <p>
     *     This collection is unmodifiable.
     * </p>
     * <p>
     *     Notice that the collection type is {@link #Allele} rather than {@link #CopyNumberTriStateAllele}.
     *     This is done to be be able to use this list directly with some of the VCF framework class that expect
     *     this type of lists.
     * </p>
     */
    public static final List<CopyNumberTriStateAllele> ALTERNATIVE_ALLELES =
            Collections.unmodifiableList(Arrays.asList(DEL, DUP));

    private CopyNumberTriStateAllele(final CopyNumberTriState state, final String name, final boolean isRef) {
        super(Utils.nonNull(name), isRef);
        this.state = state;
    }

    /**
     * Adds the alternative allele ALT meta-data lines to a vcf-header.
     * @param header the header to add the lines to.
     * @throws IllegalArgumentException if {@code header} is {@code null}.
     */
    public static void addHeaderLinesTo(final VCFHeader header) {
        Utils.nonNull(header);
        header.addMetaDataLine(new VCFSimpleHeaderLine(ALT_KEY, DEL.getBaseString(), Utils.nonNull(DEL_VCF_DESCRIPTION)));
        header.addMetaDataLine(new VCFSimpleHeaderLine(ALT_KEY, DUP.getBaseString(), Utils.nonNull(DUP_VCF_DESCRIPTION)));
    }

    /**
     * Returns the value whose allele is the same as the one provided.
     * @param allele the query allele.
     * @return never {@code null}.
     * @throws IllegalArgumentException if the input allele is {@code null} or there is no
     * value with such an allele.
     */
    public static CopyNumberTriStateAllele valueOf(final Allele allele) {
        Utils.nonNull(allele, "the input allele cannot be null");
        if (allele.isReference()) {
            return REF;
        } else if (allele.getDisplayString().equals(DEL.getDisplayString())) {
            return DEL;
        } else if (allele.getDisplayString().equals(DUP.getDisplayString())) {
            return DUP;
        } else {
            throw new IllegalArgumentException("the input allele does not represent a valid copy-number-tristate allele: " + allele);
        }
    }

    /**
     * Returns the value whose state is the same as the one provided.
     * @param state the query state.
     * @return never {@code null}.
     * @throws IllegalArgumentException if the input state is {@code null} or there is no
     * value with such a state.
     */
    public static CopyNumberTriStateAllele valueOf(final CopyNumberTriState state) {
        switch (Utils.nonNull(state)) {
            case NEUTRAL: return REF;
            case DELETION: return DEL;
            case DUPLICATION: return DUP;
            default:
                throw new GATKException.ShouldNeverReachHereException("unexpected/unsupported copy-number-tristate state: " + state);
        }
    }

    /**
     * Returns the index of the allele in the standard order (its {@link #ALL_ALLELES}' index).
     * @return 0 to the number of alleles -1.
     */
    public int index() {
        return ALL_ALLELES.indexOf(this);
    }
}
