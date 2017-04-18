package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.utils.hmm.interfaces.AlleleMetadataProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;

/**
 * This class represents integer copy number states. It is intended to be used as states of a
 * hidden Markov model.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberState implements AlleleMetadataProducer, CallStringProducer, ScalarProducer,
        Serializable {

    private static final long serialVersionUID = -3641032600858976498L;

    /**
     * Integer value of the represented copy number state
     */
    private final int copyNumber;

    /**
     * A string representation of this copy number state (used for creating human-readable copy number call lists)
     */
    private final String callString;

    /**
     * An allele representation of this copy number state (used for VCF creation)
     */
    private final Allele allele;

    public IntegerCopyNumberState(final int copyNumber) {
        this.copyNumber = ParamUtils.isPositiveOrZero(copyNumber, "The integer copy number state" +
                " must be non-negative");
        this.callString = toCallString(copyNumber);
        allele = Allele.create(toAlleleString(copyNumber));
    }

    @Override
    public Allele toAllele() {
        return allele;
    }

    /**
     * TODO github/gatk-protected issue #855 -- this is required for VCF creation
     * @param header an instance of {@link VCFHeader}
     */
    @Override
    public void addHeaderLineTo(@Nonnull VCFHeader header) {
        throw new UnsupportedOperationException("github/gatk-protected issue #855");
    }

    @Override
    public String getCallString() {
        return callString;
    }

    @Override
    public double getScalar() {
        return copyNumber;
    }

    public int getCopyNumber() { return copyNumber; }

    private static String toCallString(final int copyNumberState) {
        return String.valueOf(copyNumberState);
    }

    private static String toAlleleString(final int copyNumberState) {
        return "<" + String.valueOf(copyNumberState) + ">";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        return copyNumber == ((IntegerCopyNumberState) o).copyNumber;
    }

    @Override
    public int hashCode() {
        return copyNumber;
    }
}
