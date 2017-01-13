package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.math3.util.FastMath;
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

    private final int copyNumber;

    private final String callString;

    private final Allele allele;

    /**
     * These values are often required, we precompute them
     */
    private final double logCopyNumber, logCopyNumberSquared;

    public IntegerCopyNumberState(final int copyNumber) {
        this.copyNumber = ParamUtils.isPositiveOrZero(copyNumber, "The integer copy number state" +
                " must be non-negative");
        this.logCopyNumber = FastMath.log(copyNumber);
        this.logCopyNumberSquared = logCopyNumber * logCopyNumber;
        this.callString = toCallString(copyNumber);
        allele = Allele.create(toAlleleString(copyNumber));
    }

    @Override
    public Allele toAllele() {
        return allele;
    }

    @Override
    public void addHeaderLineTo(@Nonnull VCFHeader header) {

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

    public double getLogCopyNumber() { return logCopyNumber; }

    public double getLogCopyNumberSquared() { return logCopyNumberSquared; }

    public static String toCallString(final int copyNumberState) {
        return String.valueOf(copyNumberState);
    }

    public static String toAlleleString(final int copyNumberState) {
        return "<" + String.valueOf(copyNumberState) + ">";
    }

    /**
     * The equality comparison is only made based on the integer copy number; the boolean reference status of the
     * state is ignored
     *
     * @param o another object
     * @return boolean
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        IntegerCopyNumberState that = (IntegerCopyNumberState) o;

        return copyNumber == that.copyNumber;

    }

    @Override
    public int hashCode() {
        return copyNumber;
    }
}
