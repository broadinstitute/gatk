package org.broadinstitute.hellbender.tools.recalibration.covariates;

public class RepeatLengthCovariate extends RepeatCovariate {

    protected String getCovariateValueFromUnitAndLength(final byte[] repeatFromUnitAndLength, final int repeatLength) {
        return String.format("%d", repeatLength);
    }

    @Override
    public int maximumKeyValue() {
        // max possible values of covariate: for repeat unit, length is up to MAX_STR_UNIT_LENGTH,
        // so we have 4^MAX_STR_UNIT_LENGTH * MAX_REPEAT_LENGTH possible values
        return (1+MAX_REPEAT_LENGTH);
    }

}
