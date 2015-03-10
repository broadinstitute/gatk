package org.broadinstitute.hellbender.tools.recalibration.covariates;

public class RepeatUnitCovariate extends RepeatCovariate {

    protected String getCovariateValueFromUnitAndLength(final byte[] repeatFromUnitAndLength, final int repeatLength) {
        return new String(repeatFromUnitAndLength);
    }

    @Override
    public int maximumKeyValue() {
        // max possible values of covariate: for repeat unit, length is up to MAX_STR_UNIT_LENGTH,
        // so we have 4^MAX_STR_UNIT_LENGTH * MAX_REPEAT_LENGTH possible values
        return (1<<(2*MAX_STR_UNIT_LENGTH))  +1;
    }
}
