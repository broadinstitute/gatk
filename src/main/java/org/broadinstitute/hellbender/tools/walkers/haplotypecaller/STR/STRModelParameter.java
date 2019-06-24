package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import java.util.Arrays;
import java.util.stream.Stream;

/**
 * Created by valentin on 11/29/16.
 */
public class STRModelParameter {

    public enum Name {
        PI, TAU, DEL, INS;

        public boolean isAValidName(final String name) {
            if (name == null) {
                throw new IllegalArgumentException("the input name cannot be null");
            } else if (Stream.of(values()).anyMatch(x -> x.name().equalsIgnoreCase(name))) {
                return true;
            } else {
                return false;
            }
        }
    }

    public final double maximum;

    public final double intercept;

    public final double repeatCountCoefficient;

    private volatile double[] values;

    STRModelParameter(final double maximum, final double intercept, final double repeatCountCoefficient) {
        this.maximum = maximum;
        this.intercept = intercept;
        this.repeatCountCoefficient = repeatCountCoefficient;
        values = new double[0];
    }

    public double valueFor(final int inputRepeatCount) {
        final int repeatCountDifference = Math.abs(inputRepeatCount);
        if (repeatCountDifference >= values.length) {
            values = extendsValues(values, repeatCountDifference << 1);
        }
        return values[inputRepeatCount];
    }

    private double[] extendsValues(final double[] previousValues, final int extendTo) {
        if (previousValues.length > extendTo) {
            return previousValues;
        } else {
            final double[] result = Arrays.copyOfRange(previousValues, 0, extendTo + 1);
            for (int i = previousValues.length; i < result.length; ++i) {
                result[i] = maximum / (1.0 + Math.exp(- intercept - repeatCountCoefficient * i ));
            }
            return result;
        }
    }
}
