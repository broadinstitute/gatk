package org.broadinstitute.hellbender.utils.pairhmm;

import it.unimi.dsi.fastutil.doubles.AbstractDoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A sequence of doubles specificator that can be used as a user argument.
 */
public class DoubleSequence {

    private static final Pattern DOUBLE_PATTERN = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?");

    private final String input;

    private final String description;

    private DoubleList values;

    /**
     * Currently the only supported format for the string descriptor is:
     * start:step:end.
     */
    public static final Pattern STEP_DECRIPTION_PATTERN = Pattern.compile(
            String.format("^(%s):(%s):(%s)$", DOUBLE_PATTERN, DOUBLE_PATTERN, DOUBLE_PATTERN));


    public DoubleSequence(final String decriptor) {
        Utils.nonNull(decriptor);
        final Matcher matcher = STEP_DECRIPTION_PATTERN.matcher(decriptor);
        if (!matcher.matches()) {
            throw new CommandLineException.BadArgumentValue("invalid double sequence specificatior: " + decriptor);
        }
        final double left = Double.parseDouble(matcher.group(1));
        final double step = Double.parseDouble(matcher.group(2));
        final double right = Double.parseDouble(matcher.group(3));
        final double[] valueArray = MathUtils.doubles(left, right, step);
        values = new AbstractDoubleList() {

            @Override
            public int size() {
                return valueArray.length;
            }

            @Override
            public double getDouble(int index) {
                return valueArray[index];
            }

            @Override
            public double[] toDoubleArray() {
                return valueArray.clone();
            }
        };
        description = String.format("starting at %s to %s in increments/decrements of %s", matcher.group(1), matcher.group(3), matcher.group(2));

        input = decriptor;
    }

    public double get(final int index) {
        return values.get(index);
    }

    public long size() {
        return values.size();
    }

    public String toString() {
        return input;
    }

    public double[] toDoubleArray() {
        return values.toDoubleArray();
    }

    public String getDescription() {
        return description;
    }
}

