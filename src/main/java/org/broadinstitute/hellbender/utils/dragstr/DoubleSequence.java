package org.broadinstitute.hellbender.utils.dragstr;

import it.unimi.dsi.fastutil.doubles.AbstractDoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A sequence of doubles specificator that can be used as a user argument.
 * <p>
 *     Currently it support the format start:step:end where the seqeuence generated is
 *     start, start + step, start + 2*step,..., start + N*step where N is the largerst N that keep the number less
 *     than end.
 * </p>
 */
public final class DoubleSequence {

    /**
     * Double literal pattern.
     */
    private static final Pattern DOUBLE_PATTERN = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?");

    /**
     * Currently the only supported format for the string descriptor is:
     * start:step:end.
     */
    private static final Pattern STEP_DECRIPTION_PATTERN = Pattern.compile(
            String.format("^(%s):(%s):(%s)$", DOUBLE_PATTERN, DOUBLE_PATTERN, DOUBLE_PATTERN));

    /**
     * Holds a copy of the spec that was used to generate this sequence.
     */
    private final String spec;

    /**
     * The sequence double values.
     */
    private DoubleList values;
    
    /**
     * Creates a sequence given its spec.
     * @param spec never {@code null}.
     */
    public DoubleSequence(final String spec) {
        Utils.nonNull(spec);
        final Matcher matcher = STEP_DECRIPTION_PATTERN.matcher(spec);
        if (!matcher.matches()) {
            throw new CommandLineException.BadArgumentValue("invalid double sequence specificatior: " + spec);
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

        this.spec = spec;
    }

    public double get(final int index) {
        return values.get(index);
    }

    public long size() {
        return values.size();
    }

    public String toString() {
        return spec;
    }

    public double[] toDoubleArray() {
        return values.toDoubleArray();
    }
}

