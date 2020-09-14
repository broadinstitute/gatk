package org.broadinstitute.hellbender.utils.dragstr;

import it.unimi.dsi.fastutil.doubles.AbstractDoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * User argument to specify a sequence of doubles with 3 values in the format {@code "start:step:limit"}.
 *
 * <p>The <b>start</b> is the first values in the sequence.
 * The <p>indicates the increment or decrement between consecutive values of the sequence</p>. And finally
 * the <b>limit</b> indicates when the sequence stops;
 * if increasing, when the next value exceeded the limit, if decreases when next value subceeds the limit</p>.
 * The sequence is inclusive as far as the limit.
 */
public final class DoubleSequence {

    private static final Pattern DOUBLE_PATTERN = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?");

    private final String input;

    private final String description;

    private final DoubleList values;

    /**
     * Currently the only supported format for the string descriptor is:
     * start:step:end.
     */
    public static final Pattern STEP_DECRIPTION_PATTERN = Pattern.compile(
            String.format("^(%s):(%s):(%s)$", DOUBLE_PATTERN, DOUBLE_PATTERN, DOUBLE_PATTERN));


    public DoubleSequence(final String descriptor) {
        Utils.nonNull(descriptor);
        final Matcher matcher = STEP_DECRIPTION_PATTERN.matcher(descriptor);
        if (!matcher.matches()) {
            throw new CommandLineException.BadArgumentValue("invalid double sequence specificatior: " + descriptor);
        }
        final double start = parseDouble(matcher.group(1), "start");
        final double step = parseDouble(matcher.group(2), "step");
        final double limit = parseDouble(matcher.group(3), "limit");

        final double[] valueArray;
        try {
            valueArray = MathUtils.doubles(start, limit, step);
        } catch (final RuntimeException ex) {
            throw new UserException.BadInput(ex.getMessage() + " (" + descriptor + ")");
        }
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

        input = descriptor;
    }

    private double parseDouble(final String str, final String role) {
        try {
            return Double.parseDouble(str);
        } catch (final NumberFormatException ex) {
            throw new UserException.BadInput(role + " is not a double formatted string: " + str, ex);
        }
    }

    /**
     * Returns an element of the sequence by its 0-base index.
     * @param index
     * @return
     */
    public double get(final int index) {
        Utils.validIndex(index, values.size());
        return values.get(index);
    }

    /**
     * Number of elements in the sequence.
     * @return 0 or greater.
     */
    public long size() {
        return values.size();
    }

    /**
     * Return the original spec of the sequence
     * @return
     */
    public String toString() {
        return input;
    }

    /**
     * Returns the minimum value of the sequence.
     * @return always a finite.
     */
    public double min() {
        final int size = values.size();
        if (size > 2) {
            return Math.min(values.get(0), values.get(size - 1));
        } else {
            return values.isEmpty() ? Double.NaN : values.get(0);
        }
    }

    /**
     * Returns the maximum value of the sequence
     * @return always a finite.
     */
    public double max() {
        final int size = values.size();
        if (size > 2) {
            return Math.max(values.get(0), values.get(size - 1));
        } else {
            return values.isEmpty() ? Double.NaN : values.get(0);
        }
    }

    /**
     * Returns the step in the sequence.
     * @return unspecified if the sequence has less than 2 elements, otherwise a finite.
     */
    public double step() {
        return values.size() < 2 ? Double.NaN : values.getDouble(1) - values.getDouble(0);
    }

    /**
     * Returns a new double array containing the elements of the sequence in the same order.
     * @return never {@code}, but perhaps an empty array.
     */
    public double[] toDoubleArray() {
        return values.toDoubleArray();
    }

    public String getDescription() {
        return description;
    }
}

