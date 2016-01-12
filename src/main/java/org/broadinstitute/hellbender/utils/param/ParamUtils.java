package org.broadinstitute.hellbender.utils.param;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.math3.exception.NotFiniteNumberException;
import org.apache.commons.math3.util.MathUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * This class should eventually be merged into Utils, which is in hellbender, and then this class should be deleted.
 *
 * Double.NaN will generally cause a check to fail.  Note that any comparison with a NaN yields false.
 *
 * Created by lichtens on 8/25/15.
 */
public class ParamUtils {
    private ParamUtils () {}

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     *
     * <p>Note that min can be greater than max, but that will guarantee that an exception is thrown.</p>
     *
     * @param val value to check
     * @param min minimum value for val
     * @param max maximum value for val
     * @param message the text message that would be pass to the exception thrown when {@code o == null}.
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double inRange(final double val, final double min, final double max, final String message) {
        if ((val >= min) && (val <= max)){
            return val;
        } else {
            throw new IllegalArgumentException(message);
        }
    }

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param min minimum value for val
     * @param max maximum value for val
     * @param message the text message that would be pass to the exception thrown when val gt min or val lt max.
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long inRange(final long val, final double min, final double max, final String message) {
        if ((val >= min) && (val <= max)){
            return val;
        } else {
            throw new IllegalArgumentException(message);
        }
    }

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param min minimum value for val (inclusive)
     * @param max maximum value for val (inclusive)
     * @param message the text message that would be pass to the exception thrown when val gt min or val lt max.
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long inRange(final long val, final long min, final long max, final String message) {
        if ((val >= min) && (val <= max)){
            return val;
        } else {
            throw new IllegalArgumentException(message);
        }
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long isPositiveOrZero(final long val, final String message) {
        if (!(val >= 0)){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isPositiveOrZero(final double val, final String message) {
        if (!(val >= 0)){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is greater than zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long isPositive(final long val, final String message) {
        if (!(val > 0)){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is greater than zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isPositive(final double val, final String message) {
        if (!(val > 0)){
            throw new IllegalArgumentException(message);
        }
        return val;
    }


    /**
     * Checks that the  input is not infinity nor NaN or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isFinite(final double val, final String message) {
        try {
            MathUtils.checkFinite(val);
        } catch (final NotFiniteNumberException ne) {
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * log the value in base b
     *
     * @param val -- The value to be logged.
     * @param base -- The base of the logging
     * @return the logged value.  If val is zero, returns Double.NEGATIVE_INFINITY.  If base is zero, returns -0.0.  If any input
     *  is negative, returns Double.NaN.
     */
    public static double logb(final double val, final double base) {
        return Math.log(val)/Math.log(base);
    }

    /**
     * log the value in base 2
     *
     * @param val -- The value to be logged.
     * @return the logged value.  If val is zero, returns Double.NEGATIVE_INFINITY.  If val
     *  is negative, returns Double.NaN.
     */
    public static double log2(final double val) {
        return Math.log(val)/Math.log(2.0);
    }

    /**
     *  Writes a list of strings to the given filename.  The file will be overwritten if it exists.
     *
     * @param stringIterator Strings that will be written (one to a line)
     * @param outputFile File to write
     */
    public static void writeStringListToFile(final Iterable<String> stringIterator, final File outputFile) {
        Utils.nonNull(stringIterator, "String list cannot be null");
        try (final FileWriter writer = new FileWriter(outputFile.getAbsolutePath())) {
            for (String str : stringIterator) {
                writer.write(str);
            }
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }
    }

    /**
     *  Writes a double array to the given filename.  The file will be overwritten if it exists.
     *
     * @param data doubles that will be written (one to a line)
     * @param outputFile File to write
     */
    public static void writeValuesToFile(final double[] data, final File outputFile) {
        Utils.nonNull(data, "input data cannot be null");
        try (final FileWriter writer = new FileWriter(outputFile.getAbsolutePath())) {
            for (double d : data) {
                writer.write(String.valueOf(d) + System.lineSeparator());
            }
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }
    }

    /**
     *  Reads a double array from the given filename.
     *
     * @param inputFile doubles that will be read (one to a line)
     */
    public static double[] readValuesFromFile(final File inputFile) {
        Utils.nonNull(inputFile, "input file cannot be null");
        final List<Double> listDouble = new ArrayList<>();
        try (final BufferedReader reader = new BufferedReader(new FileReader(inputFile.getAbsolutePath()))) {
            String line = reader.readLine();
            while (line != null) {
                listDouble.add(Double.parseDouble(line));
                line = reader.readLine();
            }
            return Doubles.toArray(listDouble);
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }
    }
}
