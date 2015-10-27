package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

/**
 * Represents a named dataset consisting of numeric datapoints.
 * @param <N>   type of the datapoints
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Data<N extends Number> {
    private final String name;
    private final List<N> data;

    /**
     * Constructs a Data object given a name and a list of numeric datapoints.
     * @param name  dataset name
     * @param data  List of numeric datapoints
     * @throws IllegalArgumentException if {@code data} is empty, 
     *                                  or if either {@code name} or {@code data} are {@code null}
     */
    public Data(final String name, final List<N> data) {
        Utils.nonNull(name, "The name of the dataset cannot be null.");
        Utils.nonNull(data, "The dataset cannot be null.");
        if (data.size() == 0) {
            throw new IllegalArgumentException("The dataset cannot be empty.");
        }
        this.name = name;
        this.data = new ArrayList<>(data);
    }

    /**
     * Constructs a Data object given a name, a file of numeric datapoints, and a parsing function.
     * @param name  dataset name
     * @param file  file of numeric datapoints
     *              a single column of datapoints with one datapoint per line is assumed; no headers, etc. are allowed
     * @param parse function to parse each line in file to its numeric value
     */
    public Data(final String name, final File file, final Function<String, N> parse) {
        this(name, loadData(file, parse));
    }

    /**
     * Returns the dataset name.
     * @return  dataset name
     */
    protected String name() {
        return name;
    }

    /**
     * Returns an unmodifiable view of the List of datapoints held internally.
     * @return  an unmodifiable view of the List of datapoints held internally
     */
    public List<N> values() {
        return Collections.unmodifiableList(data);
    }

    /**
     * Returns the value at the specified index in the List of datapoints held internally.
     * @param index index
     * @return      value at {@code index} in the List of datapoints held internally
     */
    public N value(final int index) {
        return data.get(index);
    }

    /**
     * Returns the size of the list of the datapoints held internally.
     * @return      size of the list of the datapoints held internally
     */
    public int size() {
        return data.size();
    }

    /**
     * Helper function for loading numeric datapoints from a file.  A single column of datapoints
     * with one datapoint per line is assumed; no headers, etc. are allowed.
     * @param file  file containing datapoints
     * @param parse function to parse each line in file to its numeric value
     * @param <T>   type of datapoints
     * @return      List of datapoints
     */
    private static <T> List<T> loadData(final File file, final Function<String, T> parse) {
        final List<T> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                list.add(parse.apply(line));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
        return list;
    }
}
