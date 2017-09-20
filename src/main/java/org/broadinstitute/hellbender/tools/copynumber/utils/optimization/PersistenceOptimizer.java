package org.broadinstitute.hellbender.tools.copynumber.utils.optimization;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given 1-dimensional data, finds all local minima sorted by decreasing topological persistence.
 * This can be useful for finding the most significant local minima in noisy data and can be used
 * to smooth out less significant local minima that may be induced by noise.
 * Algorithm is adapted from
 * <a href="https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html">
 *     https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html</a>.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PersistenceOptimizer {
    private static final int NO_COLOR = -1;

    //represents a connected component bounded by leftIndex and rightIndex with minimum value minValue at minIndex
    //(may also contain additional local minima)
    private static final class Component {
        private int leftIndex;
        private int rightIndex;
        private final int minIndex;
        private final double minValue;

        //Components always initially contain only a local minimum
        private Component(final int minIndex,
                          final double minValue) {
            this.leftIndex = minIndex;
            this.rightIndex = minIndex;
            this.minIndex = minIndex;
            this.minValue = minValue;
        }
    }

    private static final class ExtremaPair {
        private final int minIndex;
        private final int maxIndex;                     //these could be used to find local maxima, if desired
        private final double persistence;

        private ExtremaPair(final double[] data,
                            final int datumIndex1,
                            final int datumIndex2) {
            if (data[datumIndex1] > data[datumIndex2]) {
                minIndex = datumIndex2;
                maxIndex = datumIndex1;
            } else if (data[datumIndex1] < data[datumIndex2]) {
                minIndex = datumIndex1;
                maxIndex = datumIndex2;
            } else if (datumIndex1 < datumIndex2) {     //if values are equal, sort by index
                minIndex = datumIndex1;
                maxIndex = datumIndex2;
            } else {
                minIndex = datumIndex2;
                maxIndex = datumIndex1;
            }
            this.persistence = data[maxIndex] - data[minIndex];
        }
    }

    private final double[] data;
    private final List<Integer> minimaIndices;
    private final List<Double> persistences;

    /**
     * Identifies the local minima of {@code data} based on topological persistence upon construction.
     * Note that for a region where the data takes constant values that occurs either to the right
     * of a local maximum or at the beginning of the data, the leftmost point in the region
     * is considered a local minimum.
     * @param data  values of a one-dimensional function evaluated at points in left-to-right order
     */
    public PersistenceOptimizer(final double[] data) {
        Utils.nonNull(data);
        Utils.validateArg(data.length > 0, "Data must contain at least one element.");

        this.data = Arrays.copyOf(data, data.length);   //defensive copy
        final List<Integer> sortedIndices = IntStream.range(0, this.data.length).boxed()
                .sorted(Comparator.comparingDouble(i -> this.data[i]))
                .collect(Collectors.toList());
        final List<ExtremaPair> extremaPairs = findExtremaPairs(this.data, sortedIndices);
        minimaIndices = extremaPairs.stream().map(p -> p.minIndex).collect(Collectors.toList());
        minimaIndices.add(0, sortedIndices.get(0));     //add index of global minimum
        persistences = extremaPairs.stream().map(p -> p.persistence).collect(Collectors.toList());
        persistences.add(0, this.data[sortedIndices.get(sortedIndices.size() - 1)] - this.data[sortedIndices.get(0)]);  //add global persistence (maximum value - minimum value)
    }

    /**
     * Returns an unmodifiable list of the indices of the local minima, sorted first by decreasing topological persistence
     * and then by increasing index.  The first element is the index of the global minimum.
     */
    public List<Integer> getMinimaIndices() {
        return Collections.unmodifiableList(minimaIndices);
    }

    /**
     * Returns the corresponding unmodifiable list of the topological persistences of the local minima given by
     * {@link PersistenceOptimizer#getMinimaIndices()}.
     */
    public List<Double> getPersistences() {
        return Collections.unmodifiableList(persistences);
    }

    //find extrema pairs via a watershed algorithm;
    //note that sortedIndices must be the stably sorted indices of data
    private static List<ExtremaPair> findExtremaPairs(final double[] data,
                                                      final List<Integer> sortedIndices) {
        final List<Component> components = new ArrayList<>(data.length);
        final int[] colors = new int[data.length];
        Arrays.fill(colors, NO_COLOR);
        final List<ExtremaPair> extremaPairs = new ArrayList<>(data.length);

        if (data.length == 1) {
            return Collections.emptyList();
        } else {
            for (final int index : sortedIndices) {                             //watershed algorithm starts from lowest data point and iterates upward
                if (index == 0) {                                               //first point has no left neighbor
                    if (colors[index + 1] == NO_COLOR) {                        //local minimum, create new component
                        createComponent(data, components, colors, index);
                    } else {                                                    //extend component on right to include this point
                        extendComponent(components, colors, colors[index + 1], index);
                    }
                    continue;
                } else if (index == data.length - 1) {                          //last point has no right neighbor
                    if (colors[index - 1] == NO_COLOR) {                        //local minimum, create new component
                        createComponent(data, components, colors, index);
                    } else {                                                    //extend component on left to include this point
                        extendComponent(components, colors, colors[index - 1], index);
                    }
                    continue;
                }

                //otherwise, look both left and right
                final int leftColor = colors[index - 1];
                final int rightColor = colors[index + 1];
                if (leftColor == NO_COLOR && rightColor == NO_COLOR) {          //local minimum, create new component
                    createComponent(data, components, colors, index);
                } else if (leftColor != NO_COLOR && rightColor == NO_COLOR) {   //extend component on left to include this point
                    extendComponent(components, colors, leftColor, index);
                } else if (leftColor == NO_COLOR && rightColor != NO_COLOR) {   //extend component on right to include this point
                    extendComponent(components, colors, rightColor, index);
                } else {                                                        //local maximum, merge components
                    if (components.get(rightColor).minValue < components.get(leftColor).minValue) {
                        extremaPairs.add(new ExtremaPair(data, components.get(leftColor).minIndex, index));
                    } else {
                        extremaPairs.add(new ExtremaPair(data, components.get(rightColor).minIndex, index));
                    }
                    mergeComponents(components, colors, leftColor, rightColor, index);
                }
            }
        }

        //return pairs sorted by descending persistence (explicitly-typed lambda required here for use of reversed() to compile)
        return extremaPairs.stream().sorted(Comparator.comparingDouble((ExtremaPair p) -> p.persistence).reversed()).collect(Collectors.toList());
    }

    //create a new component and add a new color when a local minimum is found at index
    private static void createComponent(final double[] data,
                                        final List<Component> components,
                                        final int[] colors,
                                        final int index) {
        colors[index] = components.size();
        components.add(new Component(index, data[index]));
    }

    //given components and colors, extend the component at componentIndex to include the point at index
    //and assign the color of that component to this point
    private static void extendComponent(final List<Component> components,
                                        final int[] colors,
                                        final int componentIndex,
                                        final int index) {
        if (index + 1 == components.get(componentIndex).leftIndex) {
            components.get(componentIndex).leftIndex = index;
        } else if (index - 1 == components.get(componentIndex).rightIndex) {
            components.get(componentIndex).rightIndex = index;
        }
        colors[index] = componentIndex;
    }

    //merge components by assigning the color of the component with the lower minimum value
    //to the boundary of the merged component when a local maximum is encountered at index
    private static void mergeComponents(final List<Component> components,
                                        final int[] colors,
                                        final int leftColor,
                                        final int rightColor,
                                        final int index) {
        final int indexToKeep;
        final int indexToMerge;
        if (components.get(leftColor).minValue < components.get(rightColor).minValue) {
            indexToKeep = leftColor;
            indexToMerge = rightColor;
        } else if (components.get(leftColor).minValue > components.get(rightColor).minValue) {
            indexToKeep = rightColor;
            indexToMerge = leftColor;
        } else if (leftColor < rightColor) {       //if values are equal, sort by index
            indexToKeep = leftColor;
            indexToMerge = rightColor;
        } else {
            indexToKeep = rightColor;
            indexToMerge = leftColor;
        }
        colors[components.get(indexToMerge).leftIndex] = indexToKeep;
        colors[components.get(indexToMerge).rightIndex] = indexToKeep;
        if (components.get(indexToKeep).minIndex > components.get(indexToMerge).minIndex) {
            components.get(indexToKeep).leftIndex = components.get(indexToMerge).leftIndex;
        } else {
            components.get(indexToKeep).rightIndex = components.get(indexToMerge).rightIndex;
        }
        colors[index] = colors[index - 1];          //local maximum takes color of component to left by convention
    }
}