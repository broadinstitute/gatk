package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import static htsjdk.samtools.util.StringUtil.intValuesToString;

/**
 * In multiple locations we need to know what cycles are output, as of now we output all non-skip cycles, but rather than sprinkle
 * this knowledge throughout the parser code, instead OutputMapping provides all the data a client might want about the
 * cycles to be output including what ReadType they are.
 *
 * @author jburke@broadinstitute.org
 */
public class OutputMapping {
    /**
     * This class represents the mapping from Raw Cycles to TwoDIndices into output data structures and ClusterData.  This class
     * also contains ReadStructure.Substructure that describes which reads/cycles should be output.
     * <p>
     * For each cycle # (1-Based) there is a corresponding element in the cycleToOutputIndex array where
     * cycleToOutputIndex[cycle#].arrayIndex indicates the readNumber that cycle will be found on and cycleToOutputIndex[cycle#].elementIndex
     * indicates the array inde on that readNumber that the cycle occupies.  There are also various intermediate byte[][]
     * structures (in BclData, QseqReadData, etc...) where the top level array corresponds with the readNumber and the second-level
     * array corresponds with cycleNumber, cycleToOutputIndex is used to index into these arrays.
     */
    private final TwoDIndex[] cycleToOutputIndex;

    /**
     * The collection of ReadDescriptors and information about those read descriptors that describe all the
     * cycles that IlluminaDataProvider should output in a ClusterData object
     */
    private final ReadStructure.Substructure outputSubstructure;

    /**
     * The original read structure without any skips
     */
    private final ReadStructure outputReadStructure;

    /**
     * Create an OutputMapping from a readStructure, currently the outputSubstructure just references the readStructure.nonSkips
     * Substructure
     *
     * @param readStructure The readStructure for the given run that we want an OutputMapping for
     */
    public OutputMapping(final ReadStructure readStructure) {
        this.outputSubstructure = readStructure.nonSkips;
        this.cycleToOutputIndex = makeCycleToOutputIndexArray(readStructure);
        this.outputReadStructure = outputSubstructure.toReadStructure();
    }

    /**
     * @return The number of reads that should be found in the output clusterData
     */
    public int numOutputReads() {
        return outputSubstructure.length();
    }

    /**
     * @return An array of cycles in ascending order of all the cycles that should be output.
     */
    public int[] getOutputCycles() {
        return outputSubstructure.getCycles();
    }

    /**
     * @return An ordered array of lengths, where each element represents the size of output reads respectively
     */
    public int[] getOutputReadLengths() {
        return outputSubstructure.getDescriptorLengths();
    }

    /**
     * @return The total number of cycles that will be output
     */
    public int getTotalOutputCycles() {
        return outputSubstructure.getTotalCycles();
    }

    /**
     * @return An iterator over the read descriptors that describe the reads to be output
     */
    public Iterable<ReadDescriptor> getOutputDescriptors() {
        return outputSubstructure;
    }

    public ReadStructure getOutputReadStructure() {
        return outputReadStructure;
    }

    /**
     * An index into two dimensional arrays or other two dimensional structures where arrayIndex equals the index
     * into the top level structure and elementIndex is equal to the index into the bottom level structure
     */
    static class TwoDIndex {
        //index into the "outer structure" i.e. if we have an array T[][], we would index T[majorIndex][minorIndex]
        public final int majorIndex;

        //index into the "inner structure", see majorIndex
        public final int minorIndex;

        public TwoDIndex(final int majorIndex, final int minorIndex) {
            this.majorIndex = majorIndex;
            this.minorIndex = minorIndex;
        }

        @Override
        public String toString() {
            return "TwoDIndex(majorIndex == " + majorIndex + ", minorIndex == " + minorIndex + ")";
        }

        @Override
        public boolean equals(final Object thatObj) {
            if (thatObj == null || !(thatObj instanceof TwoDIndex)) {
                return false;
            }

            final TwoDIndex that = (TwoDIndex) thatObj;
            return this.majorIndex == that.majorIndex && this.minorIndex == that.minorIndex;
        }
    }

    /**
     * Create an array where each index corresponds to a cycle # (with Cycle 0 = null) and each element in
     * an index into a ClusterData's reads/bases and reads/qualities for the run described by readStructure
     *
     * @param readStructure The readStructure describing the run concerned
     * @return An array of TwoDArrayIndexes
     */
    private TwoDIndex[] makeCycleToOutputIndexArray(final ReadStructure readStructure) {
        int totalCycles = readStructure.totalCycles;
        final TwoDIndex[] cycleToOutputIndex = new TwoDIndex[totalCycles + 1];

        final int[] outputCycles = getOutputCycles();
        final int[] outputLengths = getOutputReadLengths();
        int outputCycleIndex = 0;
        int arrIndex = 0;
        int elementIndex = 0;
        for (int i = 1; i <= totalCycles && outputCycleIndex < outputCycles.length; i++) {
            if (outputCycles[outputCycleIndex] == i) {
                if (elementIndex >= outputLengths[arrIndex]) {
                    elementIndex = 0;
                    ++arrIndex;
                }

                cycleToOutputIndex[i] = new TwoDIndex(arrIndex, elementIndex);
                ++elementIndex;
                ++outputCycleIndex;
            }
        }

        if (outputCycleIndex != outputCycles.length) {
            throw new IlluminaParserException("Error in read structure outputCycles (" + intValuesToString(outputCycles) + ") and total cycles (" + totalCycles + ") OutputCycleIndex(" + outputCycleIndex + ")");
        }

        return cycleToOutputIndex;
    }
}
