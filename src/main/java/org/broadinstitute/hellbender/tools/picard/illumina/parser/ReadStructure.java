package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CoordMath;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static htsjdk.samtools.util.CoordMath.getEnd;
import static java.lang.Integer.parseInt;
import static java.util.Collections.unmodifiableList;
import static java.util.regex.Pattern.compile;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadType.valueOf;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadType.values;

/**
 * Describes the intended logical output structure of clusters of an Illumina run.
 * (e.g. If the input data consists of 80 base
 * clusters and we provide a read structure of "36T8B36T" then those bases should be split into 3 reads:
 * read one should be 36 cycles of template,
 * read two should be 8 cycles of barcode,
 * read three should be another 36 cycle template read.)
 * Note: In future releases, ReadStructures will be specified by clients of IlluminaDataProvider(currently
 * read structures are detected by IlluminaDataProviderFactory via the structure of QSeq files). When working with
 * QSeq formats, while the individual reads need not fall on QSeq end file boundaries the total number of cycles
 * should equal the total number of cycles found in all end files for one tile.  (e.g. if I have 80 reads and
 * 3 end files per tile, those end files should have a total of 80 reads in them regardless of how many reads
 * appear in each individual file)
 *
 * @author jburke@broadinstitute.org
 */
public final class ReadStructure {
    public static final String PARAMETER_DOC =
            "A description of the logical structure of clusters in an Illumina Run, i.e. a description of the structure IlluminaBasecallsToSam " +
                    "assumes the  data to be in. It should consist of integer/character pairs describing the number of cycles and the type of those " +
                    "cycles (B for Barcode, T for Template, and S for skip).  E.g. If the input data consists of 80 base clusters and we provide a " +
                    "read structure of \"36T8B8S28T\" then, before being converted to SAM records those bases will be split into 4 reads where " +
                    "read one consists of 36 cycles of template, read two consists of 8 cycles of barcode, read three will be an 8 base read of " +
                    "skipped cycles and read four is another 28 cycle template read.  The read consisting of skipped cycles would NOT be included " +
                    "in output SAM/BAM file read groups.";
    public final List<ReadDescriptor> descriptors;
    public final int totalCycles;
    public final int[] readLengths;

    public final Substructure barcodes;
    public final Substructure templates;
    public final Substructure skips;

    //nonSkips include barcode and template indices in the order they appear in the descriptors list
    public final Substructure nonSkips;

    /**
     * Characters representing valid ReadTypes
     */
    private static final String ValidTypeChars;

    /**
     * ValidTypeChars except characters are separated by commas
     */
    private static final String ValidTypeCharsWSep;

    static {
        String validTypes = "";
        String vtWSep = "";

        boolean written = false;
        for (ReadType rt : values()) {
            if (written) {
                vtWSep += ",";
            }

            validTypes += rt.name();
            vtWSep += rt.name();
        }

        ValidTypeChars = validTypes;
        ValidTypeCharsWSep = vtWSep;
    }

    private static final String ReadStructureMsg = "Read structure must be formatted as follows: " +
            "<number of bases><type><number of bases><type>...<number of bases> where number of bases is a " +
            "positive (NON-ZERO) integer and type is one of the following characters " + ValidTypeCharsWSep +
            " (e.g. 76T8B68T would denote a paired-end run with a 76 base first end an 8 base barcode followed by a 68 base second end).";
    private static final Pattern FullPattern = compile("^((\\d+[" + ValidTypeChars + "]{1}))+$");
    private static final Pattern SubPattern = compile("(\\d+)([" + ValidTypeChars + "]{1})");

    /**
     * Copies collection into descriptors (making descriptors unmodifiable) and then calculates relevant statistics about descriptors.
     *
     * @param collection A collection of ReadDescriptors that describes this ReadStructure
     */
    public ReadStructure(final List<ReadDescriptor> collection) {
        if (collection.size() == 0) { //If this changes, change hashcode
            throw new IllegalArgumentException("ReadStructure does not support 0 length clusters!");
        }

        final List<Range> allRanges = new ArrayList<Range>(collection.size());
        this.descriptors = unmodifiableList(collection);
        int cycles = 0;

        final List<Integer> nonSkipIndicesList = new ArrayList<Integer>();
        final List<Integer> barcodeIndicesList = new ArrayList<Integer>();
        final List<Integer> templateIndicesList = new ArrayList<Integer>();
        final List<Integer> skipIndicesList = new ArrayList<Integer>();
        readLengths = new int[collection.size()];

        int currentCycleIndex = 0;   // Current cycle in the entire read structure
        int descIndex = 0;
        for (final ReadDescriptor desc : descriptors) {
            if (desc.length == 0 || desc.length < 0) {
                throw new IllegalArgumentException("ReadStructure only supports ReadDescriptor lengths > 0, found(" + desc.length + ")");
            }

            final int endIndexOfRange = getEnd(currentCycleIndex, desc.length);
            allRanges.add(new Range(currentCycleIndex, endIndexOfRange));
            currentCycleIndex = endIndexOfRange + 1;

            readLengths[descIndex] = desc.length;
            cycles += desc.length;
            switch (desc.type) {
                case B:
                    nonSkipIndicesList.add(descIndex);
                    barcodeIndicesList.add(descIndex);
                    break;
                case T:
                    nonSkipIndicesList.add(descIndex);
                    templateIndicesList.add(descIndex);
                    break;
                case S:
                    skipIndicesList.add(descIndex);
                    break;

                default:
                    throw new IllegalArgumentException("Unsupported ReadType (" + desc.type + ") encountered by IlluminaRunConfiugration!");
            }
            ++descIndex;
        }

        this.totalCycles = cycles;
        this.barcodes = new Substructure(barcodeIndicesList, allRanges);
        this.templates = new Substructure(templateIndicesList, allRanges);
        this.skips = new Substructure(skipIndicesList, allRanges);
        this.nonSkips = new Substructure(nonSkipIndicesList, allRanges);
    }

    /**
     * Converts readStructureString into a List<ReadDescriptor> and calls the primary constructor using this List as it's argument.
     *
     * @param readStructureString A string of the format <number of bases><type><number of bases><type>...<number of bases><type> describing
     *                            this read structure
     */
    public ReadStructure(final String readStructureString) {
        this(readStructureStringToDescriptors(readStructureString));
    }

    public int getNumDescriptors() {
        return descriptors.size();
    }

    /**
     * Converts this object into a String using rules complementary to the single string constructor above.
     *
     * @return A string of the form <number of bases><type><number of bases><type>...<number of bases><type> with one
     * <number of bases><type> per ReadDescriptor in descriptors.
     */
    @Override
    public String toString() {
        String out = "";
        for (final ReadDescriptor rd : descriptors) {
            out += rd.toString();
        }
        return out;
    }

    /**
     * Converts readStructureString into a List<ReadDescriptor>
     *
     * @param readStructure A string of the format <number of bases><type><number of bases><type>...<number of bases><type> describing
     *                      a read structure
     * @return A List<ReadDescriptor> corresponding to the input string
     */
    private final static List<ReadDescriptor> readStructureStringToDescriptors(final String readStructure) {
        final Matcher fullMatcher = FullPattern.matcher(readStructure);
        if (!fullMatcher.matches()) {
            throw new IllegalArgumentException(readStructure + " cannot be parsed as a ReadStructure! " + ReadStructureMsg);
        }


        final Matcher subMatcher = SubPattern.matcher(readStructure);
        final List<ReadDescriptor> descriptors = new ArrayList<ReadDescriptor>();
        while (subMatcher.find()) {
            final ReadDescriptor rd = new ReadDescriptor(parseInt(subMatcher.group(1)), valueOf(subMatcher.group(2)));
            descriptors.add(rd);
        }

        return descriptors;
    }

    @Override
    public boolean equals(final Object thatObj) {
        if (this == thatObj) return true;
        if (this.getClass() != thatObj.getClass()) return false;

        final ReadStructure that = (ReadStructure) thatObj;
        if (this.descriptors.size() != that.descriptors.size()) {
            return false;
        }

        for (int i = 0; i < this.descriptors.size(); i++) {
            if (!this.descriptors.get(i).equals(that.descriptors.get(i))) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int hashCode() {
        int res = descriptors.get(0).hashCode();
        for (int i = 1; i < descriptors.size(); i++) {
            res *= descriptors.get(i).hashCode();
        }

        return res;
    }

    /**
     * Represents a subset of ReadDescriptors in the containing ReadStructure, they ARE NOT necessarily contiguous
     * in the containing ReadStrucure but they ARE in the order they appear in the containing ReadStructure
     */
    public class Substructure implements Iterable<ReadDescriptor> {
        /**
         * Total number of descriptors == readTypeIndices.length
         */
        private final int numDescriptors;

        /**
         * The indices into the ReadStructure for this Substructure
         */
        private final int[] descriptorIndices;

        /**
         * The length of each individual ReadDescriptor in this substructure
         */
        private final int[] descriptorLengths;

        /**
         * Ranges of cycle indexes (cycle # - 1) covered by each descriptor
         */
        private final Range[] cycleIndexRanges;

        /**
         * The total number of cycles covered by this Substructure
         */
        private final int totalCycles;

        /**
         * Indices into the ReadStructure.descriptors for this specific substructure, indices
         * must be in the order they appear in the descriptors list (but the indices do NOT have to be continuous)
         *
         * @param descriptorIndices A list of indices into ReadStructure.descriptors of the enclosing ReadStructure
         * @param allRanges         A list of ranges for all reads (not just those in this substructure) in the same order as ReadStructure.descriptors
         */
        public Substructure(final List<Integer> descriptorIndices, final List<Range> allRanges) {
            this.numDescriptors = descriptorIndices.size();

            this.descriptorIndices = new int[numDescriptors];
            this.descriptorLengths = new int[numDescriptors];
            for (int i = 0; i < descriptorIndices.size(); i++) {
                this.descriptorIndices[i] = descriptorIndices.get(i);
                this.descriptorLengths[i] = descriptors.get(this.descriptorIndices[i]).length;
            }

            this.cycleIndexRanges = new Range[numDescriptors];
            for (int i = 0; i < numDescriptors; i++) {
                this.cycleIndexRanges[i] = allRanges.get(this.descriptorIndices[i]);
            }

            int totalLength = 0;
            for (final int length : descriptorLengths) {
                totalLength += length;
            }
            totalCycles = totalLength;
        }

        public ReadDescriptor get(final int index) {
            return descriptors.get(descriptorIndices[index]);
        }

        public boolean isEmpty() {
            return numDescriptors == 0;
        }

        public int length() {
            return numDescriptors;
        }

        public int getTotalCycles() {
            return totalCycles;
        }

        public int[] getIndices() {
            return descriptorIndices;
        }

        public int[] getDescriptorLengths() {
            return descriptorLengths;
        }

        public Range[] getCycleIndexRanges() {
            return cycleIndexRanges;
        }

        public Iterator<ReadDescriptor> iterator() {
            return new IndexedIterator(descriptorIndices);
        }

        public int[] getCycles() {
            int[] cycles = new int[totalCycles];
            int cycleIndex = 0;
            for (final Range range : cycleIndexRanges) {
                for (int i = range.start; i <= range.end; i++) {
                    cycles[cycleIndex++] = i + 1;
                }
            }
            return cycles;
        }

        /**
         * Create a ReadStructure from this substructure composed of only the descriptors contained in this substructure, Any
         * ReadDescriptors not in this substructure are treated as if they don't exist (e.g. if you have a readStructure
         * (36T8S8B36T) and this substructure consists of all the non-skipped reads than toReadStructure would return
         * (36T8B36T) in ReadStructure form
         */
        public ReadStructure toReadStructure() {
            final List<ReadDescriptor> descriptors = new ArrayList<ReadDescriptor>(numDescriptors);
            for (final ReadDescriptor rd : this) {
                descriptors.add(rd);
            }
            return new ReadStructure(descriptors);
        }
    }

    /**
     * An iterator over a Substructure's ReadDescriptors
     */
    private class IndexedIterator implements Iterator<ReadDescriptor> {
        private int index;
        private int[] indices;

        public IndexedIterator(final int[] indices) {
            this.indices = indices;
            this.index = 0;
        }

        public boolean hasNext() {
            return index < indices.length;
        }

        public ReadDescriptor next() {
            return descriptors.get(indices[index++]);
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
