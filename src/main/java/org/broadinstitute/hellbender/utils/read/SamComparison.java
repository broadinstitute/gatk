package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

// TODO: this class should be removed and replaced with the Picard version once
// https://github.com/broadinstitute/picard/pull/913 is merged into Picard and picked up by GATK

/**
 * Rudimentary SAM comparer. Compares headers, and if headers are compatible enough, compares SAMRecords,
 * looking only at basic alignment info. Summarizes the number of alignments that match, mismatch, are missing, etc.
 * Supplementary and secondary alignments are skipped. The inputs are assumed to be valid (as per SamFileValidator),
 * otherwise undefined behavior may occur.
 */
public final class SamComparison {
    private final SamReader reader1;
    private final SamReader reader2;

    private boolean sequenceDictionariesDiffer;
    private final boolean headersAreEqual;
    private final boolean alignmentsAreEqual;

    private int mappingsMatch = 0;
    private int unmappedBoth = 0;
    private int unmappedLeft = 0;
    private int unmappedRight = 0;
    private int mappingsDiffer = 0;
    private int missingLeft = 0;
    private int missingRight = 0;

    private enum AlignmentComparison {
        UNMAPPED_BOTH, UNMAPPED_LEFT, UNMAPPED_RIGHT, MAPPINGS_DIFFER, MAPPINGS_MATCH
    }

    /**
     * Note: the caller must make sure the SamReaders are closed properly.
     */
    public SamComparison(final SamReader reader1, final SamReader reader2) {
        this.reader1 = reader1;
        this.reader2 = reader2;
        this.headersAreEqual = compareHeaders();
        this.alignmentsAreEqual = compareAlignments();
    }

    public void printReport() {
        System.out.println("Match\t" + mappingsMatch);
        System.out.println("Differ\t" + mappingsDiffer);
        System.out.println("Unmapped_both\t" + unmappedBoth);
        System.out.println("Unmapped_left\t" + unmappedLeft);
        System.out.println("Unmapped_right\t" + unmappedRight);
        System.out.println("Missing_left\t" + missingLeft);
        System.out.println("Missing_right\t" + missingRight);
    }

    private boolean compareAlignments() {
        if (!compareValues(reader1.getFileHeader().getSortOrder(), reader2.getFileHeader().getSortOrder(), "Sort Order")) {
            System.out.println("Cannot compare alignments if sort orders differ.");
            return false;
        }
        switch (reader1.getFileHeader().getSortOrder()) {
            case coordinate:
                if (sequenceDictionariesDiffer) {
                    System.out.println("Cannot compare coordinate-sorted SAM files because sequence dictionaries differ.");
                    return false;
                }
                return compareCoordinateSortedAlignments();
            case queryname:
                return compareQueryNameSortedAlignments();
            case unsorted:
                return compareUnsortedAlignments();
            default:
                // unreachable
                return false;
        }
    }

    private boolean compareCoordinateSortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator itLeft =
                new SecondaryOrSupplementarySkippingIterator(reader1.iterator());
        final SecondaryOrSupplementarySkippingIterator itRight =
                new SecondaryOrSupplementarySkippingIterator(reader2.iterator());

        // Save any reads which haven't been matched during in-order scan.
        final Map<PrimaryAlignmentKey, SAMRecord> leftUnmatched = new LinkedHashMap<>();
        final Map<PrimaryAlignmentKey, SAMRecord> rightUnmatched = new LinkedHashMap<>();

        while (itLeft.hasCurrent()) {
            if (!itRight.hasCurrent()) {
                // Exhausted right side.  See if any of the remaining left reads match
                // any of the saved right reads.
                for (; itLeft.hasCurrent(); itLeft.advance()) {
                    final SAMRecord left = itLeft.getCurrent();
                    final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(left);
                    final SAMRecord right = rightUnmatched.remove(leftKey);
                    if (right == null) {
                        ++missingRight;
                    } else {
                        tallyAlignmentRecords(left, right);
                    }
                }
                break;
            }
            // Don't assume stability of order beyond the coordinate.  Therefore grab all the
            // reads from the left that has the same coordinate.
            final SAMRecord left = itLeft.getCurrent();
            final Map<PrimaryAlignmentKey, SAMRecord> leftCurrentCoordinate = new LinkedHashMap<>();
            final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(left);
            leftCurrentCoordinate.put(leftKey, left);
            while (itLeft.advance()) {
                final SAMRecord nextLeft = itLeft.getCurrent();
                if (compareAlignmentCoordinates(left, nextLeft) == 0) {
                    final PrimaryAlignmentKey nextLeftKey = new PrimaryAlignmentKey(nextLeft);
                    leftCurrentCoordinate.put(nextLeftKey, nextLeft);
                } else {
                    break;
                }
            }
            // Advance the right iterator until it is >= the left reads that have just been grabbed
            while (itRight.hasCurrent() && compareAlignmentCoordinates(left, itRight.getCurrent()) > 0) {
                final SAMRecord right = itRight.getCurrent();
                final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(right);
                rightUnmatched.put(rightKey, right);
                itRight.advance();
            }
            // For each right read that has the same coordinate as the current left reads,
            // see if there is a matching left read.  If so, process and discard.  If not,
            // save the right read for later.
            for (; itRight.hasCurrent() && compareAlignmentCoordinates(left, itRight.getCurrent()) == 0; itRight.advance()) {
                final SAMRecord right = itRight.getCurrent();
                final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(right);
                final SAMRecord matchingLeft = leftCurrentCoordinate.remove(rightKey);
                if (matchingLeft != null) {
                    tallyAlignmentRecords(matchingLeft, right);
                } else {
                    rightUnmatched.put(rightKey, right);
                }
            }

            // Anything left in leftCurrentCoordinate has not been matched
            for (final SAMRecord samRecord : leftCurrentCoordinate.values()) {
                final PrimaryAlignmentKey recordKey = new PrimaryAlignmentKey(samRecord);
                leftUnmatched.put(recordKey, samRecord);
            }
        }
        // The left iterator has been exhausted.  See if any of the remaining right reads
        // match any of the saved left reads.
        for (; itRight.hasCurrent(); itRight.advance()) {
            final SAMRecord right = itRight.getCurrent();
            final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(right);
            final SAMRecord left = leftUnmatched.remove(rightKey);
            if (left != null) {
                tallyAlignmentRecords(left, right);
            } else {
                ++missingLeft;
            }
        }

        // Look up reads that were unmatched from left, and see if they are in rightUnmatched.
        // If found, remove from rightUnmatched and tally.
        for (final Map.Entry<PrimaryAlignmentKey, SAMRecord> leftEntry : leftUnmatched.entrySet()) {
            final PrimaryAlignmentKey leftKey = leftEntry.getKey();
            final SAMRecord left = leftEntry.getValue();
            final SAMRecord right = rightUnmatched.remove(leftKey);
            if (right == null) {
                ++missingRight;
                continue;
            }
            tallyAlignmentRecords(left, right);
        }

        // Any elements remaining in rightUnmatched are guaranteed not to be in leftUnmatched.
        missingLeft += rightUnmatched.size();

        return allVisitedAlignmentsEqual();
    }

    private int compareAlignmentCoordinates(final SAMRecord left, final SAMRecord right) {
        final String leftReferenceName = left.getReferenceName();
        final String rightReferenceName = right.getReferenceName();
        if (leftReferenceName == null && rightReferenceName == null) {
            return 0;
        } else if (leftReferenceName == null) {
            return 1;
        } else if (rightReferenceName == null) {
            return -1;
        }
        final int leftReferenceIndex = reader1.getFileHeader().getSequenceIndex(leftReferenceName);
        final int rightReferenceIndex = reader1.getFileHeader().getSequenceIndex(rightReferenceName);

        if (leftReferenceIndex != rightReferenceIndex) {
            return leftReferenceIndex - rightReferenceIndex;
        }
        return left.getAlignmentStart() - right.getAlignmentStart();
    }

    private boolean compareQueryNameSortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(reader1.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(reader2.iterator());

        while (it1.hasCurrent()) {
            if (!it2.hasCurrent()) {
                missingRight += countRemaining(it1);
            }
            final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(it1.getCurrent());
            final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(it2.getCurrent());
            final int cmp = leftKey.compareTo(rightKey);
            if (cmp < 0) {
                ++missingRight;
                it1.advance();
            } else if (cmp > 0) {
                ++missingLeft;
                it2.advance();
            } else {
                tallyAlignmentRecords(it1.getCurrent(), it2.getCurrent());
                it1.advance();
                it2.advance();
            }
        }
        if (it2.hasCurrent()) {
            missingLeft += countRemaining(it2);
        }
        return allVisitedAlignmentsEqual();
    }

    /**
     * For unsorted alignments, assume nothing about the order. Determine which records to compare solely on the
     * basis of their PrimaryAlignmentKey.
     */
    private boolean compareUnsortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(reader1.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(reader2.iterator());

        final Map<PrimaryAlignmentKey, SAMRecord> leftUnmatched = new LinkedHashMap<>();
        for (; it1.hasCurrent(); it1.advance()) {
            final SAMRecord left = it1.getCurrent();
            final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(left);
            leftUnmatched.put(leftKey, left);
        }

        for (; it2.hasCurrent(); it2.advance()) {
            final SAMRecord right = it2.getCurrent();
            final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(right);
            final SAMRecord left = leftUnmatched.remove(rightKey);
            if (left != null) {
                tallyAlignmentRecords(left, right);
            } else {
                ++missingLeft;
            }
        }

        missingRight += leftUnmatched.size();

        return allVisitedAlignmentsEqual();
    }

    /**
     * Check the alignments tallied thus far for any kind of disparity.
     */
    private boolean allVisitedAlignmentsEqual() {
        return !(missingLeft > 0 || missingRight > 0 || mappingsDiffer > 0 || unmappedLeft > 0 || unmappedRight > 0);
    }

    private int countRemaining(final SecondaryOrSupplementarySkippingIterator it) {
        int i;
        for (i = 0; it.hasCurrent(); ++i) {
            it.advance();
        }
        return i;
    }

    private AlignmentComparison compareAlignmentRecords(final SAMRecord s1, final SAMRecord s2) {
        if (s1.getReadUnmappedFlag() && s2.getReadUnmappedFlag()) {
            return AlignmentComparison.UNMAPPED_BOTH;
        } else if (s1.getReadUnmappedFlag()) {
            return AlignmentComparison.UNMAPPED_LEFT;
        } else if (s2.getReadUnmappedFlag()) {
            return AlignmentComparison.UNMAPPED_RIGHT;
        } else if (alignmentsMatch(s1, s2)) {
            return AlignmentComparison.MAPPINGS_MATCH;
        } else {
            return AlignmentComparison.MAPPINGS_DIFFER;
        }
    }

    private boolean alignmentsMatch(final SAMRecord s1, final SAMRecord s2) {
        return (s1.getReferenceName().equals(s2.getReferenceName()) &&
                s1.getAlignmentStart() == s2.getAlignmentStart() &&
                s1.getReadNegativeStrandFlag() == s1.getReadNegativeStrandFlag());
    }

    /**
     * Compare the mapping information for two SAMRecords.
     */
    private void tallyAlignmentRecords(final SAMRecord s1, final SAMRecord s2) {
        if (!s1.getReadName().equals(s2.getReadName())) {
            throw new GATKException("Read names do not match: " + s1.getReadName() + " : " + s2.getReadName());
        }
        final AlignmentComparison comp = compareAlignmentRecords(s1, s2);
        switch (comp) {
            case UNMAPPED_BOTH:
                ++unmappedBoth;
                break;
            case UNMAPPED_LEFT:
                ++unmappedLeft;
                break;
            case UNMAPPED_RIGHT:
                ++unmappedRight;
                break;
            case MAPPINGS_DIFFER:
                ++mappingsDiffer;
                break;
            case MAPPINGS_MATCH:
                ++mappingsMatch;
                break;
            default:
                // unreachable
                throw new GATKException.ShouldNeverReachHereException("Unhandled comparison type: " + comp);
        }
    }

    private boolean compareHeaders() {
        final SAMFileHeader h1 = reader1.getFileHeader();
        final SAMFileHeader h2 = reader2.getFileHeader();
        boolean ret = compareValues(h1.getVersion(), h2.getVersion(), "File format version");
        ret = compareValues(h1.getCreator(), h2.getCreator(), "File creator") && ret;
        ret = compareValues(h1.getAttribute("SO"), h2.getAttribute("SO"), "Sort order") && ret;
        if (!compareSequenceDictionaries(h1, h2)) {
            ret = false;
            sequenceDictionariesDiffer = true;
        }
        ret = compareReadGroups(h1, h2) && ret;
        ret = compareProgramRecords(h1, h2) && ret;
        return ret;
    }

    private boolean compareProgramRecords(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMProgramRecord> l1 = h1.getProgramRecords();
        final List<SAMProgramRecord> l2 = h2.getProgramRecords();
        if (!compareValues(l1.size(), l2.size(), "Number of program records")) {
            return false;
        }
        boolean ret = true;
        for (int i = 0; i < l1.size(); ++i) {
            ret = compareProgramRecord(l1.get(i), l2.get(i)) && ret;
        }
        return ret;
    }

    private boolean compareProgramRecord(final SAMProgramRecord programRecord1, final SAMProgramRecord programRecord2) {
        if (programRecord1 == null && programRecord2 == null) {
            return true;
        }
        if (programRecord1 == null) {
            reportDifference("null", programRecord2.getProgramGroupId(), "Program Record");
            return false;
        }
        if (programRecord2 == null) {
            reportDifference(programRecord1.getProgramGroupId(), "null", "Program Record");
            return false;
        }
        boolean ret = compareValues(programRecord1.getProgramGroupId(), programRecord2.getProgramGroupId(),
                "Program Name");
        final String[] attributes = {"VN", "CL"};
        for (final String attribute : attributes) {
            ret = compareValues(programRecord1.getAttribute(attribute), programRecord2.getAttribute(attribute),
                    attribute + " Program Record attribute") && ret;
        }
        return ret;
    }

    private boolean compareReadGroups(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMReadGroupRecord> l1 = h1.getReadGroups();
        final List<SAMReadGroupRecord> l2 = h2.getReadGroups();
        if (!compareValues(l1.size(), l2.size(), "Number of read groups")) {
            return false;
        }
        boolean ret = true;
        for (int i = 0; i < l1.size(); ++i) {
            ret = compareReadGroup(l1.get(i), l2.get(i)) && ret;
        }
        return ret;
    }

    private boolean compareReadGroup(final SAMReadGroupRecord samReadGroupRecord1, final SAMReadGroupRecord samReadGroupRecord2) {
        boolean ret = compareValues(samReadGroupRecord1.getReadGroupId(), samReadGroupRecord2.getReadGroupId(),
                "Read Group ID");
        ret = compareValues(samReadGroupRecord1.getSample(), samReadGroupRecord2.getSample(),
                "Sample for read group " + samReadGroupRecord1.getReadGroupId()) && ret;
        ret = compareValues(samReadGroupRecord1.getLibrary(), samReadGroupRecord2.getLibrary(),
                "Library for read group " + samReadGroupRecord1.getReadGroupId()) && ret;
        final String[] attributes = {"DS", "PU", "PI", "CN", "DT", "PL"};
        for (final String attribute : attributes) {
            ret = compareValues(samReadGroupRecord1.getAttribute(attribute), samReadGroupRecord2.getAttribute(attribute),
                    attribute + " for read group " + samReadGroupRecord1.getReadGroupId()) && ret;
        }
        return ret;
    }

    private boolean compareSequenceDictionaries(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMSequenceRecord> s1 = h1.getSequenceDictionary().getSequences();
        final List<SAMSequenceRecord> s2 = h2.getSequenceDictionary().getSequences();
        if (s1.size() != s2.size()) {
            reportDifference(s1.size(), s2.size(), "Length of sequence dictionaries");
            return false;
        }
        boolean ret = true;
        for (int i = 0; i < s1.size(); ++i) {
            ret = compareSequenceRecord(s1.get(i), s2.get(i), i + 1) && ret;
        }
        return ret;
    }

    private boolean compareSequenceRecord(final SAMSequenceRecord sequenceRecord1, final SAMSequenceRecord sequenceRecord2, final int which) {
        if (!sequenceRecord1.getSequenceName().equals(sequenceRecord2.getSequenceName())) {
            reportDifference(sequenceRecord1.getSequenceName(), sequenceRecord2.getSequenceName(),
                    "Name of sequence record " + which);
            return false;
        }
        boolean ret = compareValues(sequenceRecord1.getSequenceLength(), sequenceRecord2.getSequenceLength(), "Length of sequence " +
                sequenceRecord1.getSequenceName());
        ret = compareValues(sequenceRecord1.getSpecies(), sequenceRecord2.getSpecies(), "Species of sequence " +
                sequenceRecord1.getSequenceName()) && ret;
        ret = compareValues(sequenceRecord1.getAssembly(), sequenceRecord2.getAssembly(), "Assembly of sequence " +
                sequenceRecord1.getSequenceName()) && ret;
        ret = compareValues(sequenceRecord1.getAttribute("M5"), sequenceRecord2.getAttribute("M5"), "MD5 of sequence " +
                sequenceRecord1.getSequenceName()) && ret;
        ret = compareValues(sequenceRecord1.getAttribute("UR"), sequenceRecord2.getAttribute("UR"), "URI of sequence " +
                sequenceRecord1.getSequenceName()) && ret;
        return ret;
    }

    private <T> boolean compareValues(final T v1, final T v2, final String label) {
        boolean eq = Objects.equals(v1, v2);
        if (eq) {
            return true;
        } else {
            reportDifference(v1, v2, label);
            return false;
        }
    }

    private void reportDifference(final String s1, final String s2, final String label) {
        System.out.println(label + " differs.");
        System.out.println("File 1: " + s1);
        System.out.println("File 2: " + s2);
    }

    private void reportDifference(Object o1, Object o2, final String label) {
        if (o1 == null) {
            o1 = "null";
        }
        if (o2 == null) {
            o2 = "null";
        }
        reportDifference(o1.toString(), o2.toString(), label);
    }

    public int getMappingsMatch() {
        return mappingsMatch;
    }

    public int getUnmappedBoth() {
        return unmappedBoth;
    }

    public int getUnmappedLeft() {
        return unmappedLeft;
    }

    public int getUnmappedRight() {
        return unmappedRight;
    }

    public int getMappingsDiffer() {
        return mappingsDiffer;
    }

    public int getMissingLeft() {
        return missingLeft;
    }

    public int getMissingRight() {
        return missingRight;
    }

    /**
     * Returns true if the headers are equal, false otherwise.
     */
    public boolean areHeadersEqual(){
        return headersAreEqual;
    }

    /**
     * Returns true if the alignments are equal, false otherwise.
     */
    public boolean areAlignmentsEqual(){
        return alignmentsAreEqual;
    }

    public boolean areEqual() {
        return headersAreEqual && alignmentsAreEqual;
    }
}
