/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Rudimentary SAM comparer.  Compares headers, and if headers are compatible enough, compares SAMRecords,
 * looking only at basic alignment info.  Summarizes the number of alignments that match, mismatch, are missing, etc.
 */
@CommandLineProgramProperties(
        usage = "USAGE: CompareSAMs <SAMFile1> <SAMFile2>\n" +
                "Compares the headers of the two input SAM or BAM files, and, if possible, the SAMRecords. " +
                "For SAMRecords, compares only the readUnmapped flag, reference name, start position and strand. " +
                "Reports the number of SAMRecords that match, differ in alignment, are mapped in only one input, " +
                "or are missing in one of the files",
        usageShort = "Compares two input SAM or BAM files",
        programGroup = ReadProgramGroup.class
)
public class CompareSAMs extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.", common = true, optional = true)
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    private final SamReader[] samReaders = new SamReader[2];
    private boolean sequenceDictionariesDiffer;
    private int mappingsMatch = 0;
    private int unmappedBoth = 0;
    private int unmappedLeft = 0;
    private int unmappedRight = 0;
    private int mappingsDiffer = 0;
    private int missingLeft = 0;
    private int missingRight = 0;
    private boolean areEqual;

    @Override
    protected Object doWork() {
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        for (int i = 0; i < samFiles.size(); ++i) {
            samReaders[i] = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(i));
        }
        areEqual = compareHeaders();
        areEqual = compareAlignments() && areEqual;
        printReport();
        if (areEqual) {
            System.out.println("SAM files match.");
        } else {
            System.out.println("SAM files differ.");
        }
        CloserUtil.close(samReaders);
        return null;
    }

    private void printReport() {
        System.out.println("Match\t" + mappingsMatch);
        System.out.println("Differ\t" + mappingsDiffer);
        System.out.println("Unmapped_both\t" + unmappedBoth);
        System.out.println("Unmapped_left\t" + unmappedLeft);
        System.out.println("Unmapped_right\t" + unmappedRight);
        System.out.println("Missing_left\t" + missingLeft);
        System.out.println("Missing_right\t" + missingRight);
    }

    private boolean compareAlignments() {
        if (!compareValues(samReaders[0].getFileHeader().getSortOrder(), samReaders[1].getFileHeader().getSortOrder(),
                "Sort Order")) {
            System.out.println("Cannot compare alignments if sort orders differ.");
            return false;
        }
        switch (samReaders[0].getFileHeader().getSortOrder()) {
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
                new SecondaryOrSupplementarySkippingIterator(samReaders[0].iterator());
        final SecondaryOrSupplementarySkippingIterator itRight =
                new SecondaryOrSupplementarySkippingIterator(samReaders[1].iterator());

        // Save any reads which haven't been matched during in-order scan.
        final Map<String, SAMRecord> leftUnmatched = new HashMap<String, SAMRecord>();
        final Map<String, SAMRecord> rightUnmatched = new HashMap<String, SAMRecord>();

        boolean ret = true;

        while (itLeft.hasCurrent()) {
            if (!itRight.hasCurrent()) {
                // Exhausted right side.  See if any of the remaining left reads match
                // any of the saved right reads.
                for (; itLeft.hasCurrent(); itLeft.advance()) {
                    final SAMRecord left = itLeft.getCurrent();
                    final SAMRecord right = rightUnmatched.remove(left.getReadName());
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
            final Map<String, SAMRecord> leftCurrentCoordinate = new HashMap<String, SAMRecord>();
            leftCurrentCoordinate.put(left.getReadName(), left);
            while (itLeft.advance()) {
                final SAMRecord nextLeft = itLeft.getCurrent();
                if (compareAlignmentCoordinates(left, nextLeft) == 0) {
                    leftCurrentCoordinate.put(nextLeft.getReadName(), nextLeft);
                } else {
                    break;
                }
            }
            // Advance the right iterator until it is >= the left reads that have just been grabbed
            while (itRight.hasCurrent() && compareAlignmentCoordinates(left, itRight.getCurrent()) > 0) {
                final SAMRecord right = itRight.getCurrent();
                rightUnmatched.put(right.getReadName(), right);
                itRight.advance();
            }
            // For each right read that has the same coordinate as the current left reads,
            // see if there is a matching left read.  If so, process and discard.  If not,
            // save the right read for later.
            for (; itRight.hasCurrent() && compareAlignmentCoordinates(left, itRight.getCurrent()) == 0; itRight.advance()) {
                final SAMRecord right = itRight.getCurrent();
                final SAMRecord matchingLeft = leftCurrentCoordinate.remove(right.getReadName());
                if (matchingLeft != null) {
                    ret = tallyAlignmentRecords(matchingLeft, right) && ret;
                } else {
                    rightUnmatched.put(right.getReadName(), right);
                }
            }

            // Anything left in leftCurrentCoordinate has not been matched
            for (final SAMRecord samRecord : leftCurrentCoordinate.values()) {
                leftUnmatched.put(samRecord.getReadName(), samRecord);
            }
        }
        // The left iterator has been exhausted.  See if any of the remaining right reads
        // match any of the saved left reads.
        for (; itRight.hasCurrent(); itRight.advance()) {
            final SAMRecord right = itRight.getCurrent();
            final SAMRecord left = leftUnmatched.remove(right.getReadName());
            if (left != null) {
                tallyAlignmentRecords(left, right);
            } else {
                ++missingLeft;
            }
        }

        // Look up reads that were unmatched from left, and see if they are in rightUnmatched.
        // If found, remove from rightUnmatched and tally.
        for (final Map.Entry<String, SAMRecord> leftEntry : leftUnmatched.entrySet()) {
            final String readName = leftEntry.getKey();
            final SAMRecord left = leftEntry.getValue();
            final SAMRecord right = rightUnmatched.remove(readName);
            if (right == null) {
                ++missingRight;
                continue;
            }
            tallyAlignmentRecords(left, right);
        }

        // Any elements remaining in rightUnmatched are guaranteed not to be in leftUnmatched.
        missingLeft += rightUnmatched.size();

        if (ret) {
            if (missingLeft > 0 || missingRight > 0 || mappingsDiffer > 0 || unmappedLeft > 0 || unmappedRight > 0) {
                ret = false;
            }
        }
        return ret;
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
        final int leftReferenceIndex = samReaders[0].getFileHeader().getSequenceIndex(leftReferenceName);
        final int rightReferenceIndex = samReaders[0].getFileHeader().getSequenceIndex(rightReferenceName);

        if (leftReferenceIndex != rightReferenceIndex) {
            return leftReferenceIndex - rightReferenceIndex;
        }
        return left.getAlignmentStart() - right.getAlignmentStart();
    }

    private boolean compareQueryNameSortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(samReaders[0].iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(samReaders[1].iterator());

        boolean ret = true;
        while (it1.hasCurrent()) {
            if (!it2.hasCurrent()) {
                missingRight += countRemaining(it1);
                return false;
            }
            final int cmp = it1.getCurrent().getReadName().compareTo(it2.getCurrent().getReadName());
            if (cmp < 0) {
                ++missingRight;
                it1.advance();
                ret = false;
            } else if (cmp > 0) {
                ++missingLeft;
                it2.advance();
                ret = false;
            } else {
                if (!tallyAlignmentRecords(it1.getCurrent(), it2.getCurrent())) {
                    ret = false;
                }
                it1.advance();
                it2.advance();
            }
        }
        if (it2.hasCurrent()) {
            missingLeft += countRemaining(it2);
            return false;
        }
        return ret;
    }

    private boolean compareUnsortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(samReaders[0].iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(samReaders[1].iterator());
        boolean ret = true;
        for (; it1.hasCurrent(); it1.advance(), it2.advance()) {
            if (!it2.hasCurrent()) {
                missingRight += countRemaining(it1);
                return false;
            }
            final SAMRecord s1 = it1.getCurrent();
            final SAMRecord s2 = it2.getCurrent();
            if (!compareValues(s1.getReadName(), s2.getReadName(), "Read names")) {
                System.out.println("Read names cease agreeing in unsorted SAM files .  Comparison aborting.");
            }
            ret = tallyAlignmentRecords(s1, s2) && ret;
        }

        if (it2.hasCurrent()) {
            missingLeft += countRemaining(it2);
            return false;
        }
        return ret;
    }

    private int countRemaining(final SecondaryOrSupplementarySkippingIterator it) {
        int i;
        for (i = 0; it.hasCurrent(); ++i) {
            it.advance();
        }
        return i;
    }

    private boolean tallyAlignmentRecords(final SAMRecord s1, final SAMRecord s2) {
        if (!s1.getReadName().equals(s2.getReadName())) {
            throw new GATKException("Read names do not match: " + s1.getReadName() + " : " + s2.getReadName());
        }
        if (s1.getReadUnmappedFlag() && s2.getReadUnmappedFlag()) {
            ++unmappedBoth;
            return true;
        }
        if (s1.getReadUnmappedFlag()) {
            ++unmappedLeft;
            return false;
        }
        if (s2.getReadUnmappedFlag()) {
            ++unmappedRight;
            return false;
        }
        final boolean ret = (s1.getReferenceName().equals(s2.getReferenceName()) &&
                s1.getAlignmentStart() == s2.getAlignmentStart() &&
                s1.getReadNegativeStrandFlag() == s1.getReadNegativeStrandFlag());
        if (!ret) {
            ++mappingsDiffer;
        } else {
            ++mappingsMatch;
        }
        return ret;
    }

    private boolean compareHeaders() {
        final SAMFileHeader h1 = samReaders[0].getFileHeader();
        final SAMFileHeader h2 = samReaders[1].getFileHeader();
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
        System.out.println(samFiles.get(0) + ": " + s1);
        System.out.println(samFiles.get(1) + ": " + s2);
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

    public boolean areEqual() {
        return areEqual;
    }

}
