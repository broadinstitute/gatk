package org.broadinstitute.hellbender.utils.test;

import com.google.common.collect.Sets;
import htsjdk.samtools.*;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.tools.picard.sam.SortSam;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.SamComparison;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Collection of utilities for making common assertions about SAM files for unit testing purposes.
 */
public final class SamAssertionUtils {

    private static SamReader getReader(final File sam, final ValidationStringency validationStringency, final File reference) {
        return SamReaderFactory.makeDefault().validationStringency(validationStringency).referenceSequence(reference).open(sam);
    }

    /**
     *  causes an exception if the given sam files aren't equal
     *  @param actualSam the actual file
     *  @param expectedSam the expected file
     *  @param validationStringency how stringently do we validate the files
     *  @param reference is allowed to be null
     */
    public static void assertSamsEqual(final File actualSam, final File expectedSam, final ValidationStringency validationStringency, final File reference) throws IOException {
        final String equalStringent = samsEqualStringent(actualSam, expectedSam, validationStringency, reference);
        Assert.assertNull(equalStringent, "SAM file " + actualSam.getPath() + " differs from expected output:" + expectedSam.getPath() + " " + equalStringent);
    }

    /**
     * causes an exception if the given sam files aren't equal
     *  @param actualSam the actual file
     *  @param expectedSam the expected file
     *  @param validationStringency how stringently do we validate the files
     */
    public static void assertSamsEqual(final File actualSam, final File expectedSam, final ValidationStringency validationStringency) throws IOException {
        assertSamsEqual(actualSam, expectedSam, validationStringency, null);
    }

    /**
     * causes an exception if the given sam files aren't equal
     *  @param actualSam the actual file
     *  @param expectedSam the expected file
     *  @param reference is allowed to be null
     */
    public static void assertSamsEqual(final File actualSam, final File expectedSam, final File reference) throws IOException {
        assertSamsEqual(actualSam, expectedSam, ValidationStringency.DEFAULT_STRINGENCY, reference);
    }

    /**
     * causes an exception if the given sam files aren't equal
     *  @param actualSam the actual file
     *  @param expectedSam the expected file
     */
    public static void assertSamsEqual(final File actualSam, final File expectedSam) throws IOException {
        assertSamsEqual(actualSam, expectedSam, ValidationStringency.DEFAULT_STRINGENCY, null);
    }

    /**
     * causes an exception if the given sam isn't valid
     * @param reference is allowed to be null
     */
    public static void assertSamValid(final File sam, final ValidationStringency validationStringency, final File reference) throws IOException {
        assertCRAMContentsIfCRAM(sam);
        try (final SamReader samReader = getReader(sam, validationStringency, reference)) {
            final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
            validator.setIgnoreWarnings(true);
            validator.setVerbose(true, 1000);
            validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
            final boolean validated = validator.validateSamFileVerbose(samReader, null);
            Assert.assertTrue(validated, "SAM file validation failed");
        }
    }

    /**
     * causes an exception if the given sam isn't valid
     */
    public static void assertSamValid(final File sam, final ValidationStringency validationStringency) throws IOException {
        assertSamValid(sam, validationStringency, null);
    }

    /**
     * causes an exception if the given sam isn't valid
     * @param reference is allowed to be null
     * the default ValidationStringency value for this method is LENIENT
     */
    public static void assertSamValid(final File sam, final File reference) throws IOException {
        assertSamValid(sam, ValidationStringency.LENIENT, reference);
    }

    /**
     * causes an exception if the given sam isn't valid
     * the default ValidationStringency value for this method is LENIENT
     */
    public static void assertSamValid(final File sam) throws IOException {
        assertSamValid(sam, ValidationStringency.LENIENT, null);
    }

    /**
     * Compares SAM/BAM files in a stringent way but not by byte identity (allow reorder of attributes).
     * Returns null if the files are considered equals and returns a String describing the reason for comparison failure.
     * The lenient comparison only checks headers and alignment info {@see SamComparison}. Compares headers, and if headers are compatible enough, compares SAMRecords,
     * looking only at basic alignment info.
     */
    public static String samsEqualLenient(final File actualSam, final File expectedSam, final ValidationStringency validation, final File reference) throws IOException {
        assertCRAMContentsIfCRAM(actualSam);
        assertCRAMContentsIfCRAM(expectedSam);
        try(final SamReader reader1 = getReader(actualSam, validation, reference);
            final SamReader reader2 = getReader(expectedSam, validation, reference)) {

            final SamComparison comparison = new SamComparison(reader1, reader2);
            return comparison.areEqual() ? null : "SamComparison fails";
        }
    }

    /**
     * Compares SAM/BAM files in a stringent way but not by byte identity (allow reorder of attributes)
     * Comparing by MD5s is too strict and comparing by SamComparison is too lenient. So we need this method.
     *
     * This differs from a byte-to-byte comparison:
     * - @PG and @CO lines in headers are ignored in the comparison.
     * - each read in the actual file are allowed to have a superset of the attributes of the corresponding read in the expected set
     * @return null if equal or message string if not equal.
     */
    public static String samsEqualStringent(final File actualSam, final File expectedSam, final ValidationStringency validation, final File reference) throws IOException {
        if (sameMD5s(actualSam, expectedSam)) {
            return null;
        }

        //  verify that CRAM files have CRAM contents
        assertCRAMContentsIfCRAM(actualSam);
        assertCRAMContentsIfCRAM(expectedSam);

        String msg = equalHeadersIgnoreCOandPG(actualSam, expectedSam, validation, reference);
        if (msg != null) { return msg; }

        //At this point we know that the files are not byte-wise identical, but are equal according to SamComparison and their headers are equal
        //So we iterate over reads and compare them one by one.
        return compareReads(actualSam, expectedSam, validation, reference);
    }


    private static boolean sameMD5s(final File actualSam, final File expectedSam) throws IOException {
        final String fileMD5_1 = Utils.calculateFileMD5(actualSam);
        final String fileMD5_2 = Utils.calculateFileMD5(expectedSam);
        return fileMD5_1.equals(fileMD5_2);
    }

    private static String compareReads(final File actualSam, final File expectedSam, final ValidationStringency validation, final File reference) throws IOException {
        try(final SamReader reader1 = getReader(actualSam, validation, reference);
            final SamReader reader2 = getReader(expectedSam, validation, reference)) {
            final SAMRecordIterator it1 = reader1.iterator();
            final SAMRecordIterator it2 = reader2.iterator();
            while (it1.hasNext() && it2.hasNext()) {
                final SAMRecord read1 = it1.next();
                final SAMRecord read2 = it2.next();
                final String eqMessage = readsEqualAllowAddingAttributes(read1, read2);
                if (eqMessage != null){
                    return eqMessage;
                }
            }
            if (it1.hasNext() || it2.hasNext()) {
                //at least one has no more records (because the while loop is done) and at least one does have records. So we're not equal.
                return "Not the same number of reads";
            }
            return null;
        }
    }

    private static String equalHeadersIgnoreCOandPG(final File actualSam, final File expectedSam, final ValidationStringency validation, final File reference) throws IOException {
        try(final SamReader reader1 = getReader(actualSam, validation, reference);
            final SamReader reader2 = getReader(expectedSam, validation, reference)){

            final SAMFileHeader h1 = reader1.getFileHeader();
            final SAMFileHeader h2 = reader2.getFileHeader();
            String msg;

            //Note: we allow the versions to differ

            msg = compareValues(h1.getCreator(), h2.getCreator(), "File creator");
            if (msg != null) { return msg; }

            msg = compareValues(h1.getAttribute("SO"), h2.getAttribute("SO"), "Sort order");
            if (msg != null) { return msg; }

            if (! Objects.equals(h1.getSequenceDictionary(), h2.getSequenceDictionary())){
                return "Different Sequence dictionaries";
            }

            msg = compareReadGroups(h1, h2);
            if (msg != null) { return msg; }

            return msg;
        }
    }

    private static String compareReadGroups(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMReadGroupRecord> l1 = h1.getReadGroups();
        final List<SAMReadGroupRecord> l2 = h2.getReadGroups();
        final String msg = compareValues(l1.size(), l2.size(), "Number of read groups");
        if (msg != null){ return msg; }

        for (int i = 0; i < l1.size(); ++i) {
            if (! Objects.equals(l1.get(i), l2.get(i))){
                 return "Read group records different:" + l1.get(i) + " vs " + l2.get(i);
            }
        }
        return null;
    }

    private static <T> String compareValues(final T v1, final T v2, final String label) {
        boolean eq = Objects.equals(v1, v2);
        if (eq) {
            return null;
        } else {
            final String s1 = String.valueOf(v1);
            final String s2 = String.valueOf(v2);
            return label + " differs. File 1: " + s1 + " File 2: " + s2;
        }
    }

    /**
     * Compares the reads but ignores order of attributes.
     * Also allows actualRead to have a superset of attributes of expectedRead.
     */
    private static String readsEqualAllowAddingAttributes(final SAMRecord actualRead, final SAMRecord expectedRead) {
        final String actualName = actualRead.getReadName();
        final String expectedName = expectedRead.getReadName();

        String msg;

        msg = compareValues(actualName, expectedName, "name");
        if (msg != null){ return msg; }

        final String readNames = "actualName:" + actualName + " expectedName:" + expectedName;

        msg = compareValues(SAMFlag.getFlags(actualRead.getFlags()), SAMFlag.getFlags(expectedRead.getFlags()), readNames + " getFlags");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getInferredInsertSize(), expectedRead.getInferredInsertSize(), readNames + " getInferredInsertSize");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMappingQuality(), expectedRead.getMappingQuality(), readNames + " getMappingQuality");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMateReferenceIndex(), expectedRead.getMateReferenceIndex(), readNames + "getMateReferenceIndex");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMateAlignmentStart(), expectedRead.getMateAlignmentStart(), readNames + "getMateAlignmentStart");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getReferenceIndex(), expectedRead.getReferenceIndex(), readNames + " getReferenceIndex");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getAlignmentStart(), expectedRead.getAlignmentStart(), readNames + " getAlignmentStart");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getCigar(), expectedRead.getCigar(), readNames + " getCigar");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getReferenceName(), expectedRead.getReferenceName(), readNames + " getReferenceName");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMateReferenceName(), expectedRead.getMateReferenceName(), readNames + " getMateReferenceName");
        if (msg != null){ return msg; }

        if (!Arrays.equals(actualRead.getReadBases(), expectedRead.getReadBases())){
            return "getReadBases different actualRead:" + actualName + " expectedRead:" + expectedName + " (" + Arrays.toString(actualRead.getReadBases()) + " vs " + Arrays.toString(expectedRead.getReadBases()) + ")";
        }
        if (!Arrays.equals(actualRead.getBaseQualities(), expectedRead.getBaseQualities())){
            return "getBaseQualities different actualRead:" + actualName + " expectedRead:" + expectedName + " (" + Arrays.toString(actualRead.getBaseQualities()) + " vs " + Arrays.toString(expectedRead.getBaseQualities()) + ")";
        }
        return compareReadAttibutes(actualRead, expectedRead);
    }

    private static String compareReadAttibutes(final SAMRecord actualRead, final SAMRecord expectedRead) {
        final String actualName = actualRead.getReadName();
        final String expectedName = expectedRead.getReadName();
        final String readNames = "actualName:" + actualName + " expectedName:" + expectedName;

        String msg;
        final List<SAMRecord.SAMTagAndValue> actualAttributes = actualRead.getAttributes();
        final List<SAMRecord.SAMTagAndValue> expectedAttributes = expectedRead.getAttributes();

        //We want to compare attributes regardless of order, so we put them in a map
        final Map<String, Object> actualAttributesByName = new HashMap<>();
        final Map<String, Object> expectedAttributesByName = new HashMap<>();

        for (final SAMRecord.SAMTagAndValue samTagAndValue : actualAttributes) {
            actualAttributesByName.put(samTagAndValue.tag, samTagAndValue.value);
        }
        for (final SAMRecord.SAMTagAndValue samTagAndValue : expectedAttributes) {
            expectedAttributesByName.put(samTagAndValue.tag, samTagAndValue.value);
        }

        final Sets.SetView<String> attrDiff = Sets.difference(expectedAttributesByName.keySet(), actualAttributesByName.keySet());
        if (!attrDiff.isEmpty()){
            final StringBuilder sb= new StringBuilder();
            sb.append("expected read contains attributes that actual read lacks: " + readNames + " " + attrDiff + "\n");
            for (final String attr : attrDiff) {
                sb.append(attr + " " + expectedAttributesByName.get(attr) + "\n");
            }
            return sb.toString();
        }

        for (int i = 0; i < expectedAttributesByName.size(); i++) {
            final String expectedTag = expectedAttributes.get(i).tag;
            final Object expectedValue = expectedAttributesByName.get(expectedTag);
            final Object actualValue = actualAttributesByName.get(expectedTag);

            msg = compareValues(actualValue, expectedValue, readNames + " attribute " + expectedTag);
            if (msg != null){ return msg; }
        }
        return null;
    }

    /**
     * Compares the two given bam files, optionally sorting them before comparison.
     * The sorting is helpful to compare files that are different but equivalent (eg read pairs with same coordinates get reordered).
     */
    public static void assertEqualBamFiles(
            final File resultFile,
            final File expectedFile,
            final boolean compareBamFilesSorted,
            final ValidationStringency stringency) throws IOException {
        assertEqualBamFiles(resultFile, expectedFile, null, compareBamFilesSorted, stringency);
    }

    /**
     * Compares the two given bam files, optionally sorting them before comparison.
     * The sorting is helpful to compare files that are different but equivalent (eg read pairs with same coordinates get reordered).
     */
    public static void assertEqualBamFiles(
            final File resultFile,
            final File expectedFile,
            final File reference,
            final boolean compareBamFilesSorted,
            final ValidationStringency stringency) throws IOException {

        if (compareBamFilesSorted) {
            final File resultFileSorted= BaseTest.createTempFile("resultsFileSorted", "." + FilenameUtils.getExtension(resultFile.getName()));
            final File expectedFileSorted = BaseTest.createTempFile("expectedFileSorted", "." + FilenameUtils.getExtension(expectedFile.getName()));

            sortSam(resultFile, resultFileSorted, reference, stringency);
            sortSam(expectedFile, expectedFileSorted, reference, stringency);

            assertSamsEqual(resultFileSorted, expectedFileSorted, stringency, reference);
        } else {
            assertSamsEqual(resultFile, expectedFile, stringency, reference);
        }
    }

    /**
     * Validate/assert that the contents are CRAM if the extension is .cram
     */
    public static void assertCRAMContentsIfCRAM(final File putativeCRAMFile) {
        if (IOUtils.isCramFile(putativeCRAMFile)) {
            assertCRAMContents(putativeCRAMFile);
        }
    }

    /**
     * Unconditionally validate/assert that the contents are CRAM
     */
    public static void assertCRAMContents(final File putativeCRAMFile) {
        Assert.assertTrue(ReadUtils.hasCRAMFileContents(putativeCRAMFile));
    }

    private static void sortSam(final File input, final File output, final File reference, final ValidationStringency stringency) {
        final SortSam sort = new SortSam();
        sort.INPUT = input;
        sort.OUTPUT = output;
        sort.REFERENCE_SEQUENCE = reference;
        sort.SORT_ORDER = SAMFileHeader.SortOrder.coordinate;
        sort.VALIDATION_STRINGENCY = stringency;
        sort.runTool();
    }
}
