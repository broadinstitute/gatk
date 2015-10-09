package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
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
     *  @param reference is allowed to be null
     */
    public static void assertSamsEqual(final File sam1, final File sam2, final ValidationStringency validationStringency, final File reference) throws IOException {
        try (final SamReader reader1 = getReader(sam1, validationStringency, reference);
             final SamReader reader2 = getReader(sam2, validationStringency, reference)) {
            final SamComparison comparison = new SamComparison(reader1, reader2);
            final boolean equal = comparison.areEqual();
            Assert.assertTrue(equal, "SAM file " + sam1.getPath() + " differs from expected output:" + sam2.getPath());
        }
    }

    /**
     * causes an exception if the given sam files aren't equal
     */
    public static void assertSamsEqual(final File sam1, final File sam2, final ValidationStringency validationStringency) throws IOException {
        assertSamsEqual(sam1, sam2, validationStringency, null);
    }

    /**
     * causes an exception if the given sam files aren't equal
     * @param reference is allowed to be null
     * the default ValidationStringency value for this method is DEFAULT_STRINGENCY
     */
    public static void assertSamsEqual(final File sam1, final File sam2, final File reference) throws IOException {
        assertSamsEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY, reference);
    }

    /**
     * causes an exception if the given sam files aren't equal
     * the default ValidationStringency value for this method is DEFAULT_STRINGENCY
     */
    public static void assertSamsEqual(final File sam1, final File sam2) throws IOException {
        assertSamsEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY, null);
    }

    /**
     * causes an exception if the given sam files are equal
     * @param reference is allowed to be null
     */
    public static void assertSamsNonEqual(final File sam1, final File sam2, final ValidationStringency validationStringency, final File reference) throws IOException {
        try (final SamReader reader1 = getReader(sam1, validationStringency, reference);
             final SamReader reader2 = getReader(sam2, validationStringency, reference)) {
            final SamComparison comparison = new SamComparison(reader1, reader2);
            final boolean equal = comparison.areEqual();
            Assert.assertFalse(equal, "SAM files are expected to differ, but they do not");
        }
    }

    /**
     * causes an exception if the given sam files are equal
     */
    public static void assertSamsNonEqual(final File sam1, final File sam2, final ValidationStringency validationStringency) throws IOException {
        assertSamsNonEqual(sam1, sam2, validationStringency, null);
    }

    /**
     * causes an exception if the given sam files are equal
     * @param reference is allowed to be null
     * the default ValidationStringency value for this method is DEFAULT_STRINGENCY
     */
    public static void assertSamsNonEqual(final File sam1, final File sam2, final File reference) throws IOException {
        assertSamsNonEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY, reference);
    }

    /**
     * causes an exception if the given sam files are equal
     * the default ValidationStringency value for this method is DEFAULT_STRINGENCY
     */
    public static void assertSamsNonEqual(final File sam1, final File sam2) throws IOException {
        assertSamsNonEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY, null);
    }

    /**
     * causes an exception if the given sam isn't valid
     * @param reference is allowed to be null
     */
    public static void assertSamValid(final File sam, final ValidationStringency validationStringency, final File reference) throws IOException {
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

}
