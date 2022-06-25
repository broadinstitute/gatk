package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.Logger;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class CopyNumberTestUtils {

    /**
     * This test utility method asserts if two TSV files (each containing a header followed by a matrix of doubles)
     * are approximately equivalent, by checking if corresponding double values are equal up to some number of decimal
     * points represented by the delta parameter. Note that the headers of both files need to be bitwise identical,
     * and the number of columns and rows of each matrix encoded in TSV files must be equal as well.
     *
     * @param f1 first file to compare
     * @param f2 second file to compare
     * @param delta the absolute tolerable difference between the actual and the expected values
     * @param logger the Logger object to use
     */
    public static void assertFilesEqualUpToAllowedDeltaForDoubleValues(final File f1, final File f2, final double delta,
                                                                       final Logger logger) {
        try {
            IOUtil.assertFilesEqual(f1, f2);
        } catch (final SAMException e) {
            logger.warn(String.format("Files %s and %s failed exact-match test, attempting comparison of doubles at the %6.3e level...%n",
                    f1, f2, delta));
            try {
                final List<String> lines1 = Files.lines(f1.toPath()).collect(Collectors.toList());
                final List<String> lines2 = Files.lines(f2.toPath()).collect(Collectors.toList());
                Assert.assertEquals(lines1.size(), lines2.size(),
                        String.format("Files %s and %s do not have the same number of lines.", f1, f2));

                for (int i = 0; i < lines1.size(); i++) {
                    final String line1 = lines1.get(i);
                    final String line2 = lines2.get(i);
                    if (line1.equals(line2)) {
                        continue;
                    }

                    final List<String> splitLine1 = Arrays.asList(line1.split("\t"));
                    final List<String> splitLine2 = Arrays.asList(line2.split("\t"));
                    Assert.assertEquals(splitLine1.size(), splitLine2.size(),
                            String.format("Line %d does not have the same number of fields in files %s and %s.", i, f1, f2));

                    for (int j = 0; j < splitLine1.size(); j++) {
                        final String field1 = splitLine1.get(j);
                        final String field2 = splitLine2.get(j);
                        if (field1.equals(field2)) {
                            continue;
                        }
                        try {
                            final double double1 = Double.parseDouble(field1);
                            final double double2 = Double.parseDouble(field2);
                            Assert.assertEquals(double1, double2, delta,
                                    String.format("Field %d in line %d in files %s and %s is not equivalent at the %6.3e level: %f != %f.",
                                            j, i, f1, f2, delta, double1, double2));
                        } catch (final NumberFormatException nfe) {
                            Assert.fail(String.format("Non-double field %d in line %d in files %s and %s is not identical: %s != %s.",
                                    j, i, f1, f2, field1, field2));
                        }
                    }
                }
            } catch (final IOException ioe) {
                Assert.fail(String.format("Encountered IOException when trying to compare %s and %s: %s", f1, f2, ioe));
            }
        }
    }
}
