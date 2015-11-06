package org.broadinstitute.hellbender.utils.io;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import org.broadinstitute.hellbender.utils.test.BaseTest;


/**
 * Tests selected functionality in the CommandLineExecutable class
 */
public class ListFileUtilsUnitTest extends BaseTest {


    @Test
    public void testUnpackSet() throws Exception {
        Set<String> expected = new HashSet<>(Arrays.asList(publicTestDir + "exampleBAM.bam"));
        Set<String> actual = ListFileUtils.unpackSet(Arrays.asList(publicTestDir + "exampleBAM.bam"));
        Assert.assertEquals(actual, expected);

        File tempListFile = createTempListFile("testUnpackSet",
                "#",
                publicTestDir + "exampleBAM.bam",
                "#" + publicTestDir + "foo.bam",
                "      # " + publicTestDir + "bar.bam"
        );
        actual = ListFileUtils.unpackSet(Arrays.asList(tempListFile.getAbsolutePath()));
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name="includeMatchingTests")
    public Object[][] getIncludeMatchingTests() {
        return new Object[][] {
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), true, new HashSet<>(Arrays.asList("a")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), false, new HashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), true, Collections.EMPTY_SET },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), false, new HashSet<>(Arrays.asList("ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), true, new HashSet<>(Arrays.asList("a")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), false, new HashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), true, new HashSet<>(Arrays.asList("a", "ab")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), false, new HashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), true, Collections.EMPTY_SET },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), false, new HashSet<>(Arrays.asList("ab", "abc") )},
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), true, Collections.EMPTY_SET },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), false, new HashSet<>(Arrays.asList("a", "ab", "abc") )}
        };
    }

    @Test(dataProvider = "includeMatchingTests")
    public void testIncludeMatching(Set<String> values, Collection<String> filters, boolean exactMatch, Set<String> expected) {
        Set<String> actual = ListFileUtils.includeMatching(values, ListFileUtils.IDENTITY_STRING_CONVERTER, filters, exactMatch);
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name="excludeMatchingTests")
    public Object[][] getExcludeMatchingTests() {
        return new Object[][] {
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), true, new HashSet<>(Arrays.asList("ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), false, Collections.EMPTY_SET },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), true, new HashSet<>(Arrays.asList("a", "ab", "abc") )},
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), false, new HashSet<>(Arrays.asList("a")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), true, new HashSet<>(Arrays.asList("ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), false, Collections.EMPTY_SET },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), true, new HashSet<>(Arrays.asList("abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), false, Collections.EMPTY_SET },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), true, new HashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), false, new HashSet<>(Arrays.asList("a")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), true, new HashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new HashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), false, Collections.EMPTY_SET }
        };
    }

    @Test(dataProvider = "excludeMatchingTests")
    public void testExcludeMatching(Set<String> values, Collection<String> filters, boolean exactMatch, Set<String> expected) {
        Set<String> actual = ListFileUtils.excludeMatching(values, ListFileUtils.IDENTITY_STRING_CONVERTER, filters, exactMatch);
        Assert.assertEquals(actual, expected);
    }

    private static File createTempListFile(final String tempFilePrefix, final String... lines) {
        try {
            final File tempListFile = createTempFile(tempFilePrefix, ".list");

            try (final PrintWriter out = new PrintWriter(tempListFile)) {
                for (final String line : lines) {
                    out.println(line);
                }
            }
            return tempListFile;
        } catch (IOException ex) {
            throw new UserException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

}
