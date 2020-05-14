package org.broadinstitute.hellbender.engine;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import org.apache.commons.lang3.SystemUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.*;

public class GATKPathSpecifierUnitTest extends GATKBaseTest {

    final static String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

    @DataProvider
    public Object[][] validPathSpecifiers() {
        return new Object[][] {
                // Paths specifiers that ar syntactically valid as either a file name or a URI and can be
                // represented internally as a URI, but which may fail hasFileSystemProvider or isPath

                // input String, expected resulting URI String, expected hasFileSystemProvider, expected isPath

                //********************************
                // Local (non-URI) file references
                //********************************

                {"localFile.bam",                   "file://" + getCWDAsURIPathString() + "localFile.bam", true, true},
                // absolute reference to a file in the root of the current file system (Windows accepts the "/" as root)
                {"/localFile.bam",                  "file://" + getRootDirectoryAsURIPathString() + "localFile.bam", true, true},
                // absolute reference to a file in the root of the current file system, where root is specified using the
                // default FS separator
                {FS_SEPARATOR + "localFile.bam",  "file://" + getRootDirectoryAsURIPathString() + "localFile.bam", true, true},
                // absolute reference to a file
                {FS_SEPARATOR + joinWithFSSeparator("path", "to", "localFile.bam"),
                        "file://" + getRootDirectoryAsURIPathString() + "path/to/localFile.bam",  true, true},
                // absolute reference to a file that contains a URI excluded character in the path ("#") which without
                // encoding will be treated as a fragment delimiter
                {FS_SEPARATOR + joinWithFSSeparator("project", "gvcf-pcr", "23232_1#1", "1.g.vcf.gz"),
                        "file://" + getRootDirectoryAsURIPathString() + "project/gvcf-pcr/23232_1%231/1.g.vcf.gz", true, true},
                // relative reference to a file on the local file system (relative to the current working directory)
                {joinWithFSSeparator("path", "to", "localFile.bam"),
                        "file://" + getCWDAsURIPathString() + "path/to/localFile.bam", true, true},
                // Windows also accepts "/" as a valid root specifier
                {"/", "file://" + getRootDirectoryAsURIPathString(), true, true},
                {".", "file://" + getCWDAsURIPathString() + "./", true, true},
                {"../.", "file://" + getCWDAsURIPathString() + ".././", true, true},
                // an empty path is equivalent to accessing the current directory of the default file system
                {"", "file://" + getCWDAsURIPathString(), true, true},

                //***********************************************************
                // Local file references using a URI with a "file://" scheme.
                //***********************************************************

                {"file:localFile.bam",              "file:localFile.bam",           true, false}, // absolute, opaque (not hierarchical)
                {"file:/localFile.bam",             "file:/localFile.bam",          true, true},  // absolute, hierarchical
                {"file://localFile.bam",            "file://localFile.bam",         true, false}, // file URLs can't have an authority ("localFile.bam")
                {"file:///localFile.bam",           "file:///localFile.bam",        true, true},  // empty authority
                {"file:path/to/localFile.bam",      "file:path/to/localFile.bam",   true, false},
                {"file:/path/to/localFile.bam",     "file:/path/to/localFile.bam",  true, true},
                // "path" looks like an authority, and will be accepted on Windows since it will be interpreted as a UNC authority
                {"file://path/to/localFile.bam",    "file://path/to/localFile.bam", true, SystemUtils.IS_OS_WINDOWS},
                // "localhost" is accepted as a special case authority for "file://" Paths on Windows; but not Linux
                {"file://localhost/to/localFile.bam","file://localhost/to/localFile.bam", true, SystemUtils.IS_OS_WINDOWS},
                {"file:///path/to/localFile.bam",   "file:///path/to/localFile.bam",    true, true},  // empty authority

                //*****************************************************************************
                // Valid URIs which are NOT valid NIO paths (no installed file system provider)
                //*****************************************************************************

                {"gs://file.bam",                   "gs://file.bam",                    true, true},
                {"gs://bucket/file.bam",            "gs://bucket/file.bam",             true, true},
                {"gs:///bucket/file.bam",           "gs:///bucket/file.bam",            true, false},
                {"gs://auth/bucket/file.bam",       "gs://auth/bucket/file.bam",        true, true},
                {"gs://hellbender/test/resources/", "gs://hellbender/test/resources/",  true, true},
                {"gcs://abucket/bucket",            "gcs://abucket/bucket",             false, false},
                {"gendb://somegdb",                 "gendb://somegdb",                  false, false},
                {"chr1:1-100",                      "chr1:1-100",                       false, false},

                //**********************************************************************************************
                // Valid URIs which ARE valid NIO URIs (there *IS* an installed file system provider), but are
                // not actually resolvable as paths because the scheme-specific part is not valid for one reason
                // or another.
                //**********************************************************************************************

                // uri must have a path: jimfs:file.bam
                {"jimfs:file.bam",      "jimfs:file.bam", true, false},
                // java.lang.AssertionError: java.net.URISyntaxException: Expected scheme-specific part at index 6: jimfs:
                {"jimfs:/file.bam",     "jimfs:/file.bam", true, false},
                // java.lang.AssertionError: uri must have a path: jimfs://file.bam
                {"jimfs://file.bam",    "jimfs://file.bam", true, false},
                // java.lang.AssertionError: java.net.URISyntaxException: Expected scheme-specific part at index 6: jimfs:
                {"jimfs:///file.bam",   "jimfs:///file.bam", true, false},
                // java.nio.file.FileSystemNotFoundException: jimfs://root
                {"jimfs://root/file.bam","jimfs://root/file.bam", true, false},

                //*****************************************************************************************
                // Reference that contain characters that require URI-encoding. If the input string is presented
                // with no scheme, it will be be automatically encoded by PathSpecifier, otherwise it
                // must already be URI-encoded.
                //*****************************************************************************************

                // relative (non-URI) reference to a file on the local file system that contains a URI fragment delimiter
                // is automatically URI-encoded
                {joinWithFSSeparator("project", "gvcf-pcr", "23232_1#1", "1.g.vcf.gz"),
                        "file://" + getCWDAsURIPathString() + "project/gvcf-pcr/23232_1%231/1.g.vcf.gz", true, true},
                // URI references with fragment delimiter is not automatically URI-encoded
                {"file:project/gvcf-pcr/23232_1#1/1.g.vcf.gz",  "file:project/gvcf-pcr/23232_1#1/1.g.vcf.gz", true, false},
                {"file:/project/gvcf-pcr/23232_1#1/1.g.vcf.gz", "file:/project/gvcf-pcr/23232_1#1/1.g.vcf.gz", true, false},
                {"file:///project/gvcf-pcr/23232_1%231/1.g.vcf.g", "file:///project/gvcf-pcr/23232_1%231/1.g.vcf.g", true, true},
        };
    }

    @Test(dataProvider = "validPathSpecifiers")
    public void testPathSpecifier(final String referenceString, final String expectedURIString, final boolean hasFileSystemProvider, final boolean isPath) {
        final PathURI pathURI = new GATKPathSpecifier(referenceString);
        Assert.assertNotNull(pathURI);
        Assert.assertEquals(pathURI.getURI().toString(), expectedURIString);
    }

    @Test(dataProvider = "validPathSpecifiers")
    public void testIsNIO(final String referenceString, final String expectedURIString, final boolean hasFileSystemProvider, final boolean isPath) {
        final PathURI pathURI = new GATKPathSpecifier(referenceString);
        Assert.assertEquals(pathURI.hasFileSystemProvider(), hasFileSystemProvider);
    }

    @Test(dataProvider = "validPathSpecifiers")
    public void testIsPath(final String referenceString, final String expectedURIString, final boolean hasFileSystemProvider, final boolean isPath) {
        final PathURI pathURI = new GATKPathSpecifier(referenceString);
        if (isPath) {
            Assert.assertEquals(pathURI.isPath(), isPath, pathURI.getToPathFailureReason());
        } else {
            Assert.assertEquals(pathURI.isPath(), isPath);
        }
    }

    @Test(dataProvider = "validPathSpecifiers")
    public void testToPath(final String referenceString, final String expectedURIString, final boolean hasFileSystemProvider, final boolean isPath) {
        final PathURI pathURI = new GATKPathSpecifier(referenceString);
        if (isPath) {
            final Path path = pathURI.toPath();
            Assert.assertEquals(path != null, isPath, pathURI.getToPathFailureReason());
        } else {
            Assert.assertEquals(pathURI.isPath(), isPath);
        }
    }

    @DataProvider
    public Object[][] invalidPathSpecifiers() {
        return new Object[][] {
                // the nul character is rejected on all of the supported platforms in both local
                // filenames and URIs, so use it to test PathSpecifier constructor failure on all platforms
                {"\0"},
        };
    }

    @Test(dataProvider = "invalidPathSpecifiers", expectedExceptions = {IllegalArgumentException.class})
    public void testPathSpecifierInvalid(final String referenceString) {
        new GATKPathSpecifier(referenceString);
    }

    @DataProvider
    public Object[][] invalidPath() {
        return new Object[][] {
                // valid references that are not valid as a path

                {"file:/project/gvcf-pcr/23232_1#1/1.g.vcf.gz"},    // not encoded
                {"file:project/gvcf-pcr/23232_1#1/1.g.vcf.gz"},     // scheme-specific part is not hierarchical

                // The hadoop file system provider explicitly throws an NPE if no host is specified and HDFS is not
                // the default file system
                //{"hdfs://nonexistent_authority/path/to/file.bam"},  // unknown authority "nonexistent_authority"
                {"hdfs://userinfo@host:80/path/to/file.bam"},           // UnknownHostException "host"

                {"unknownscheme://foobar"},
                {"gendb://adb"},
                {"gcs://abucket/bucket"},

                // URIs with schemes that are backed by an valid NIO provider, but for which the
                // scheme-specific part is not valid.
                {"file://nonexistent_authority/path/to/file.bam"},  // unknown authority "nonexistent_authority"
        };
    }

    @Test(dataProvider = "invalidPath")
    public void testIsPathInvalid(final String invalidPathString) {
        final PathURI htsURI = new GATKPathSpecifier(invalidPathString);
        Assert.assertFalse(htsURI.isPath());
    }

    @Test(dataProvider = "invalidPath", expectedExceptions = {
            IllegalArgumentException.class, FileSystemNotFoundException.class,ProviderNotFoundException.class})
    public void testToPathInvalid(final String invalidPathString) {
        final PathURI htsURI = new GATKPathSpecifier(invalidPathString);
        htsURI.toPath();
    }

    @Test
    public void testInstalledNonDefaultFileSystem() throws IOException {
        // create a jimfs file system and round trip through PathSpecifier/stream
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path outputPath = jimfs.getPath("alternateFileSystemTest.txt");
            doStreamRoundTrip(outputPath.toUri().toString());
        }
    }

    @DataProvider
    public Object[][] inputStreamSpecifiers() {
        return new Object[][]{
                // references that can be resolved to an actual test file that can be read

                //"src/test/resources/org/broadinstitute/hellbender/tools/testTextFile.txt"
                // relative (file) reference to a local file
                {joinWithFSSeparator("src", "test", "resources", "org", "broadinstitute", "hellbender","tools", "testTextFile.txt"), "Test file."},

                // absolute reference to a local file
                {getCWDAsFileReference() + FS_SEPARATOR + joinWithFSSeparator("src", "test", "resources", "org", "broadinstitute", "hellbender", "tools", "testTextFile.txt"), "Test file."},

                // URI reference to a local file, where the path is absolute
                {"file://" + getCWDAsURIPathString() + "src/test/resources/org/broadinstitute/hellbender/tools/testTextFile.txt", "Test file."},

                // reference to a local file with an embedded fragment delimiter ("#") in the name; if the file
                // scheme is included, the rest of the path must already be encoded; if no file scheme is
                // included, the path is encoded by the PathSpecifier class
                {joinWithFSSeparator("src", "test", "resources", "org", "broadinstitute", "hellbender","tools",  "testDirWith#InName", "testTextFile.txt"), "Test file."},
                {"file://" + getCWDAsURIPathString() + "src/test/resources/org/broadinstitute/hellbender/tools/testDirWith%23InName/testTextFile.txt", "Test file."},
        };
    }

    @Test(dataProvider = "inputStreamSpecifiers")
    public void testGetInputStream(final String referenceString, final String expectedFileContents) throws IOException {
        final PathURI htsURI = new GATKPathSpecifier(referenceString);

        try (final InputStream is = htsURI.getInputStream();
             final DataInputStream dis = new DataInputStream(is)) {
            final byte[] actualFileContents = new byte[expectedFileContents.length()];
            dis.readFully(actualFileContents);

            Assert.assertEquals(new String(actualFileContents), expectedFileContents);
        }
    }

    @DataProvider
    public Object[][] outputStreamSpecifiers() {
        return new Object[][]{
                // output URIs that can be resolved to an actual test file
                {IOUtils.createTempPath("testOutputStream", ".txt").toString()},
                {"file://" + getLocalFileAsURIPathString(IOUtils.createTempPath("testOutputStream", ".txt"))},
        };
    }

    @Test(dataProvider = "outputStreamSpecifiers")
    public void testGetOutputStream(final String referenceString) throws IOException {
        doStreamRoundTrip(referenceString);
    }

    @Test
    public void testStdIn() throws IOException {
        final PathURI htsURI = new GATKPathSpecifier(
                SystemUtils.IS_OS_WINDOWS ?
                        "-" :
                        "/dev/stdin");
        try (final InputStream is = htsURI.getInputStream();
             final DataInputStream dis = new DataInputStream(is)) {
            final byte[] actualFileContents = new byte[0];
            dis.readFully(actualFileContents);

            Assert.assertEquals(new String(actualFileContents), "");
        }
    }

    @Test
    public void testStdOut() throws IOException {
        if (SystemUtils.IS_OS_WINDOWS) {
            // stdout is not addressable as a device in the file system namespace on Windows, so skip
            throw new SkipException(("No stdout test on Windows"));
        } else {
            final PathURI pathURI = new GATKPathSpecifier("/dev/stdout");
            try (final OutputStream os = pathURI.getOutputStream();
                 final DataOutputStream dos = new DataOutputStream(os)) {
                dos.write("some stuff".getBytes());
            }
        }
    }

    @DataProvider(name = "getExtensionTestCases")
    public Object[][] getExtensionTestCases() {
        return new Object[][] {
                // input, extension
                {"localFile.bam", ".bam"},
                {"localFile.BAM", ".BAM"},
                {"/localFile.bam", ".bam"},
                {"gs://bucket/aFile.bam", ".bam"},
                {"gs://hellbender/test/resources/aFile.adam", ".adam"},
                {"gs://hellbender/test/resources/aFile.fasta", ".fasta"},
                {"http://bucket/aFile.bam?query=param", ".bam"},

                // getExtension() returns ".gz", but this case also satisfies hasExtension(".fasta.gz")
                {"aFile.fasta.gz", ".gz"},
                // basename is ".fasta"!
                {".fasta.gz", ".gz"},
        };
    }

    @Test(dataProvider = "getExtensionTestCases")
    public void testGetExtension(final String spec, final String expectedExtension) {
        final GATKPathSpecifier pathSpec = new GATKPathSpecifier(spec);
        final String actualExtension = pathSpec.getExtension();

        Assert.assertEquals(actualExtension, expectedExtension);
        // verify that hasExtension(getExtension()) is always true
        Assert.assertTrue(pathSpec.hasExtension(actualExtension));
    }

    @DataProvider(name="negativeGetExtensionTestCases")
    public Object[][] negativeGetExtensionTestCases() {
        return new Object[][]{
                // no extensions
                {""},
                {"/"},
                {"."},
                {"localFile"},
                {"localFile."},
                {"/localFile."},
                {"gs://hellbender/test/resources"},
                {"gs://hellbender/test/resources?query=param"},
                {"gs://hellbender/test/resources/"},
                {"gs://hellbender/test/resources/?query=param"},
        };
    }

    @Test(dataProvider = "negativeGetExtensionTestCases", expectedExceptions={IllegalArgumentException.class})
    public void testNegativeGetExtension(final String spec) {
        new GATKPathSpecifier(spec).getExtension();
    }

    @DataProvider(name = "hasExtensionTestCases")
    public Object[][] hasExtensionTestCases() {
        return new Object[][]{
                // input, extension satisfies "hasExtension"
                {"localFile.bam", ".bam", true },
                {"localFile.BAM", ".BAM", true },
                {"localFile.BAM", ".bam", true },
                {"localFile.bam", ".BAM", true },
                {"/localFile.bam", ".bam", true },
                {"gs://bucket/aFile.bam", ".bam", true },
                {"gs://hellbender/test/resources/aFile.adam", ".adam", true },
                {"gs://hellbender/test/resources/aFile.fasta", ".fasta", true },
                {"http://bucket/aFile.bam?query=param", ".bam", true },

                {"aFile.fasta.gz", ".gz", true },
                {"aFile.fasta.gz", ".fasta.gz", true },
                // basename is ".fasta"!
                {".fasta.gz", ".gz", true },
                {".fasta.gz", ".fasta.gz", true },

                // no extensions
                {"/", ".ext", false },
                {".", ".ext", false },
                {"localFile", ".a", false }, // extension must have length > 1
                {"localFile.", ".a", false },
                {"gs://hellbender/test/resources", ".fasta", false },
                {"gs://hellbender/test/resources?query=param", ".fasta", false },
                {"gs://hellbender/test/resources/", ".fasta", false },
                {"gs://hellbender/test/resources/?query=param", ".fasta", false },
        };
    }

    @Test(dataProvider = "hasExtensionTestCases")
    public void testHasExtension(final String spec, final String extension, final boolean expectedResult) {
        Assert.assertEquals(new GATKPathSpecifier(spec).hasExtension(extension), expectedResult);
    }

    @DataProvider(name = "getBaseNameTestCases")
    public Object[][] getBaseNameTestCases() {
        return new Object[][] {
                // input, baseName
                {"localFile.bam", "localFile"},
                {"localFile.BAM", "localFile"},
                {"/localFile.bam", "localFile"},
                {"gs://bucket/aFile.bam", "aFile"},
                {"gs://hellbender/test/resources/aFile.adam", "aFile"},
                {"gs://hellbender/test/resources/aFile.fasta", "aFile"},
                {"http://bucket/aFile.bam?query=param", "aFile"},

                // This case satisfies hasExtension(".fasta.gz"), but getExtension() returns ".gz".
                {"aFile.fasta.gz", "aFile.fasta"},
                // basename is ".fasta"!
                {".fasta.gz", ".fasta",},
        };
    }

    @Test(dataProvider = "getBaseNameTestCases")
    public void testGetBaseName(final String spec, final String baseName) {
        Assert.assertEquals(new GATKPathSpecifier(spec).getBaseName(), baseName);
    }

    @DataProvider(name="negativeGetBaseNameTestCases")
    public Object[][] negativeGetBaseNameTestCases() {
        return new Object[][]{
                // no extensions
                {"/"},
                {"."},
                {"/."},
                {"/name/.fasta"},
                {"localFile"},
                {"gs://hellbender/test/resources"},
                {"gs://hellbender/test/resources?query=param"},
                {"gs://hellbender/test/resources/"},
                {"gs://hellbender/test/resources/?query=param"},
        };
    }

    @Test(dataProvider = "negativeGetBaseNameTestCases", expectedExceptions = {IllegalArgumentException.class})
    public void testNegativeGetBaseName(final String spec) {
        new GATKPathSpecifier(spec).getBaseName();
    }

    @DataProvider(name="isFastaTestCases")
    public Object[][] isFastaTestCases() {
        final String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";
        return new Object[][] {
                { twoBitRefURL, false },
                { "file://" + twoBitRefURL, false },
                { hg38Reference, true }, // gzipped
                { "file://" + hg38Reference, true }, // gzipped
                { GCS_b37_CHR20_21_REFERENCE_2BIT, false },
                { GCS_b37_CHR20_21_REFERENCE, true },
                // dummy query params at the end to make sure URI.getPath does the right thing
                { GCS_b37_CHR20_21_REFERENCE + "?query=param", true}
        };
    }

    @Test(dataProvider = "isFastaTestCases")
    public void testIsFasta(final String referenceSpec, final boolean expectedIsFasta) {
        Assert.assertEquals(new GATKPathSpecifier(referenceSpec).isFasta(), expectedIsFasta);
    }

    /**
     * Return the string resulting from joining the individual components using the local default
     * file system separator.
     *
     * This is used to create test inputs that are local file references, as would be presented by a
     * user on the platform on which these tests are running.
     */
    private String joinWithFSSeparator(String... parts) {
        return String.join(FileSystems.getDefault().getSeparator(), parts);
    }

    private void doStreamRoundTrip(final String referenceString) throws IOException {
        final String expectedFileContents = "Test contents";

        final PathURI pathURI = new GATKPathSpecifier(referenceString);
        try (final OutputStream os = pathURI.getOutputStream();
             final DataOutputStream dos = new DataOutputStream(os)) {
            dos.write(expectedFileContents.getBytes());
        }

        // read it back in and make sure it matches expected contents
        try (final  InputStream is = pathURI.getInputStream();
             final DataInputStream dis = new DataInputStream(is)) {
            final byte[] actualFileContents = new byte[expectedFileContents.length()];
            dis.readFully(actualFileContents);

            Assert.assertEquals(new String(actualFileContents), expectedFileContents);
        }
    }

    /**
     * Get an absolute reference to the current working directory using local file system syntax and
     * the local file system separator. Used to construct valid local, absolute file references as test inputs.
     *
     * Returns /Users/user/=
     */
    private static String getCWDAsFileReference() {
        return new File(".").getAbsolutePath();
    }

    /**
     * Get the current working directory as a locally valid, hierarchical URI string. Used to
     * construct expected URI string values for test inputs that are local file references.
     */
    private String getCWDAsURIPathString() {
        return getLocalFileAsURIPathString(Paths.get("."));
    }

    /**
     * Get just the path part of the URI representing the current working directory. Used
     * to construct expected URI string values for test inputs that specify a file in the
     * root of the local file system.
     *
     * This will return a string of the form "/" on *nix and "d:/" on Windows (its a URI string).
     */
    private String getRootDirectoryAsURIPathString() {
        return getLocalFileAsURIPathString(Paths.get(FS_SEPARATOR));
    }

    /**
     * Get just the path part of the URI representing a file on the local file system. This
     * uses java.io.File to get a locally valid file reference, which is then converted to
     * a URI.
     */
    private String getLocalFileAsURIPathString(final Path localPath) {
        return localPath.toUri().normalize().getPath();
    }

}
