package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeSuite;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;


/**
 * This is the base test class for all of our test cases.  All test cases should extend from this
 * class; it sets up the logger, and resolves the location of directories that we rely on.
 */
public abstract class BaseTest {

    static {
        // set properties for the local Spark runner
        System.setProperty("dataflow.spark.test.reuseSparkContext", "true");
        SparkContextFactory.enableTestSparkContext();
    }

    @BeforeSuite
    public void setTestVerbosity(){
        LoggingUtils.setLoggingLevel(Log.LogLevel.WARNING);
    }

    public static final Logger logger = LogManager.getLogger("org.broadinstitute.gatk");

    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    public static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";

    private static final String publicTestDirRelative = "src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";
    public static final String publicTestDirRoot = publicTestDir.replace(publicTestDirRelative, "");

    public static final String moduleTestDataDir = publicTestDir + "org/broadinstitute/hellbender/";
    public static final String toolTestDataDir = moduleTestDataDir + "tools/";
    public static final String bqsrTestDataDir = toolTestDataDir + "BQSR/";

    public static final String GCS_GATK_TEST_RESOURCES = "gs://hellbender/test/resources/";

    public static final String GCS_b37_REFERENCE_2BIT = GCS_GATK_TEST_RESOURCES + "benchmark/human_g1k_v37.2bit";
    public static final String GCS_b37_CHR20_21_REFERENCE_2BIT = GCS_GATK_TEST_RESOURCES + "human_g1k_v37.20.21.2bit";

    /**
     * LARGE FILES FOR TESTING (MANAGED BY GIT LFS)
     */
    public static final String largeFileTestDir = new File(publicTestDir, "large").getAbsolutePath() + "/";

    // All of chromosomes 20 and 21 from the b37 reference
    public static final String b37_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.fasta";

    public static final String b37_2bit_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.2bit";

    // All of chromosomes 20 and 21 from the b38 reference
    public static final String b38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";

    // ~600,000 reads from chromosomes 20 and 21 of an NA12878 WGS bam aligned to b37, plus ~50,000 unmapped reads
    public static final String NA12878_20_21_WGS_bam = largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";

    // Variants from a DBSNP 138 VCF overlapping the reads in NA12878_20_21_WGS_bam
    public static final String dbsnp_138_b37_20_21_vcf = largeFileTestDir + "dbsnp_138.b37.20.21.vcf";

    // Variants from a DBSNP 138 VCF form the first 65Mb of chr1
    public static final String dbsnp_138_b37_1_65M_vcf = largeFileTestDir + "dbsnp_138.b37.1.1-65M.vcf";

    public static final String WGS_B37_CH20_1M_1M1K_BAM = "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
    public static final String DBSNP_138_B37_CH20_1M_1M1K_VCF = "dbsnp_138.b37.excluding_sites_after_129.ch20.1m-1m1k.vcf";

    /**
     * END OF LARGE FILES FOR TESTING
     */

    public static final String NA12878_chr17_1k_BAM = publicTestDir + "NA12878.chr17_69k_70k.dictFix.bam";
    public static final String NA12878_chr17_1k_CRAM = publicTestDir + "NA12878.chr17_69k_70k.dictFix.cram";
    public static final String v37_chr17_1Mb_Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";

    public static final String hg19_chr1_1M_Reference = publicTestDir + "Homo_sapiens_assembly19_chr1_1M.fasta";
    public static final String hg19_chr1_1M_dbSNP = publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf";

    // the following file has been modified such that the first chromosome length is 1M; this is sometimes
    // required due to sequence dictionary validation, since a reference FASTA with only 1M bases is used
    public static final String hg19_chr1_1M_dbSNP_modified = publicTestDir + "HSA19.dbsnp135.chr1_1M.exome_intervals.modified.vcf";

    public static final String hg19_chr1_1M_exampleVCF = publicTestDir + "joint_calling.chr1_1M.1kg_samples.10samples.noINFO.vcf";
    public static final String hg19MiniReference = publicTestDir + "hg19mini.fasta";

    public static final String exampleFASTA = publicTestDir + "exampleFASTA.fasta";
    public static final String exampleReference = hg19MiniReference;
    public static final String hg19MiniIntervalFile = publicTestDir + "hg19mini.interval_list";

    public static final String hg19MicroReference = publicTestDir + "hg19micro.fasta";

    /**
     * BQSR Test Files
     */

    public static final String BQSR_WGS_B37_CH20_21_10M_100_CRAM = bqsrTestDataDir +
            "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.cram";

    public CachingIndexedFastaSequenceFile hg19ReferenceReader;
    public GenomeLocParser hg19GenomeLocParser;

    // used to seed the genome loc parser with a sequence dictionary
    protected SAMFileHeader hg19Header;

    /**
     * name of the google cloud project that stores the data and will run the code
     * @return HELLBENDER_TEST_PROJECT env. var if defined, throws otherwise.
     */
    public static String getGCPTestProject() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_PROJECT");
    }

    /**
     * API key for HELLBENDER_TEST_PROJECT
     * @return HELLBENDER_TEST_APIKEY env. var if defined, throws otherwise.
     */
    public static String getGCPTestApiKey() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_APIKEY");
    }

    /**
     * A writable GCS path where java files can be cached and temporary test files can be written
     * @return HELLBENDER_TEST_STAGING env. var if defined, throws otherwise.
     */
    public static String getGCPTestStaging() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_STAGING");
    }

    /**
     *  A GCS path where the test inputs are stored.
     *
     *  The value of HELLBENDER_TEST_INPUTS should end in a "/" (for example, "gs://hellbender/test/resources/")
     *  
     *  @return HELLBENDER_TEST_INPUTS env. var if defined, throws otherwise.
     */
    public static String getGCPTestInputPath() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_INPUTS");
    }

    /**
     * A local path where the service account credentials are stored
     * @return GOOGLE_APPLICATION_CREDENTIALS env. var if defined, throws otherwise.
     */
    public static String getGoogleServiceAccountKeyPath() {
      return getNonNullEnvironmentVariable("GOOGLE_APPLICATION_CREDENTIALS");
    }

    protected static String getNonNullEnvironmentVariable(String envVarName) {
        String value = System.getenv(envVarName);
        if (null == value) {
            throw new UserException("For this test, please define environment variable \""+envVarName+"\"");
        }
        return value;
    }

    @BeforeClass
    public void initGenomeLocParser() throws FileNotFoundException {
        hg19ReferenceReader = new CachingIndexedFastaSequenceFile(new File(hg19MiniReference));
        hg19Header = new SAMFileHeader();
        hg19Header.setSequenceDictionary(hg19ReferenceReader.getSequenceDictionary());
        hg19GenomeLocParser = new GenomeLocParser(hg19ReferenceReader);
    }

    protected List<GenomeLoc> intervalStringsToGenomeLocs( String... intervals) {
        return intervalStringsToGenomeLocs(Arrays.asList(intervals));
    }

    protected List<GenomeLoc> intervalStringsToGenomeLocs( List<String> intervals ) {
        List<GenomeLoc> locs = new ArrayList<>();
        for (String interval: intervals)
            locs.add(hg19GenomeLocParser.parseGenomeLoc(interval));
        return Collections.unmodifiableList(locs);
    }

    /**
     * Returns the location of the resource directory for the command line program.
     */
    public String getToolTestDataDir(){
        return "src/test/resources/" +this.getClass().getPackage().getName().replace(".","/") +"/" + getTestedClassName() + "/";
    }

    /**
     * Returns the name of the class being tested.
     * The default implementation takes the simple name of the test class and removes the trailing "Test".
     * Override if needed.
     */
    public String getTestedClassName(){
        if (getClass().getSimpleName().contains("IntegrationTest"))
            return getClass().getSimpleName().replaceAll("IntegrationTest$", "");
        else if (getClass().getSimpleName().contains("UnitTest"))
            return getClass().getSimpleName().replaceAll("UnitTest$", "");
        else
            return getClass().getSimpleName().replaceAll("Test$", "");
    }

    /**
     * @param fileName the name of a file
     * @return a File resolved using getToolTestDataDir as the parent and fileName
     */
    public File getTestFile(String fileName) {
        return new File(getToolTestDataDir(), fileName);
    }

    /**
     * Simple generic utility class to creating TestNG data providers:
     *
     * 1: inherit this class, as in
     *
     *      private class SummarizeDifferenceTest extends TestDataProvider {
     *         public SummarizeDifferenceTest() {
     *           super(SummarizeDifferenceTest.class);
     *         }
     *         ...
     *      }
     *
     * Provide a reference to your class to the TestDataProvider constructor.
     *
     * 2: Create instances of your subclass.  Return from it the call to getTests, providing
     * the class type of your test
     *
     * <code>
     * {@literal @}DataProvider(name = "summaries")
     * public Object[][] createSummaries() {
     *   new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
     *   new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
     *   return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
     * }
     * </code>
     *
     * This class magically tracks created objects of this
     */
    public static class TestDataProvider {
        private static final Map<Class<?>, List<Object>> tests = new LinkedHashMap<>();
        protected String name;

        /**
         * Create a new TestDataProvider instance bound to the class variable C
         */
        public TestDataProvider(Class<?> c, String name) {
            if ( ! tests.containsKey(c) )
                tests.put(c, new ArrayList<>());
            tests.get(c).add(this);
            this.name = name;
        }

        public TestDataProvider(Class<?> c) {
            this(c, "");
        }

        public void setName(final String name) {
            this.name = name;
        }

        /**
         * Return all of the data providers in the form expected by TestNG of type class C
         */
        public static Object[][] getTests(Class<?> c) {
            List<Object[]> params2 = new ArrayList<>();
            for ( Object x : tests.get(c) ) params2.add(new Object[]{x});
            return params2.toArray(new Object[][]{});
        }

        @Override
        public String toString() {
            return "TestDataProvider("+name+")";
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(final String name, final String extension) {
        return IOUtils.createTempFile(name, extension);
    }

    /**
     * Return a File object representing a file with the given name and extension that is guaranteed not to exist.
     * @param fileNameWithExtension
     * @return File object representing a file that is guaranteed not to exist
     */
    public static File getSafeNonExistentFile(final String fileNameWithExtension) {
        final File tempDir = createTempDir("nonExistentFileHolder");
        final File nonExistingFile = new File(tempDir, fileNameWithExtension);
        return nonExistingFile;
    }

    /**
     * Return a Path object representing a file with the given name and extension that is guaranteed not to exist.
     * @param fileNameWithExtension
     * @return Path object representing a file that is guaranteed not to exist
     */
    public static Path getSafeNonExistentPath(final String fileNameWithExtension) {
        return getSafeNonExistentFile(fileNameWithExtension).toPath();
    }

    /**
     * Creates an empty temp directory which will be deleted on exit after tests are complete
     *
     * @param prefix prefix for the directory name
     * @return an empty directory starting with prefix which will be deleted after the program exits
     */
    public static File createTempDir(final String prefix){
        final File dir = IOUtils.tempDir(prefix, "");
        IOUtils.deleteRecursivelyOnExit(dir);
        return dir;
    }

    /**
     * Log this message so that it shows up inline during output as well as in html reports
     */
    public static void log(final String message) {
        Reporter.log(message, true);
    }

    private static final double DEFAULT_FLOAT_TOLERANCE = 1e-1;


    /**
     * Checks whether two double array contain the same values or not.
     * @param actual actual produced array.
     * @param expected expected array.
     * @param tolerance maximum difference between double value to be consider equivalent.
     */
    protected static void assertEqualsDoubleArray(final double[] actual, final double[] expected, final double tolerance) {
        if (expected == null)
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual.length,expected.length,"array length");
        }
        for (int i = 0; i < actual.length; i++)
            Assert.assertEquals(actual[i],expected[i],tolerance,"array position " + i);
    }

    /**
     * Checks whether two long arrays contain the same values or not.
     * @param actual actual produced array.
     * @param expected expected array.
     */
    protected static void assertEqualsLongArray(final long[] actual, final long[] expected) {
        if (expected == null)
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual.length, expected.length,"array length ");
        }
        for (int i = 0; i < actual.length; i++)
            Assert.assertEquals(actual[i],expected[i],"array position " + i);
    }

    public static void assertEqualsDoubleSmart(final Object actual, final Double expected, final double tolerance) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double) (Double) actual, (double) expected, tolerance);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected) {
        assertEqualsDoubleSmart(actual, expected, DEFAULT_FLOAT_TOLERANCE);
    }

    public static <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new LinkedHashSet<>(actual);
        final Set<T> expectedSet = new LinkedHashSet<>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance) {
        assertEqualsDoubleSmart(actual, expected, tolerance, null);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance, final String message) {
        if ( Double.isNaN(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isNaN(actual), "expected is nan, actual is not");
        else if ( Double.isInfinite(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isInfinite(actual), "expected is infinite, actual is not");
        else {
            final double delta = Math.abs(actual - expected);
            final double ratio = Math.abs(actual / expected - 1.0);
            Assert.assertTrue(delta < tolerance || ratio < tolerance, "expected = " + expected + " actual = " + actual
                    + " not within tolerance " + tolerance
                    + (message == null ? "" : "message: " + message));
        }
    }


    /**
     * captures {@link System#out} while runnable is executing
     * @param runnable a code block to execute
     * @return everything written to {@link System#out} by runnable
     */
    public static String captureStdout(Runnable runnable){
        return captureSystemStream(runnable, System.out, System::setOut);
    }


    /**
     * captures {@link System#err} while runnable is executing
     * @param runnable a code block to execute
     * @return everything written to {@link System#err} by runnable
     */
    public static String captureStderr(Runnable runnable){
        return captureSystemStream(runnable, System.err, System::setErr);
    }

    private static String captureSystemStream(Runnable runnable,  PrintStream stream, Consumer<? super PrintStream> setter){
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        setter.accept(new PrintStream(out));
        try {
            runnable.run();
        } finally{
            setter.accept(stream);
        }
        return out.toString();
    }

    public static void assertContains(String actual, String expectedSubstring){
        Assert.assertTrue(actual.contains(expectedSubstring),  expectedSubstring +" was not found in " + actual+ ".");
    }

    public static <T> void assertCondition(Iterable<T> actual, Iterable<T> expected, BiConsumer<T,T> assertion){
        final Iterator<T> iterActual = actual.iterator();
        final Iterator<T> iterExpected = expected.iterator();
        while(iterActual.hasNext() && iterExpected.hasNext()){
            assertion.accept(iterActual.next(), iterExpected.next());
        }
        if (iterActual.hasNext()){
            Assert.fail("actual is longer than expected with at least one additional element: " + iterActual.next());
        }
        if (iterExpected.hasNext()){
            Assert.fail("actual is shorter than expected, missing at least one element: " + iterExpected.next());
        }
    }

}

