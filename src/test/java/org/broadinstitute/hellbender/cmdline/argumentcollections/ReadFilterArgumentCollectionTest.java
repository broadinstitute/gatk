package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

// - walker and spark tool integration tests

public class ReadFilterArgumentCollectionTest {

    class TestArgCollection {
        @ArgumentCollection
        public ReadFilterArgumentCollection rfc = new ReadFilterArgumentCollection();
        public TestArgCollection(){};
        public TestArgCollection (Consumer<ReadFilterArgumentCollection> t) {
            t.accept(this.rfc);
        }
    };

    @DataProvider(name="filtersWithArguments")
    public Object[][] badArguments(){
        return new Object[][]{
                { ReadFilterArgumentCollection.CommandLineReadFilter.LIBRARY, "--library", "fakeLibrary" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.MAPPING_QUALITY, "--mappingQuality", "25" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.PLATFORM, "--PLFilterName", "fakePlatform" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.PLATFORM_UNIT, "--blackListedLanes", "fakeUnit" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_GROUP, "--readGroup", "fakeGroup" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_GROUP_BLACK_LISTED, "--blackList", "tg:sub"},
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_NAME, "--readName", "fakeRead" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_LENGTH, "--maxReadLength", "10" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.SAMPLE, "--sample", "fakeSample" }
        };
    }

    // verify that all filters that require arguments and have no defaults complain if the
    // arguments are not set on the command line
    @Test(dataProvider = "filtersWithArguments", expectedExceptions = UserException.CommandLineException.class)
    public void testRequiresArguments(
            ReadFilterArgumentCollection.CommandLineReadFilter filter,
            String argName,   //unused
            String argValue)  //unused
    {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        String[] args = { "--readFilter", filter.name() }; // no args, just enable filters

        clp.parseArguments(System.out, args);

        // we need to render the filter to trigger execution of the validation code
        testArgs.rfc.getMergedCommandLineFilters(Arrays.asList())
                .stream()
                .map(f-> f.getFilterInstance(new SAMFileHeader(), testArgs.rfc))
                .collect(Collectors.toList());
    }

    // verify that all filters that require arguments throw if the arguments are passed but the
    // corresponding filter is not enabled
    @Test(dataProvider = "filtersWithArguments", expectedExceptions = UserException.CommandLineException.class)
    public void testArgumentsRequireFilter(
            ReadFilterArgumentCollection.CommandLineReadFilter filter,
            String argName,
            String argValue)
    {
        TestArgCollection tArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(tArgs);
        String[] args = { argName, argValue }; // no read filter set

        // dependsOn errors are caught by the command line parser
        clp.parseArguments(System.out, args);
    }

    @DataProvider(name="validationFailures")
    public Object[][] validationFailures(){
        return new Object[][]{
                { ReadFilterArgumentCollection.CommandLineReadFilter.LIBRARY, "--library", "" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.MAPPING_QUALITY, "--mappingQuality", "1000" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_GROUP, "--readGroup", "" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_GROUP_BLACK_LISTED, "--blackList", "zzz:sub"},
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_NAME, "--readName", "" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_LENGTH, "--maxReadLength", "-1" }
        };
    }

    // test that all filters that require validation fail given bad arguments
    @Test(dataProvider = "validationFailures", expectedExceptions =
            {UserException.CommandLineException.class, UserException.class})
    public void testArgumentFailsValidations(
            ReadFilterArgumentCollection.CommandLineReadFilter filter,
            String argName,
            String argValue)
    {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        String[] args = {
                "--readFilter", filter.name(),
                argName, argValue
        };

        clp.parseArguments(System.out, args);
        testArgs.rfc.getMergedCommandLineFilters(Arrays.asList())
                .stream()
                .map(f-> f.getFilterInstance(new SAMFileHeader(), testArgs.rfc))
                .collect(Collectors.toList());
    }

    @DataProvider(name="filtersWithGoodArguments")
    public Object[][] filtersWithGoodArguments(){
        return new Object[][]{
                { ReadFilterArgumentCollection.CommandLineReadFilter.LIBRARY,
                        (Consumer<SetupTest>) this::setupLibraryTest, "--library", "Foo" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.MAPPING_QUALITY,
                        (Consumer<SetupTest>) this::setupMappingQualityTest, "--mappingQuality", "255" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_GROUP,
                        (Consumer<SetupTest>) this::setupReadGroupTest, "--readGroup", "fred" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.READ_NAME,
                        (Consumer<SetupTest>) this::setupReadNameTest, "--readName", "fred" },
                { ReadFilterArgumentCollection.CommandLineReadFilter.SAMPLE,
                        (Consumer<SetupTest>) this::setupSampleTest, "--sample", "fred" }
        };
    }

    @Test(dataProvider = "filtersWithGoodArguments")
    public void testFiltersWithArguments(
            ReadFilterArgumentCollection.CommandLineReadFilter filter,
            Consumer<SetupTest> setup,
            String argName,
            String argValue)
    {
        TestArgCollection tArgs = new TestArgCollection();

        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineParser(tArgs);
        String[] args = {
                "--readFilter", filter.name(),
                argName, argValue
        };
        clp.parseArguments(System.out, args);
        ReadFilter rf = instantiateFilter(tArgs.rfc, header);

        // to ensure that the filter is actually working, verify that the test record
        // we're using fails the filter test *before* we set it up to pass the filter
        Assert.assertFalse(rf.test(read));

        // setup the header and read for this test
        setup.accept(new SetupTest(header, read, argValue));

        Assert.assertTrue(rf.test(read));
    }

    @Test
    public void testNoReadFiltersSpecified() {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        clp.parseArguments(System.out, new String[]{});
        List<ReadFilterArgumentCollection.CommandLineReadFilter> filters =
                testArgs.rfc.getMergedCommandLineFilters(Arrays.asList());
        Assert.assertEquals(filters.size(), 0);
    }

    @Test
    public void testDisableAllReadFilters() {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED.name(), "--disableAllReadFilters"});
        List<ReadFilterArgumentCollection.CommandLineReadFilter> filters =
                testArgs.rfc.getMergedCommandLineFilters(Arrays.asList(
                        ReadFilterArgumentCollection.CommandLineReadFilter.CIGAR_IS_SUPPORTED,
                        ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR
                ));
        Assert.assertEquals(filters.size(), 0);
    }

    @Test
    public void testDisableOneReadFilter() {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED.name(),
                "--RF", ReadFilterArgumentCollection.CommandLineReadFilter.HAS_MATCHING_BASES_AND_QUALS.name(),
                "-disableReadFilter", ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR.name()});

        List<ReadFilterArgumentCollection.CommandLineReadFilter> filters =
                testArgs.rfc.getMergedCommandLineFilters(Arrays.asList(
                        ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR
                ));

        Assert.assertEquals(filters.size(), 2);
        Assert.assertTrue(filters.contains(ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED));
        Assert.assertTrue(filters.contains(ReadFilterArgumentCollection.CommandLineReadFilter.HAS_MATCHING_BASES_AND_QUALS));
    }

    @Test
    public void testDisableMultipleReadFilters() {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED.name(),
                "-disableReadFilter", ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR.name()});

        List<ReadFilterArgumentCollection.CommandLineReadFilter> filters =
                testArgs.rfc.getMergedCommandLineFilters(Arrays.asList(
                        ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR,
                        ReadFilterArgumentCollection.CommandLineReadFilter.HAS_MATCHING_BASES_AND_QUALS
                ));

        Assert.assertEquals(filters.size(), 2);
        Assert.assertTrue(filters.contains(ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED));
        Assert.assertTrue(filters.contains(ReadFilterArgumentCollection.CommandLineReadFilter.HAS_MATCHING_BASES_AND_QUALS));
    }

    @Test
    public void testPreserveSpecifiedOrder() {
        TestArgCollection testArgs = new TestArgCollection();
        CommandLineParser clp = new CommandLineParser(testArgs);
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED.name(),
                "--RF", ReadFilterArgumentCollection.CommandLineReadFilter.HAS_MATCHING_BASES_AND_QUALS.name()});

        List<ReadFilterArgumentCollection.CommandLineReadFilter> filters =
                testArgs.rfc.getMergedCommandLineFilters(Arrays.asList(
                        ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR
                ));

        // defaults first, then command line, in order
        Assert.assertEquals(filters.size(), 3);
        Assert.assertEquals(filters.get(0), ReadFilterArgumentCollection.CommandLineReadFilter.GOOD_CIGAR);
        Assert.assertEquals(filters.get(1), ReadFilterArgumentCollection.CommandLineReadFilter.MAPPED);
        Assert.assertEquals(filters.get(2), ReadFilterArgumentCollection.CommandLineReadFilter.HAS_MATCHING_BASES_AND_QUALS);
    }

    @Test
    public void testMultipleFilters() {
        TestArgCollection tArgs = new TestArgCollection();

        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineParser(tArgs);
        String[] args = {
                "--readFilter", ReadFilterArgumentCollection.CommandLineReadFilter.READ_LENGTH.name(),
                "--minReadLength", "10",
                "--maxReadLength", "20",
                "--readFilter", ReadFilterArgumentCollection.CommandLineReadFilter.READ_NAME.name(),
                "--readName", "fred"
        };
        clp.parseArguments(System.out, args);
        ReadFilter rf = instantiateFilter(tArgs.rfc, header);

        Assert.assertFalse(rf.test(read));
        read.setName("fred");
        read.setBases(new byte[15]);
        Assert.assertTrue(rf.test(read));
    }

    @Test
    public void testReadLengthFilter() {
        TestArgCollection tArgs = new TestArgCollection();

        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineParser(tArgs);
        String[] args = {
                "--readFilter", ReadFilterArgumentCollection.CommandLineReadFilter.READ_LENGTH.name(),
                "--minReadLength", "10",
                "--maxReadLength", "20"
        };
        clp.parseArguments(System.out, args);
        ReadFilter rf = instantiateFilter(tArgs.rfc, header);

        read.setBases(new byte[5]);
        Assert.assertFalse(rf.test(read));
        read.setBases(new byte[25]);
        Assert.assertFalse(rf.test(read));
        read.setBases(new byte[15]);
        Assert.assertTrue(rf.test(read));
    }

    private ReadFilter instantiateFilter(final ReadFilterArgumentCollection rfc, final SAMFileHeader header) {
        final List<ReadFilterArgumentCollection.CommandLineReadFilter> filters =
                rfc.getMergedCommandLineFilters(Arrays.asList());
        return filters.stream()
                .map(f-> f.getFilterInstance(header, rfc))
                .reduce(
                        ReadFilterLibrary.ALLOW_ALL_READS,
                        (f1, f2) -> f1.and(f2)
                );
    }

    private static final int CHR_COUNT = 2;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private SAMFileHeader createHeaderWithReadGroups() {
        return ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);
    }

    /**
     * Creates a read record.
     *
     * @param header header for the new record
     * @param cigar the new record CIGAR.
     * @param group the new record group index that must be in the range \
     *              [0,{@link #GROUP_COUNT})
     * @param reference the reference sequence index (0-based)
     * @param start the start position of the read alignment in the reference
     *              (1-based)
     * @return never <code>null</code>
     */
    private GATKRead createRead(final SAMFileHeader header, final Cigar cigar, final int group,
                                final int reference, final int start ) {
        final GATKRead record = ArtificialReadUtils.createArtificialRead(header, cigar);
        record.setPosition(header.getSequence(reference).getSequenceName(), start);
        record.setReadGroup(header.getReadGroups().get(group).getReadGroupId());
        return record;
    }

    private GATKRead simpleGoodRead( final SAMFileHeader header ) {
        final String cigarString = "101M";
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        GATKRead read = createRead(header, cigar, 1, 0, 10);
        read.setMappingQuality(50);
        return read;
    }

    // object for holding test setup params to use as a type for test lambdas
    public class SetupTest {
        public SAMFileHeader hdr;
        public GATKRead read;
        public String argValue;

        public SetupTest(SAMFileHeader hdr, GATKRead read, String argValue) {
            this.hdr = hdr;
            this.read = read;
            this.argValue = argValue;
        }
    }

    // setup methods to initialize reads for individual filter tests
    private void setupLibraryTest(SetupTest setup) {
        setup.hdr.getReadGroup(setup.read.getReadGroup()).setLibrary(setup.argValue);
    }

    private void setupMappingQualityTest(SetupTest setup) {
        setup.read.setMappingQuality(Integer.parseInt(setup.argValue));
    }

    private void setupReadGroupTest(SetupTest setup) {
        setup.read.setReadGroup(setup.argValue);
    }

    private void setupReadNameTest(SetupTest setup) {
        setup.read.setName(setup.argValue);
    }

    private void setupSampleTest(SetupTest setup) {
        setup.hdr.getReadGroup(setup.read.getReadGroup()).setSample(setup.argValue);
    }

}
