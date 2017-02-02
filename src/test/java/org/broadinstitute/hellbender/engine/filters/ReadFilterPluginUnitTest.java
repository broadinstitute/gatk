package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;


public class ReadFilterPluginUnitTest {

    //Filters with arguments to verify filter test method actually filters
    @DataProvider(name="filtersWithFilteringArguments")
    public Object[][] filtersWithGoodArguments(){
        return new Object[][]{
                { LibraryReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupLibraryTest, "--library", "Foo" },
                { MappingQualityReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupMappingQualityTest, "--minimumMappingQuality", "255" },
                { MappingQualityReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupMappingQualityTest, "--minimumMappingQuality", "60" },
                { ReadGroupReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupReadGroupTest, "--keepReadGroup", "fred" },
                { ReadNameReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupReadNameTest, "--readName", "fred" },
                { SampleReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupSampleTest, "--sample", "fred" }
        };
    }

    @Test(dataProvider = "filtersWithFilteringArguments")
    public void testFilterArguments(
            final String filter,
            final Consumer<SetupTest> setup,
            final String argName,
            final String argValue)
    {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));

        String[] args = {
                "--readFilter", filter,
                argName, argValue
        };
        clp.parseArguments(System.out, args);
        ReadFilter rf = instantiateFilter(clp, header);

        // to ensure that the filter is actually working, verify that the test record
        // we're using fails the filter test *before* we set it up to pass the filter
        Assert.assertFalse(rf.test(read));

        // setup the header and read for the test
        setup.accept(new SetupTest(header, read, argValue));

        Assert.assertTrue(rf.test(read));
    }

    @DataProvider(name="filtersWithRequiredArguments")
    public Object[][] filtersWithRequiredArguments(){
        return new Object[][]{
                { LibraryReadFilter.class.getSimpleName(), "--library", "fakeLibrary" },
                { PlatformReadFilter.class.getSimpleName(), "--PLFilterName", "fakePlatform" },
                { PlatformUnitReadFilter.class.getSimpleName(), "--blackListedLanes", "fakeUnit" },
                { ReadGroupReadFilter.class.getSimpleName(), "--keepReadGroup", "fakeGroup" },
                { ReadGroupBlackListReadFilter.class.getSimpleName(), "--blackList", "tg:sub"},
                { ReadNameReadFilter.class.getSimpleName(), "--readName", "fakeRead" },
                { ReadLengthReadFilter.class.getSimpleName(), "--maxReadLength", "10" },
                { SampleReadFilter.class.getSimpleName(), "--sample", "fakeSample" }
        };
    }

    // fail if a filter with required arguments is specified without corresponding arguments
    @Test(dataProvider = "filtersWithRequiredArguments", expectedExceptions = CommandLineException.MissingArgument.class)
    public void testDependentFilterArguments(
            final String filter,
            final String argName,   //unused
            final String argValue)  //unused
    {
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        String[] args = {
                "--readFilter", filter  // no args, just enable filters
        };

        clp.parseArguments(System.out, args);
    }

    // fail if a filter's arguments are passed but the filter itself is not enabled
    @Test(dataProvider = "filtersWithRequiredArguments", expectedExceptions = CommandLineException.class)
    public void testDanglingFilterArguments(
            final String filter,
            final String argName,
            final String argValue)
    {
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));

        String[] args = { argName, argValue }; // no read filter set

        // no need to instantiate the filters - dependency errors are caught by the command line parser
        clp.parseArguments(System.out, args);
    }

    @Test
    public void testNoFiltersSpecified() {
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        clp.parseArguments(System.out, new String[]{});

        // get the command line read filters
        final GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        final List<ReadFilter> readFilters = readFilterPlugin.getAllInstances();
        Assert.assertEquals(readFilters.size(), 0);
    }

    @Test
    public void testEnableMultipleFilters() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        String[] args = {
                "--readFilter", ReadLengthReadFilter.class.getSimpleName(),
                "--minReadLength", "10",
                "--maxReadLength", "102",
                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                "--readName", "fred"
        };
        clp.parseArguments(System.out, args);
        ReadFilter rf = instantiateFilter(clp, header);

        // default read name is rejected by the readName filter
        Assert.assertFalse(rf.test(read));

        String readName = read.getName();
        read.setName("fred");
        Assert.assertTrue(rf.test(read)); // accepted

        // trigger the read length filter to reject
        read.setBases(new byte[150]);
        Assert.assertFalse(rf.test(read));
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testEnableNonExistentFilter() {
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        clp.parseArguments(System.out, new String[] {"--RF", "fakeFilter"});
    }

    @DataProvider(name = "duplicateFilters")
    public Object[][] duplicateFilters() {
        return new Object[][] {
                {new String[] {
                        "--RF", ReadLengthReadFilter.class.getSimpleName(),
                        "--RF", ReadLengthReadFilter.class.getSimpleName()}},
                {new String[] {
                        "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                        "--RF", ReadLengthReadFilter.class.getSimpleName(),
                        "--RF", ReadLengthReadFilter.class.getSimpleName()}
                }
        };
    }

    @Test(dataProvider = "duplicateFilters", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testEnableDuplicateFilter(final String[] arguments) {
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        clp.parseArguments(System.out, arguments);
    }

    @Test
    public void testDisableOneFilter() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(
                        Collections.singletonList(ReadFilterLibrary.GOOD_CIGAR)
                )));
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "--RF", ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                "-disableReadFilter", ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()});

        // get the command line read filters
        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        List<ReadFilter> readFilters = readFilterPlugin.getAllInstances();

        Assert.assertEquals(readFilters.size(), 2);
        Assert.assertEquals(readFilters.get(0).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());
        Assert.assertEquals(readFilters.get(1).getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName());

        Assert.assertEquals(readFilterPlugin.disableFilters.size(), 1);
        Assert.assertTrue(readFilterPlugin.disableFilters.contains(
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()));
        Assert.assertTrue(readFilterPlugin.isDisabledFilter(
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()));
    }

    @Test
    public void testDisableMultipleFilters() {
        List<ReadFilter> defaultFilters = new ArrayList<>();
        defaultFilters.add(ReadFilterLibrary.GOOD_CIGAR);
        defaultFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaultFilters)));
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "-disableReadFilter", ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName(),
                "-disableReadFilter", ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName()});

        // get the command line read filters
        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        List<ReadFilter> readFilters = readFilterPlugin.getAllInstances();

        Assert.assertEquals(readFilters.size(), 1);
        Assert.assertEquals(readFilters.get(0).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());

        Assert.assertEquals(readFilterPlugin.disableFilters.size(), 2);
        Assert.assertTrue(readFilterPlugin.disableFilters.contains(
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()));
        Assert.assertTrue(readFilterPlugin.disableFilters.contains(
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName()));

        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        Assert.assertEquals(
                rf.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());
    }

    @Test
    public void testDisableAllFilters() {
        List<ReadFilter> defaultFilters = new ArrayList<>();
        defaultFilters.add(ReadFilterLibrary.GOOD_CIGAR);
        defaultFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaultFilters)));
        clp.parseArguments(System.out, new String[] {
                "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "--disableAllReadFilters"});

        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        Assert.assertTrue(readFilterPlugin.disableAllReadFilters);

        List<ReadFilter> readFilters = readFilterPlugin.getAllInstances();
        Assert.assertEquals(readFilters.size(), 1); // allow all

        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        Assert.assertEquals(
                rf.getClass().getSimpleName(),
                ReadFilterLibrary.ALLOW_ALL_READS.getClass().getSimpleName());
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testDisabledDefaultWithArgsProvided() {
        // test for arguments provided for a default filter that is also disabled
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(
                        Collections.singletonList(new SampleReadFilter())
                )));
        clp.parseArguments(System.out, new String[] {
                "--sample", "fred",
                "-disableReadFilter", SampleReadFilter.class.getSimpleName()});

        // get the command line read filters
        clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testEnableDisableConflict() {
        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        clp.parseArguments(System.out, new String[] {
                "--RF", "GoodCigarReadFilter",
                "--disableReadFilter", "GoodCigarReadFilter"});
    }

    @Test
    public void testPreserveCommandLineOrder() {
        List<ReadFilter> orderedDefaults = new ArrayList<>();
        orderedDefaults.add(new WellformedReadFilter());
        orderedDefaults.add(ReadFilterLibrary.HAS_READ_GROUP);
        orderedDefaults.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);

        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
              Collections.singletonList(new GATKReadFilterPluginDescriptor(orderedDefaults)));
        clp.parseArguments(System.out, new String[] {
                "-readFilter", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "-readFilter", ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                "-readFilter", ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()});

        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);

        // get and verify the list of filters enabled on the command line (not including defaults)
        List<ReadFilter> orderedFilters = readFilterPlugin.getAllInstances();
        Assert.assertEquals(orderedFilters.size(), 3);
        Assert.assertEquals(orderedFilters.get(0).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());
        Assert.assertEquals(orderedFilters.get(1).getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName());
        Assert.assertEquals(orderedFilters.get(2).getClass().getSimpleName(),
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName());

        // Now get the final merged read filter and verify the execution order. We need to ensure that
        // getMergedReadFilter creates a composite filter that honors the filter test execution
        // order rules (tool defaults first, in order, followed by command line-specified, in order
        // listed), so validate by reaching inside the filter and navigating down the tree.
        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        String expectedOrder[] = {
                WellformedReadFilter.class.getSimpleName(),
                ReadFilterLibrary.HAS_READ_GROUP.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()
        };

        int count = ReadFilterUnitTest.verifyAndFilterOrder(rf, expectedOrder);
        Assert.assertEquals(count, 6);
    }

    @Test
    public void testPreserveToolDefaultFilterOrder() {
        List<ReadFilter> orderedDefaults = new ArrayList<>();
        orderedDefaults.add(new WellformedReadFilter());
        orderedDefaults.add(ReadFilterLibrary.MAPPED);
        orderedDefaults.add(ReadFilterLibrary.HAS_READ_GROUP);
        orderedDefaults.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        orderedDefaults.add(ReadFilterLibrary.PAIRED);
        orderedDefaults.add(ReadFilterLibrary.NONZERO_FRAGMENT_LENGTH_READ_FILTER);
        orderedDefaults.add(ReadFilterLibrary.FIRST_OF_PAIR);
        orderedDefaults.add(ReadFilterLibrary.PROPERLY_PAIRED);
        orderedDefaults.add(ReadFilterLibrary.NOT_DUPLICATE);
        orderedDefaults.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        orderedDefaults.add(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);

        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(orderedDefaults)));
        clp.parseArguments(System.out, new String[] {
                //disable one just to mix things up
                "-disableReadFilter", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "-readFilter", ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                "-readFilter", ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()});

        // Now get the final merged read filter and verify the execution order. We need to ensure that
        // getMergedReadFilter creates a composite filter that honors the filter test execution
        // order rules (tool defaults first, in order, followed by command line-specified, in order
        // listed). So reach inside the filter and navigate down the tree.
        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        String expectedOrder[] = {
                WellformedReadFilter.class.getSimpleName(),
                ReadFilterLibrary.HAS_READ_GROUP.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO.getClass().getSimpleName(),
                ReadFilterLibrary.PAIRED.getClass().getSimpleName(),
                ReadFilterLibrary.NONZERO_FRAGMENT_LENGTH_READ_FILTER.getClass().getSimpleName(),
                ReadFilterLibrary.FIRST_OF_PAIR.getClass().getSimpleName(),
                ReadFilterLibrary.PROPERLY_PAIRED.getClass().getSimpleName(),
                ReadFilterLibrary.NOT_DUPLICATE.getClass().getSimpleName(),
                ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT.getClass().getSimpleName(),
                ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT.getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()
        };

        int count = ReadFilterUnitTest.verifyAndFilterOrder(rf, expectedOrder);
        Assert.assertEquals(count, 12);
    }

    @Test
    public void testReadLengthFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)));
        String[] args = {
                "--readFilter", ReadLengthReadFilter.class.getSimpleName(),
                "--minReadLength", "10",
                "--maxReadLength", "20"
        };
        clp.parseArguments(System.out, args);
        ReadFilter rf = instantiateFilter(clp, header);

        read.setBases(new byte[5]);
        Assert.assertFalse(rf.test(read));
        read.setBases(new byte[25]);
        Assert.assertFalse(rf.test(read));

        read.setBases(new byte[15]);
        Assert.assertTrue(rf.test(read));
    }

    final private static String readgroupName = "ReadGroup1";

    @DataProvider(name="testDefaultFilters")
    public Object[][] defaultFilters() {
        return new Object[][]{
                // ReadGroupReadFilter has a required "keepReadGroup" arg; provide it
                {new String[]{"--keepReadGroup", readgroupName}},

                // ReadGroupReadFilter has a required "keepReadGroup" arg; provide it
                // *and* specify it on the command line
                {new String[]{
                        "--readFilter", "ReadGroupReadFilter",
                        "--keepReadGroup", readgroupName}
                }
        };
    }

    @Test(dataProvider = "testDefaultFilters")
    public void testToolHasDefaultRequiredArgsPositive(final String[] args) {
        CommandLineParser clp = new CommandLineArgumentParser(
            new Object(),
            Collections.singletonList(
                    new GATKReadFilterPluginDescriptor(Collections.singletonList(new ReadGroupReadFilter()))
            )
        );
        clp.parseArguments(System.out, args);

        SAMFileHeader samHeader = createHeaderWithReadGroups();
        ReadFilter rf = instantiateFilter(clp, samHeader);
        final GATKRead read = simpleGoodRead(samHeader);

        // make sure the test read actually has this read group
        Assert.assertEquals(read.getReadGroup(), readgroupName);
        Assert.assertTrue(rf.test(read));
    }

    @Test(expectedExceptions = CommandLineException.MissingArgument.class)
    public void testToolHasDefaultRequiredArgsNegative() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(
                        new GATKReadFilterPluginDescriptor(Collections.singletonList(new ReadGroupReadFilter()))
                )
        );

        // ReadGroupReadFilter has a required "keepReadGroup" arg; don't provide it and fail
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    private ReadFilter instantiateFilter(
            final CommandLineParser clp,
            final SAMFileHeader header)
    {
        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return readFilterPlugin.getMergedReadFilter(header);
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
