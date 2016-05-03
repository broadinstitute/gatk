package org.broadinstitute.hellbender.cmdline;

import htsjdk.samtools.util.CollectionUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.*;

public final class CommandLineParserTest {
    enum FrobnicationFlavor {
        FOO, BAR, BAZ
    }

    @CommandLineProgramProperties(
            summary = "Usage: frobnicate [arguments] input-file output-file\n\nRead input-file, frobnicate it, and write frobnicated results to output-file\n",
            oneLineSummary = "Read input-file, frobnicate it, and write frobnicated results to output-file",
            programGroup = QCProgramGroup.class
    )
    class FrobnicateArguments {
        @ArgumentCollection
        SpecialArgumentsCollection specialArgs = new SpecialArgumentsCollection();

        @PositionalArguments(minElements=2, maxElements=2)
        public List<File> positionalArguments = new ArrayList<>();

        @Argument(shortName="T", doc="Frobnication threshold setting.")
        public Integer FROBNICATION_THRESHOLD = 20;

        @Argument
        public FrobnicationFlavor FROBNICATION_FLAVOR;

        @Argument(doc="Allowed shmiggle types.", optional = false)
        public List<String> SHMIGGLE_TYPE = new ArrayList<>();

        @Argument
        public Boolean TRUTHINESS = false;
    }

    @CommandLineProgramProperties(
            summary = "Usage: framistat [arguments]\n\nCompute the plebnick of the freebozzle.\n",
            oneLineSummary = "ompute the plebnick of the freebozzle",
            programGroup = QCProgramGroup.class
    )
    class ArgumentsWithoutPositional {
        public static final int DEFAULT_FROBNICATION_THRESHOLD = 20;
        @Argument(shortName="T", doc="Frobnication threshold setting.")
        public Integer FROBNICATION_THRESHOLD = DEFAULT_FROBNICATION_THRESHOLD;

        @Argument
        public FrobnicationFlavor FROBNICATION_FLAVOR;

        @Argument(doc="Allowed shmiggle types.", optional = false)
        public List<String> SHMIGGLE_TYPE = new ArrayList<>();

        @Argument
        public Boolean TRUTHINESS;
    }

    class MutexArguments {
        @Argument(mutex={"M", "N", "Y", "Z"})
        public String A;
        @Argument(mutex={"M", "N", "Y", "Z"})
        public String B;
        @Argument(mutex={"A", "B", "Y", "Z"})
        public String M;
        @Argument(mutex={"A", "B", "Y", "Z"})
        public String N;
        @Argument(mutex={"A", "B", "M", "N"})
        public String Y;
        @Argument(mutex={"A", "B", "M", "N"})
        public String Z;
        
    }

    @CommandLineProgramProperties(
            summary = "[oscillation_frequency]\n\nResets oscillation frequency.\n",
            oneLineSummary = "Reset oscillation frequency.",
            programGroup = QCProgramGroup.class
    )
    public class RequiredOnlyArguments {
        @Argument(doc="Oscillation frequency.", optional = false)
        public String OSCILLATION_FREQUENCY;
    }

    @Test
    public void testRequiredOnlyUsage() {
        final RequiredOnlyArguments nr = new RequiredOnlyArguments();
        final CommandLineParser clp = new CommandLineParser(nr);
        final String out = BaseTest.captureStderr(() -> clp.usage(System.err, false)); // without common args
        final int reqIndex = out.indexOf("Required Arguments:");
        Assert.assertTrue(reqIndex > 0);
        Assert.assertTrue(out.indexOf("Optional Arguments:", reqIndex) < 0);
    }

    class AbbreviatableArgument{
        public static final String ARGUMENT_NAME = "longNameArgument";
        @Argument(fullName= ARGUMENT_NAME)
        public boolean longNameArgument;
    }

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testAbbreviationsAreRejected() {
        final AbbreviatableArgument abrv = new AbbreviatableArgument();
        final CommandLineParser clp = new CommandLineParser(abrv);
        //argument name is valid when it isn't abbreviated
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--" + AbbreviatableArgument.ARGUMENT_NAME}));

        //should throw when the abbreviated name is used
        clp.parseArguments(System.err, new String[]{"--" + AbbreviatableArgument.ARGUMENT_NAME.substring(0,5)});
    }

    @CommandLineProgramProperties(
            summary = "[oscillation_frequency]\n\nRecalibrates overthruster oscillation. \n",
            oneLineSummary = "Recalibrates the overthruster.",
            programGroup = QCProgramGroup.class
    )
    public class OptionalOnlyArguments {
        @Argument(doc="Oscillation frequency.", optional = true)
        public String OSCILLATION_FREQUENCY = "20";
    }

    @Test
    public void testOptionalOnlyUsage() {
        final OptionalOnlyArguments oo = new OptionalOnlyArguments();
        final CommandLineParser clp = new CommandLineParser(oo);
        final String out = BaseTest.captureStderr(() -> clp.usage(System.err, false)); // without common args
        final int reqIndex = out.indexOf("Required Arguments:");
        Assert.assertTrue(reqIndex < 0);
        Assert.assertTrue(out.indexOf("Optional Arguments:", reqIndex) > 0);
    }

    /**
     * Validate the text emitted by a call to usage by ensuring that required arguments are
     * emitted before optional ones.
     */
    private void validateRequiredOptionalUsage(final CommandLineParser clp, final boolean withDefault) {
        final String out = BaseTest.captureStderr(() -> clp.usage(System.err, withDefault)); // with common args
        // Required arguments should appear before optional ones
        final int reqIndex = out.indexOf("Required Arguments:");
        Assert.assertTrue(reqIndex > 0);
        Assert.assertTrue(out.indexOf("Optional Arguments:", reqIndex) > 0);
    }

    @Test
    public void testRequiredOptionalWithDefaultUsage() {
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        validateRequiredOptionalUsage(clp, true); // with common args
    }

    @Test
    public void testRequiredOptionalWithoutDefaultUsage() {
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        validateRequiredOptionalUsage(clp, false); // without common args
    }

    @Test
    public void testWithoutPositionalWithDefaultUsage() {
        final ArgumentsWithoutPositional fo = new ArgumentsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        validateRequiredOptionalUsage(clp, true); // with common args
    }

    @Test
    public void testWithoutPositionalWithoutDefaultUsage() {
        final ArgumentsWithoutPositional fo = new ArgumentsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        validateRequiredOptionalUsage(clp, false); // without commoon args
    }

    @Test
    public void testPositive() {
        final String[] args = {
                "-T","17",
                "-FROBNICATION_FLAVOR","BAR",
                "-TRUTHINESS",
                "-SHMIGGLE_TYPE","shmiggle1",
                "-SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(fo.positionalArguments.size(), 2);
        final File[] expectedPositionalArguments = { new File("positional1"), new File("positional2")};
        Assert.assertEquals(fo.positionalArguments.toArray(), expectedPositionalArguments);
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 17);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 2);
        final String[] expectedShmiggleTypes = {"shmiggle1", "shmiggle2"};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertTrue(fo.TRUTHINESS);
    }

    @Test
    public void testGetCommandLine() {
        final String[] args = {
                "-T","17",
                "-FROBNICATION_FLAVOR","BAR",
                "-TRUTHINESS",
                "-SHMIGGLE_TYPE","shmiggle1",
                "-SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(clp.getCommandLine(),
                "org.broadinstitute.hellbender.cmdline.CommandLineParserTest$FrobnicateArguments  " +
                        "positional1 positional2 --FROBNICATION_THRESHOLD 17 --FROBNICATION_FLAVOR BAR " +
                        "--SHMIGGLE_TYPE shmiggle1 --SHMIGGLE_TYPE shmiggle2 --TRUTHINESS true  --help false " +
                        "--version false");
    }

    private static class WithSensitiveValues {

        @Argument(sensitive = true)
        public String secretValue;

        @Argument
        public String openValue;
    }

    @Test
    public void testGetCommandLineWithSensitiveArgument(){
        final String supersecret = "supersecret";
        final String unclassified = "unclassified";
        final String[] args = {
                "--secretValue", supersecret,
                "--openValue", unclassified
        };
        final WithSensitiveValues sv = new WithSensitiveValues();
        final CommandLineParser clp = new CommandLineParser(sv);
        Assert.assertTrue(clp.parseArguments(System.err, args));

        final String commandLine = clp.getCommandLine();

        Assert.assertTrue(commandLine.contains(unclassified));
        Assert.assertFalse(commandLine.contains(supersecret));

        Assert.assertEquals(sv.openValue, unclassified);
        Assert.assertEquals(sv.secretValue, supersecret);
    }

    @Test
    public void testDefault() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 20);
    }

    @Test(expectedExceptions = UserException.MissingArgument.class)
    public void testMissingRequiredArgument() {
        final String[] args = {
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    class CollectionRequired{
        @Argument(optional = false)
        List<Integer> ints;
    }

    @Test(expectedExceptions = UserException.MissingArgument.class)
    public void testMissingRequiredCollectionArgument(){
        final String[] args = {};
        final CollectionRequired cr = new CollectionRequired();
        final CommandLineParser clp = new CommandLineParser(cr);
        clp.parseArguments(System.err, args);
    }

    @Test( expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadValue() {
        final String[] args = {
                "--FROBNICATION_THRESHOLD","ABC",
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadEnumValue() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","HiMom",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.MissingArgument.class)
    public void testNotEnoughOfListArgument() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testTooManyPositional() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "positional1",
                "positional2",
                "positional3",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.MissingArgument.class)
    public void testNotEnoughPositional() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test( expectedExceptions = UserException.CommandLineException.class)
    public void testUnexpectedPositional() {
        final String[] args = {
                "--T","17",
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "positional"
        };
        final ArgumentsWithoutPositional fo = new ArgumentsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testArgumentUseClash() {
        final String[] args = {
                "--FROBNICATION_FLAVOR", "BAR",
                "--FROBNICATION_FLAVOR", "BAZ",
                "--SHMIGGLE_TYPE", "shmiggle1",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }

    @Test
    public void testArgumentsFile() throws Exception {
        final File argumentsFile = BaseTest.createTempFile("clp.", ".arguments");
        try (final PrintWriter writer = new PrintWriter(argumentsFile)) {
            writer.println("-T 18");
            writer.println("--TRUTHINESS");
            writer.println("--SHMIGGLE_TYPE shmiggle0");
            writer.println("--" + SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME + " " + argumentsFile.getPath());
            //writer.println("--STRANGE_ARGUMENT shmiggle0");
        }
        final String[] args = {
                "--"+SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME, argumentsFile.getPath(),
                // Multiple arguments files are allowed
                "--"+SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME, argumentsFile.getPath(),
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS",
                "--SHMIGGLE_TYPE","shmiggle0",
                "--SHMIGGLE_TYPE","shmiggle1",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(fo.positionalArguments.size(), 2);
        final File[] expectedPositionalArguments = { new File("positional1"), new File("positional2")};
        Assert.assertEquals(fo.positionalArguments.toArray(), expectedPositionalArguments);
        Assert.assertEquals(fo.FROBNICATION_THRESHOLD.intValue(), 18);
        Assert.assertEquals(fo.FROBNICATION_FLAVOR, FrobnicationFlavor.BAR);
        Assert.assertEquals(fo.SHMIGGLE_TYPE.size(), 3);
        final String[] expectedShmiggleTypes = {"shmiggle0", "shmiggle0", "shmiggle1"};
        Assert.assertEquals(fo.SHMIGGLE_TYPE.toArray(), expectedShmiggleTypes);
        Assert.assertTrue(fo.TRUTHINESS);
    }


    /**
     * In an arguments file, should not be allowed to override an argument set on the command line
     * @throws Exception
     */
    @Test( expectedExceptions = UserException.CommandLineException.class)
    public void testArgumentsFileWithDisallowedOverride() throws Exception {
        final File argumentsFile = BaseTest.createTempFile("clp.", ".arguments");
        try (final PrintWriter writer = new PrintWriter(argumentsFile)) {
            writer.println("--T 18");
        }
        final String[] args = {
                "--T","17",
                "--"+SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME ,argumentsFile.getPath()
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.parseArguments(System.err, args);
    }
    
    @DataProvider(name="failingMutexScenarios")
    public Object[][] failingMutexScenarios() {
        return new Object[][] {
                { "no args", new String[0], false },
                { "1 of group required", new String[] {"-A","1"}, false },
                { "mutex", new String[]  {"-A","1", "-Y","3"}, false },
                { "mega mutex", new String[]  {"-A","1", "-B","2", "-Y","3", "-Z","1", "-M","2", "-N","3"}, false }
        };
    }

    @Test
    public void passingMutexCheck(){
        final MutexArguments o = new MutexArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"-A", "1", "-B", "2"}));
    }

    @Test(dataProvider="failingMutexScenarios", expectedExceptions = UserException.CommandLineException.class)
    public void testFailingMutex(final String testName, final String[] args, final boolean expected) {
        final MutexArguments o = new MutexArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        clp.parseArguments(System.err, args);
    }

    class UninitializedCollectionArguments {
        @Argument
        public List<String> LIST;
        @Argument
        public ArrayList<String> ARRAY_LIST;
        @Argument
        public HashSet<String> HASH_SET;
        @PositionalArguments
        public Collection<File> COLLECTION;

    }

    @Test
    public void testUninitializedCollections() {
        final UninitializedCollectionArguments o = new UninitializedCollectionArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"--LIST","L1", "--LIST","L2", "--ARRAY_LIST","S1", "--HASH_SET","HS1", "P1", "P2"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.LIST.size(), 2);
        Assert.assertEquals(o.ARRAY_LIST.size(), 1);
        Assert.assertEquals(o.HASH_SET.size(), 1);
        Assert.assertEquals(o.COLLECTION.size(), 2);
    }

    class UninitializedCollectionThatCannotBeAutoInitializedArguments {
        @Argument
        public Set<String> SET;
    }

    @Test(expectedExceptions = GATKException.CommandLineParserInternalException.class)
    public void testCollectionThatCannotBeAutoInitialized() {
        final UninitializedCollectionThatCannotBeAutoInitializedArguments o = new UninitializedCollectionThatCannotBeAutoInitializedArguments();
        new CommandLineParser(o);
    }

    class CollectionWithDefaultValuesArguments {
        @Argument
        public List<String> LIST = CollectionUtil.makeList("foo", "bar");
    }

    @Test
    public void testClearDefaultValuesFromListArgument() {
        final CollectionWithDefaultValuesArguments o = new CollectionWithDefaultValuesArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"--LIST","null"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.LIST.size(), 0);
    }

    @Test
    public void testClearDefaultValuesFromListArgumentAndAddNew() {
        final CollectionWithDefaultValuesArguments o = new CollectionWithDefaultValuesArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"--LIST","null", "--LIST","baz", "--LIST","frob"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.LIST, CollectionUtil.makeList("baz", "frob"));
    }

    @Test
    public void testDefaultValuesListArgument() {
        final CollectionWithDefaultValuesArguments o = new CollectionWithDefaultValuesArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"--LIST","baz", "--LIST","frob"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.LIST, CollectionUtil.makeList("foo", "bar", "baz", "frob"));
    }


    @Test
       public void testFlagNoArgument(){
        final BooleanFlags o = new BooleanFlags();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--flag1"}));
        Assert.assertTrue(o.flag1);
    }

    @Test
    public void testFlagsWithArguments(){
        final BooleanFlags o = new BooleanFlags();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--flag1", "false", "--flag2", "false"}));
        Assert.assertFalse(o.flag1);
        Assert.assertFalse(o.flag2);
    }

    class ArgsCollection {
        @Argument(fullName = "arg1")
        public int Arg1;
    }

    class ArgsCollectionHaver{

        public ArgsCollectionHaver(){}

        @ArgumentCollection
        public ArgsCollection default_args = new ArgsCollection();

        @Argument(fullName = "somenumber",shortName = "n")
        public int someNumber = 0;
    }

    @Test
    public void testArgumentCollection(){
        final ArgsCollectionHaver o = new ArgsCollectionHaver();
        final CommandLineParser clp = new CommandLineParser(o);

        String[] args = {"--arg1", "42", "--somenumber", "12"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.someNumber, 12);
        Assert.assertEquals(o.default_args.Arg1, 42);

    }

    class BooleanFlags{
        @Argument
        public Boolean flag1 = false;

        @Argument
        public boolean flag2 = true;

        @Argument
        public boolean flag3 = false;

        @Argument(mutex="flag1")
        public boolean antiflag1 = false;

        @ArgumentCollection
        SpecialArgumentsCollection special = new SpecialArgumentsCollection();
    }

    @Test
    public void testCombinationOfFlags(){
        final BooleanFlags o = new BooleanFlags();
        final CommandLineParser clp = new CommandLineParser(o);

        clp.parseArguments(System.err, new String[]{"--flag1", "false", "--flag2"});
        Assert.assertFalse(o.flag1);
        Assert.assertTrue(o.flag2);
        Assert.assertFalse(o.flag3);
    }

    class WithBadField{
        @Argument
        @ArgumentCollection
        public boolean badfield;
    }

    @Test(expectedExceptions = GATKException.CommandLineParserInternalException.class)
    public void testBadFieldCausesException(){
        WithBadField o = new WithBadField();
        final CommandLineParser clp = new CommandLineParser(o);
    }

    class PrivateArgument{
        @Argument
        private boolean privateArgument = false;

        @Argument(optional = true)
        private List<Integer> privateCollection = new ArrayList<>();

        @ArgumentCollection
        private BooleanFlags booleanFlags= new BooleanFlags();

        @PositionalArguments()
        List<Integer> positionals = new ArrayList<>();
    }

    @Test
    public void testFlagWithPositionalFollowing(){
        PrivateArgument o = new PrivateArgument();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--flag1","1","2" }));
        Assert.assertTrue(o.booleanFlags.flag1);
        Assert.assertEquals(o.positionals, Arrays.asList(1, 2));
    }

    @Test
    public void testPrivateArgument(){
        PrivateArgument o = new PrivateArgument();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--privateArgument",
                "--privateCollection", "1", "--privateCollection", "2", "--flag1"}));
        Assert.assertTrue(o.privateArgument);
        Assert.assertEquals(o.privateCollection, Arrays.asList(1,2));
        Assert.assertTrue(o.booleanFlags.flag1);
    }

    /**
     * Test that the special flags --help and --version are handled correctly
     */
    @Test
    public void testSpecialFlags(){
        final BooleanFlags o = new BooleanFlags();
        final CommandLineParser clp = new CommandLineParser(o);

        String[] helpArgs = new String[]{"--"+SpecialArgumentsCollection.HELP_FULLNAME};

        String out = BaseTest.captureStderr( () -> {
            Assert.assertFalse(clp.parseArguments(System.err, helpArgs));
            });
        BaseTest.assertContains(out, "Usage:");

        final String[] versionArgs = {"--" + SpecialArgumentsCollection.VERSION_FULLNAME};
        out = BaseTest.captureStderr(() -> {
                    Assert.assertFalse(clp.parseArguments(System.err, versionArgs));
            });
        BaseTest.assertContains(out,"Version:");

        Assert.assertTrue(clp.parseArguments(System.err,  new String[]{"--help", "false"}));
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--version","false"}));

        Assert.assertFalse(clp.parseArguments(System.err, new String[]{"--help", "true"}));
        Assert.assertFalse(clp.parseArguments(System.err, new String[]{"--version", "true"}));
    }

    // test ArgumentCollection depends on
    class ArgTarget {
        @Argument(fullName="argTarget", optional = true)
        String argTarget;
    }
    class ArgAndValueTarget {
        @Argument(fullName="argAndValueTarget", optional = true)
        String argAndValueTarget;
    }

    class TestArgs {
        @Argument(fullName="master", optional = true)
        int rf;

        @ArgumentCollection(dependsOnArgument = "master")
        ArgTarget t1 = new ArgTarget();

        @ArgumentCollection(dependsOnArgument = "master", dependsOnValue="27")
        ArgAndValueTarget t2 = new ArgAndValueTarget();
    }

    @Test
    public void testDependsOnArgument() {
        TestArgs tArgs = new TestArgs();
        String args[] = {
                "--master", "35",
                "-argTarget", "foo"
        };
        final CommandLineParser clp = new CommandLineParser(tArgs);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testFailsDependsOnArgument() {
        TestArgs tArgs = new TestArgs();
        String args[] = {
                //missing --master
                "--argTarget", "foo"
        };
        final CommandLineParser clp = new CommandLineParser(tArgs);
        clp.parseArguments(System.err, args);
    }

    @Test
    public void testDependsOnArgumentAndValue() {
        TestArgs tArgs = new TestArgs();
        String args[] = {
                "--master", "27",
                "--argAndValueTarget", "foo"
        };
        final CommandLineParser clp = new CommandLineParser(tArgs);
        clp.parseArguments(System.err, args);
    }

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testFailsDependsOnArgumentValue() {
        TestArgs tArgs = new TestArgs();
        String args[] = {
                "--master", "0",
                "--argAndValueTarget", "foo" // depends on --master 27
        };
        final CommandLineParser clp = new CommandLineParser(tArgs);
        clp.parseArguments(System.err, args);
    }

    /***************************************************************************************
     * Start of tests and helper classes for CommandLineParser.gatherArgumentValuesOfType()
     ***************************************************************************************/

    /**
     * Classes and argument collections for use with CommandLineParser.gatherArgumentValuesOfType() tests below.
     *
     * Structured to ensure that we test support for:
     *
     * -distinguishing between arguments of the target type, and arguments not of the target type
     * -distinguishing between annotated and unannotated fields of the target type
     * -gathering arguments that are a subtype of the target type
     * -gathering multi-valued arguments of the target type within Collection types
     * -gathering arguments of the target type that are not specified on the command line
     * -gathering arguments of the target type from superclasses of our tool
     * -gathering arguments of the target type from argument collections
     * -gathering arguments when the target type is itself a parameterized type (eg., FeatureInput<VariantContext>)
     */

    private static class GatherArgumentValuesTestSourceParent {
        @Argument(fullName = "parentSuperTypeTarget", shortName = "parentSuperTypeTarget", doc = "")
        private GatherArgumentValuesTargetSuperType parentSuperTypeTarget;

        @Argument(fullName = "parentSubTypeTarget", shortName = "parentSubTypeTarget", doc = "")
        private GatherArgumentValuesTargetSubType parentSubTypeTarget;

        @Argument(fullName = "parentListSuperTypeTarget", shortName = "parentListSuperTypeTarget", doc = "")
        private List<GatherArgumentValuesTargetSuperType> parentListSuperTypeTarget;

        @Argument(fullName = "parentListSubTypeTarget", shortName = "parentListSubTypeTarget", doc = "")
        private List<GatherArgumentValuesTargetSubType> parentListSubTypeTarget;

        @Argument(fullName = "uninitializedParentTarget", shortName = "uninitializedParentTarget", optional = true, doc = "")
        private GatherArgumentValuesTargetSuperType uninitializedParentTarget;

        @Argument(fullName = "parentNonTargetArgument", shortName = "parentNonTargetArgument", doc = "")
        private int parentNonTargetArgument;

        private GatherArgumentValuesTargetSuperType parentUnannotatedTarget;

        @ArgumentCollection
        private GatherArgumentValuesTestSourceParentCollection parentCollection = new GatherArgumentValuesTestSourceParentCollection();
    }

    private static class GatherArgumentValuesTestSourceChild extends GatherArgumentValuesTestSourceParent {
        @Argument(fullName = "childSuperTypeTarget", shortName = "childSuperTypeTarget", doc = "")
        private GatherArgumentValuesTargetSuperType childSuperTypeTarget;

        @Argument(fullName = "childSubTypeTarget", shortName = "childSubTypeTarget", doc = "")
        private GatherArgumentValuesTargetSubType childSubTypeTarget;

        @Argument(fullName = "childListSuperTypeTarget", shortName = "childListSuperTypeTarget", doc = "")
        private List<GatherArgumentValuesTargetSuperType> childListSuperTypeTarget;

        @Argument(fullName = "childListSubTypeTarget", shortName = "childListSubTypeTarget", doc = "")
        private List<GatherArgumentValuesTargetSubType> childListSubTypeTarget;

        @Argument(fullName = "uninitializedChildTarget", shortName = "uninitializedChildTarget", optional = true, doc = "")
        private GatherArgumentValuesTargetSuperType uninitializedChildTarget;

        @Argument(fullName = "uninitializedChildListTarget", shortName = "uninitializedChildListTarget", optional = true, doc = "")
        private List<GatherArgumentValuesTargetSuperType> uninitializedChildListTarget;

        @Argument(fullName = "childNonTargetArgument", shortName = "childNonTargetArgument", doc = "")
        private int childNonTargetArgument;

        @Argument(fullName = "childNonTargetListArgument", shortName = "childNonTargetListArgument", doc = "")
        private List<Integer> childNonTargetListArgument;

        private GatherArgumentValuesTargetSuperType childUnannotatedTarget;

        @ArgumentCollection
        private GatherArgumentValuesTestSourceChildCollection childCollection = new GatherArgumentValuesTestSourceChildCollection();
    }

    private static class GatherArgumentValuesTestSourceParentCollection implements ArgumentCollectionDefinition {
        private static final long serialVersionUID = 1L;

        @Argument(fullName = "parentCollectionSuperTypeTarget", shortName = "parentCollectionSuperTypeTarget", doc = "")
        private GatherArgumentValuesTargetSuperType parentCollectionSuperTypeTarget;

        @Argument(fullName = "parentCollectionSubTypeTarget", shortName = "parentCollectionSubTypeTarget", doc = "")
        private GatherArgumentValuesTargetSubType parentCollectionSubTypeTarget;

        @Argument(fullName = "uninitializedParentCollectionTarget", shortName = "uninitializedParentCollectionTarget", optional = true, doc = "")
        private GatherArgumentValuesTargetSuperType uninitializedParentCollectionTarget;

        @Argument(fullName = "parentCollectionNonTargetArgument", shortName = "parentCollectionNonTargetArgument", doc = "")
        private int parentCollectionNonTargetArgument;

        private GatherArgumentValuesTargetSuperType parentCollectionUnannotatedTarget;
    }

    private static class GatherArgumentValuesTestSourceChildCollection implements ArgumentCollectionDefinition {
        private static final long serialVersionUID = 1L;

        @Argument(fullName = "childCollectionSuperTypeTarget", shortName = "childCollectionSuperTypeTarget", doc = "")
        private GatherArgumentValuesTargetSuperType childCollectionSuperTypeTarget;

        @Argument(fullName = "childCollectionSubTypeTarget", shortName = "childCollectionSubTypeTarget", doc = "")
        private GatherArgumentValuesTargetSubType childCollectionSubTypeTarget;

        @Argument(fullName = "childCollectionListSuperTypeTarget", shortName = "childCollectionListSuperTypeTarget", doc = "")
        private List<GatherArgumentValuesTargetSuperType> childCollectionListSuperTypeTarget;

        @Argument(fullName = "uninitializedChildCollectionTarget", shortName = "uninitializedChildCollectionTarget", optional = true, doc = "")
        private GatherArgumentValuesTargetSuperType uninitializedChildCollectionTarget;

        @Argument(fullName = "childCollectionNonTargetArgument", shortName = "childCollectionNonTargetArgument", doc = "")
        private int childCollectionNonTargetArgument;

        private GatherArgumentValuesTargetSuperType childCollectionUnannotatedTarget;
    }

    /**
     * Our tests will search for argument values of this type, subtypes of this type, and Collections of
     * this type or its subtypes. Has a String constructor so that the argument parsing system can correctly
     * initialize it.
     */
    private static class GatherArgumentValuesTargetSuperType {
        private String value;

        public GatherArgumentValuesTargetSuperType( String s ) {
            value = s;
        }

        public String getValue() {
            return value;
        }
    }

    private static class GatherArgumentValuesTargetSubType extends GatherArgumentValuesTargetSuperType {
        public GatherArgumentValuesTargetSubType( String s ) {
            super(s);
        }
    }

    @DataProvider(name = "gatherArgumentValuesOfTypeDataProvider")
    public Object[][] gatherArgumentValuesOfTypeDataProvider() {
        // Non-Collection arguments of the target type
        final List<String> targetScalarArguments = Arrays.asList("childSuperTypeTarget", "childSubTypeTarget",
                                                                 "parentSuperTypeTarget", "parentSubTypeTarget",
                                                                 "childCollectionSuperTypeTarget", "childCollectionSubTypeTarget",
                                                                 "parentCollectionSuperTypeTarget", "parentCollectionSubTypeTarget");
        // Collection arguments of the target type
        final List<String> targetListArguments = Arrays.asList("childListSuperTypeTarget", "childListSubTypeTarget",
                                                               "parentListSuperTypeTarget", "parentListSubTypeTarget",
                                                               "childCollectionListSuperTypeTarget");
        // Arguments of the target type that we won't specify on our command line
        final List<String> uninitializedTargetArguments = Arrays.asList("uninitializedChildTarget", "uninitializedChildListTarget",
                                                                        "uninitializedParentTarget", "uninitializedChildCollectionTarget",
                                                                        "uninitializedParentCollectionTarget");
        // Arguments not of the target type
        final List<String> nonTargetArguments = Arrays.asList("childNonTargetArgument", "parentNonTargetArgument",
                                                              "childCollectionNonTargetArgument", "parentCollectionNonTargetArgument",
                                                              "childNonTargetListArgument");

        List<String> commandLineArguments = new ArrayList<>();
        List<Pair<String, String>> sortedExpectedGatheredValues = new ArrayList<>();

        for ( String targetScalarArgument : targetScalarArguments ) {
            final String argumentValue = targetScalarArgument + "Value";

            commandLineArguments.add("--" + targetScalarArgument);
            commandLineArguments.add(argumentValue);
            sortedExpectedGatheredValues.add(Pair.of(targetScalarArgument, argumentValue));
        }

        // Give each list argument multiple values
        for ( String targetListArgument : targetListArguments ) {
            for ( int argumentNum = 1; argumentNum <= 3; ++argumentNum ) {
                final String argumentValue = targetListArgument + "Value" + argumentNum;

                commandLineArguments.add("--" + targetListArgument);
                commandLineArguments.add(argumentValue);
                sortedExpectedGatheredValues.add(Pair.of(targetListArgument, argumentValue));
            }
        }

        // Make sure the uninitialized args of the target type not included on the command line are
        // represented in the expected output
        for ( String uninitializedTargetArgument : uninitializedTargetArguments ) {
            sortedExpectedGatheredValues.add(Pair.of(uninitializedTargetArgument, null));
        }

        // The non-target args are all of type int, so give them an arbitrary int value on the command line.
        // These should not be gathered at all, so are not added to the expected output.
        for ( String nonTargetArgument : nonTargetArguments ) {
            commandLineArguments.add("--" + nonTargetArgument);
            commandLineArguments.add("1");
        }

        Collections.sort(sortedExpectedGatheredValues);

        return new Object[][] {{
            commandLineArguments, sortedExpectedGatheredValues
        }};
    }

    @Test(dataProvider = "gatherArgumentValuesOfTypeDataProvider")
    public void testGatherArgumentValuesOfType( final List<String> commandLineArguments, final List<Pair<String, String>> sortedExpectedGatheredValues ) {
        GatherArgumentValuesTestSourceChild argumentSource = new GatherArgumentValuesTestSourceChild();

        // Parse the command line, and inject values into our test instance
        CommandLineParser clp = new CommandLineParser(argumentSource);
        clp.parseArguments(System.err, commandLineArguments.toArray(new String[commandLineArguments.size()]));

        // Gather all argument values of type GatherArgumentValuesTargetSuperType (or Collection<GatherArgumentValuesTargetSuperType>),
        // including subtypes.
        List<Pair<Field, GatherArgumentValuesTargetSuperType>> gatheredArguments =
                CommandLineParser.gatherArgumentValuesOfType(GatherArgumentValuesTargetSuperType.class, argumentSource);

        // Make sure we gathered the expected number of argument values
        Assert.assertEquals(gatheredArguments.size(), sortedExpectedGatheredValues.size(), "Gathered the wrong number of arguments");

        // Make sure actual gathered argument values match expected values
        List<Pair<String, String>> sortedActualGatheredArgumentValues = new ArrayList<>();
        for ( Pair<Field, GatherArgumentValuesTargetSuperType> gatheredArgument : gatheredArguments ) {
            Assert.assertNotNull(gatheredArgument.getKey().getAnnotation(Argument.class), "Gathered argument is not annotated with an @Argument annotation");

            String argumentName = gatheredArgument.getKey().getAnnotation(Argument.class).fullName();
            GatherArgumentValuesTargetSuperType argumentValue = gatheredArgument.getValue();

            sortedActualGatheredArgumentValues.add(Pair.of(argumentName, argumentValue != null ? argumentValue.getValue() : null));
        }
        Collections.sort(sortedActualGatheredArgumentValues);

        Assert.assertEquals(sortedActualGatheredArgumentValues, sortedExpectedGatheredValues,
                            "One or more gathered argument values not correct");

    }

    /**
     * Nonsensical parameterized class, just to ensure that CommandLineParser.gatherArgumentValuesOfType()
     * can gather argument values of a generic type
     *
     * @param <T> meaningless type parameter
     */
    private static class GatherArgumentValuesParameterizedTargetType<T> {
        private String value;
        private T foo;

        public GatherArgumentValuesParameterizedTargetType( String s ) {
            value = s;
            foo = null;
        }

        public String getValue() {
            return value;
        }
    }

    private static class GatherArgumentValuesParameterizedTypeSource {
        @Argument(fullName = "parameterizedTypeArgument", shortName = "parameterizedTypeArgument", doc = "")
        private GatherArgumentValuesParameterizedTargetType<Integer> parameterizedTypeArgument;

        @Argument(fullName = "parameterizedTypeListArgument", shortName = "parameterizedTypeListArgument", doc = "")
        private List<GatherArgumentValuesParameterizedTargetType<Integer>> parameterizedTypeListArgument;
    }

    @Test
    @SuppressWarnings("rawtypes")
    public void testGatherArgumentValuesOfTypeWithParameterizedType() {
        GatherArgumentValuesParameterizedTypeSource argumentSource = new GatherArgumentValuesParameterizedTypeSource();

        // Parse the command line, and inject values into our test instance
        CommandLineParser clp = new CommandLineParser(argumentSource);
        clp.parseArguments(System.err, new String[]{"--parameterizedTypeArgument", "parameterizedTypeArgumentValue",
                                                    "--parameterizedTypeListArgument", "parameterizedTypeListArgumentValue"});

        // Gather argument values of the raw type GatherArgumentValuesParameterizedTargetType, and make
        // sure that we match fully-parameterized declarations
        List<Pair<Field, GatherArgumentValuesParameterizedTargetType>> gatheredArguments =
                CommandLineParser.gatherArgumentValuesOfType(GatherArgumentValuesParameterizedTargetType.class, argumentSource);

        Assert.assertEquals(gatheredArguments.size(), 2, "Wrong number of arguments gathered");

        Assert.assertNotNull(gatheredArguments.get(0).getKey().getAnnotation(Argument.class), "Gathered argument is not annotated with an @Argument annotation");
        Assert.assertEquals(gatheredArguments.get(0).getKey().getAnnotation(Argument.class).fullName(), "parameterizedTypeArgument", "Wrong argument gathered");
        Assert.assertEquals(gatheredArguments.get(0).getValue().getValue(), "parameterizedTypeArgumentValue", "Wrong value for gathered argument");
        Assert.assertNotNull(gatheredArguments.get(1).getKey().getAnnotation(Argument.class), "Gathered argument is not annotated with an @Argument annotation");
        Assert.assertEquals(gatheredArguments.get(1).getKey().getAnnotation(Argument.class).fullName(), "parameterizedTypeListArgument", "Wrong argument gathered");
        Assert.assertEquals(gatheredArguments.get(1).getValue().getValue(), "parameterizedTypeListArgumentValue", "Wrong value for gathered argument");
    }

    /***************************************************************************************
     * End of tests and helper classes for CommandLineParser.gatherArgumentValuesOfType()
     ***************************************************************************************/
}
