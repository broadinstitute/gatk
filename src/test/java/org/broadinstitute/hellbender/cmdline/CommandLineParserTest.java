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
package org.broadinstitute.hellbender.cmdline;

import htsjdk.samtools.util.CollectionUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

public class CommandLineParserTest {
    enum FrobnicationFlavor {
        FOO, BAR, BAZ
    }

    @CommandLineProgramProperties(
            usage = "Usage: frobnicate [arguments] input-file output-file\n\nRead input-file, frobnicate it, and write frobnicated results to output-file\n",
            usageShort = "Read input-file, frobnicate it, and write frobnicated results to output-file"
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

        @Argument(doc="Allowed shmiggle types.", minElements=1, maxElements = 3)
        public List<String> SHMIGGLE_TYPE = new ArrayList<>();

        @Argument
        public Boolean TRUTHINESS = false;
    }

    @CommandLineProgramProperties(
            usage = "Usage: framistat [arguments]\n\nCompute the plebnick of the freebozzle.\n",
            usageShort = "ompute the plebnick of the freebozzle"
    )
    class ArgumentsWithoutPositional {
        public static final int DEFAULT_FROBNICATION_THRESHOLD = 20;
        @Argument(shortName="T", doc="Frobnication threshold setting.")
        public Integer FROBNICATION_THRESHOLD = DEFAULT_FROBNICATION_THRESHOLD;

        @Argument
        public FrobnicationFlavor FROBNICATION_FLAVOR;

        @Argument(doc="Allowed shmiggle types.", minElements=1, maxElements = 3)
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


    @Test
    public void testUsage() {
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, false);
    }

    @Test
    public void testUsageWithDefault() {
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, true);
    }

    @Test
    public void testUsageWithoutPositional() {
        final ArgumentsWithoutPositional fo = new ArgumentsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, false);
    }

    @Test
    public void testUsageWithoutPositionalWithDefault() {
        final ArgumentsWithoutPositional fo = new ArgumentsWithoutPositional();
        final CommandLineParser clp = new CommandLineParser(fo);
        clp.usage(System.out, true);
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
                        "--SHMIGGLE_TYPE [shmiggle1, shmiggle2] --TRUTHINESS true    --help false --version false");
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

    @Test
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
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
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
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
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
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
    public void testNotEnoughOfListArgument() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
    public void testTooManyListArgument() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
                "--SHMIGGLE_TYPE","shmiggle3",
                "--SHMIGGLE_TYPE","shmiggle4",
                "positional1",
                "positional2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
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
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
    public void testNotEnoughPositional() {
        final String[] args = {
                "--FROBNICATION_FLAVOR","BAR",
                "--TRUTHINESS","False",
                "--SHMIGGLE_TYPE","shmiggle1",
                "--SHMIGGLE_TYPE","shmiggle2",
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
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
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
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
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }

    @Test
    public void testArgumentsFile() throws Exception {
        final File argumentsFile = File.createTempFile("clp.", ".arguments");
        argumentsFile.deleteOnExit();
        final PrintWriter writer = new PrintWriter(argumentsFile);
        writer.println("-T 18");
        writer.println("--TRUTHINESS");
        writer.println("--SHMIGGLE_TYPE shmiggle0");
        writer.println("--"+SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME + " " +argumentsFile.getPath());
        //writer.println("--STRANGE_ARGUMENT shmiggle0");
        writer.close();
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
    @Test
    public void testArgumentsFileWithDisallowedOverride() throws Exception {
        final File argumentsFile = File.createTempFile("clp.", ".arguments");
        argumentsFile.deleteOnExit();
        final PrintWriter writer = new PrintWriter(argumentsFile);
        writer.println("--T 18");
        writer.close();
        final String[] args = {
                "--T","17",
                "--"+SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME ,argumentsFile.getPath()
        };
        final FrobnicateArguments fo = new FrobnicateArguments();
        final CommandLineParser clp = new CommandLineParser(fo);
        Assert.assertFalse(clp.parseArguments(System.err, args));
    }
    
    @DataProvider(name="mutexScenarios")
    public Object[][] mutexScenarios() {
        return new Object[][] {
                { "pass", new String[] {"-A","1", "-B","2"}, true },
                { "no args", new String[0], false },
                { "1 of group required", new String[] {"-A","1"}, false },
                { "mutex", new String[]  {"-A","1", "-Y","3"}, false },
                { "mega mutex", new String[]  {"-A","1", "-B","2", "-Y","3", "-Z","1", "-M","2", "-N","3"}, false }
        };
    }
    
    @Test(dataProvider="mutexScenarios")
    public void testMutex(final String testName, final String[] args, final boolean expected) {
        final MutexArguments o = new MutexArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertEquals(clp.parseArguments(System.err, args), expected);
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
        Assert.fail("Exception should have been thrown");
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
    public void testAddToDefaultValuesListArgument() {
        final CollectionWithDefaultValuesArguments o = new CollectionWithDefaultValuesArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"--LIST","baz", "--LIST","frob"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.LIST, CollectionUtil.makeList("foo", "bar", "baz", "frob"));
    }

    @Test
    public void testCommandLineToolTest() {
        final CollectionWithDefaultValuesArguments o = new CollectionWithDefaultValuesArguments();
        final CommandLineParser clp = new CommandLineParser(o);
        final String[] args = {"--LIST","baz", "--LIST","frob"};
        Assert.assertTrue(clp.parseArguments(System.err, args));
        Assert.assertEquals(o.LIST, CollectionUtil.makeList("foo", "bar", "baz", "frob"));
    }

    class ArgumentsWithArgumentCollectionAgain {
        @ArgumentCollection(doc="Doc for inner FROB")
        public ArgumentsWithoutPositional FROB = new ArgumentsWithoutPositional();
    }


    class FlagCollection{
        @Argument
        public boolean flag1 = false;

        @Argument
        public boolean flag2 = true;

        @ArgumentCollection
        SpecialArgumentsCollection special = new SpecialArgumentsCollection();

    }

    @Test
       public void testFlagNoArgument(){
        final FlagCollection o = new FlagCollection();
        final CommandLineParser clp = new CommandLineParser(o);
        Assert.assertTrue(clp.parseArguments(System.err, new String[]{"--flag1"}));
        Assert.assertTrue(o.flag1);
    }

    @Test
    public void testFlagsWithArguments(){
        final FlagCollection o = new FlagCollection();
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
        Assert.fail(); //shouldn't reach here
    }

    class PrivateArgument{
        @Argument
        private boolean privateArgument = false;

        @Argument
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
        Assert.assertEquals(o.positionals, Arrays.asList(1,2));
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
}
