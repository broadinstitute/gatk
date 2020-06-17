package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class ReadInputArgumentCollectionTest {

    @Test(expectedExceptions = CommandLineException.class)
    public void testRequiredIsRequired(){
        Object req = new Object(){
            @ArgumentCollection
            private ReadInputArgumentCollection ric = new RequiredReadInputArgumentCollection();
        };
        CommandLineParser clp = new CommandLineArgumentParser(req);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test
    public void testOptionalIsOptional(){
        Object req = new Object(){
            @ArgumentCollection
            private ReadInputArgumentCollection ric = new OptionalReadInputArgumentCollection();
        };
        CommandLineParser clp = new CommandLineArgumentParser(req);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test
    public void testManuallySpecifiedIndicesWrongNumberOfIndices() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput("file.bam");
        args.addInput("otherFile.bam");
        args.add(StandardArgumentDefinitions.READ_INDEX_LONG_NAME, "file1.bai");
        ReadInputArgumentCollection ric = getInitializedReadArgumentInputCollection(args);
        assertAllGettersThrow(ric, UserException.MissingIndex.class);
    }

    private static void assertAllGettersThrow(final ReadInputArgumentCollection ric, final Class<? extends Exception> expectedException) {
        Assert.assertThrows(expectedException, ric::getReadIndexPairs);
        Assert.assertThrows(expectedException, ric::getReadPaths);
        Assert.assertThrows(expectedException, ric::getReadPaths);
    }

    @Test
    public void testJsonAndIndexesCantBeSpecifiedTogether() {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput("file.bam");
        args.add(StandardArgumentDefinitions.READ_INDEX_LONG_NAME, "file1.bai");
        args.addInput("bundle.json");
        args.add(StandardArgumentDefinitions.READ_INDEX_LONG_NAME, "file2.bai"); //the check for the correct number of indexes/inputs takes effect
        ReadInputArgumentCollection ric = getInitializedReadArgumentInputCollection(args);
        assertAllGettersThrow(ric, UserException.class);
    }

    private static ReadInputArgumentCollection getInitializedReadArgumentInputCollection(final ArgumentsBuilder args) {
        ReadInputArgumentCollection ric = new OptionalReadInputArgumentCollection();
        CommandLineParser clp = new CommandLineArgumentParser(ric);
        clp.parseArguments(System.out, args.getArgsArray());
        return ric;
    }

    @Test
    public void testReadBundleLoading(){

    }
}