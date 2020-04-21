package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

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

    // todo this should be tested here instead
    //@Test(dataProvider = "manuallySpecifiedIndexTestData", expectedExceptions = UserException.class)
  /*  public void testManuallySpecifiedIndicesWrongNumberOfIndices(final List<Path> bams, final List<Path> indices ) {
        final List<Path> wrongIndices = new ArrayList<>();
        wrongIndices.add(indices.get(0)); // Add one index, but not the other

        final ReadsDataSource readsSource = new ReadsDataSource(bams, wrongIndices);
    }*/
}