package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

public final class FeatureInputArgumentCollectionTest extends GATKBaseTest {

    @Test(expectedExceptions = CommandLineException.class)
    public void testRequiredIsRequired(){
        Object req = new Object(){
            @ArgumentCollection
            private RequiredFeatureInputArgumentCollection ric = new RequiredFeatureInputArgumentCollection();
        };
        CommandLineParser clp = new CommandLineArgumentParser(req);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test
    public void testOptionalIsOptional(){
        Object req = new Object(){
            @ArgumentCollection
            private OptionalFeatureInputArgumentCollection ric = new OptionalFeatureInputArgumentCollection();
        };
        CommandLineParser clp = new CommandLineArgumentParser(req);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }
}