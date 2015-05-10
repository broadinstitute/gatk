package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public final class ReadInputArgumentCollectionTest {

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testRequiredIsRequired(){
        Object req = new Object(){
            @ArgumentCollection
            private ReadInputArgumentCollection ric = new RequiredReadInputArgumentCollection();
        };
        CommandLineParser clp = new CommandLineParser(req);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }

    @Test
    public void testOptionalIsOptional(){
        Object req = new Object(){
            @ArgumentCollection
            private ReadInputArgumentCollection ric = new OptionalReadInputArgumentCollection();
        };
        CommandLineParser clp = new CommandLineParser(req);
        String[] args = {};
        clp.parseArguments(System.out, args);
    }
}