package org.broadinstitute.hellbender.cmdline.argumentcollections;


import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

public final class ReferenceInputArgumentCollectionTest {
    private static class WithOptionalReferenceCollection {
        @ArgumentCollection
        ReferenceInputArgumentCollection ric = new OptionalReferenceInputArgumentCollection();
    }

    private static class WithRequiredReferenceCollection {
        @ArgumentCollection
        ReferenceInputArgumentCollection ric = new RequiredReferenceInputArgumentCollection();
    }

    @Test
    public void testOptionalIsOptional(){
        String[] args = {};
        WithOptionalReferenceCollection optional = new WithOptionalReferenceCollection();
        CommandLineParser clp = new CommandLineArgumentParser(optional);
        clp.parseArguments(System.out, args);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testRequiredIsRequired(){
        String[] args = {};
        WithRequiredReferenceCollection required = new WithRequiredReferenceCollection();
        CommandLineParser clp = new CommandLineArgumentParser(required);
        clp.parseArguments(System.out, args);
    }
}