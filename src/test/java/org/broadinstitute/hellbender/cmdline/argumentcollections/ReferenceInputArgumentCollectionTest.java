package org.broadinstitute.hellbender.cmdline.argumentcollections;


import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
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
        CommandLineParser clp = new CommandLineParser(optional);
        clp.parseArguments(System.out, args);
    }

    @Test(expectedExceptions = UserException.CommandLineException.class)
    public void testRequiredIsRequired(){
        String[] args = {};
        WithRequiredReferenceCollection required = new WithRequiredReferenceCollection();
        CommandLineParser clp = new CommandLineParser(required);
        clp.parseArguments(System.out, args);
    }
}