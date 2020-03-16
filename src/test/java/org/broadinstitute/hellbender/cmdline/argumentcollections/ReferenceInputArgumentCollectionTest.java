package org.broadinstitute.hellbender.cmdline.argumentcollections;


import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.testng.Assert;
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


    @Test
    public void testGetNullPath(){
        final ReferenceInputArgumentCollection nullPath = new ReferenceInputArgumentCollection(){
            private static final long serialVersionUID = 0L;
            @Override public GATKPathSpecifier getReferenceInputPath() { return null;}
        };

        Assert.assertNull(nullPath.getReferencePath());
    }
}