package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SequenceDictionaryValidationArgumentCollectionTest {

    private static class StandardArgumentCollection {
        @ArgumentCollection
        public SequenceDictionaryValidationArgumentCollection standard = new SequenceDictionaryValidationArgumentCollection.StandardValidationCollection();
    }

    @Test
    public void testStandardArgumentCollectionDefaultsToTrue(){
        Assert.assertTrue(new SequenceDictionaryValidationArgumentCollection.StandardValidationCollection().performSequenceDictionaryValidation());
    }

    @Test
    public void testStandardArgumentCollectionCanBeDisabled(){
        final String[] disabled = {"--"+StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME};
        StandardArgumentCollection std = new StandardArgumentCollection();
        CommandLineParser clp = new CommandLineArgumentParser(std);
        clp.parseArguments(System.out, disabled);
        Assert.assertFalse(std.standard.performSequenceDictionaryValidation());
    }

    @Test
    public void testNoValidationDefaultsToFalse(){
        Assert.assertFalse(new SequenceDictionaryValidationArgumentCollection.NoValidationCollection().performSequenceDictionaryValidation());
    }

}