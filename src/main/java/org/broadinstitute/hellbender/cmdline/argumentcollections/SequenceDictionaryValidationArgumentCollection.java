package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

/**
 * interface for argument collections that control how sequence dictionary validation should be handled
 */
public interface SequenceDictionaryValidationArgumentCollection {

    /**
     * Should sequence dictionary validation be performed
     * @return true if the tool should perform sequence dictionary validation
     */
    boolean performSequenceDictionaryValidation();


    /**
     * most tools will want to use this, it defaults to performing sequence dictionary validation but provides the option
     * to disable it
     */
    class StandardValidationCollection implements SequenceDictionaryValidationArgumentCollection {
        @Argument(fullName = StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, shortName = StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, doc = "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!", optional = true)
        private boolean disableSequenceDictionaryValidation = false;

        @Override
        public boolean performSequenceDictionaryValidation() {
            return !disableSequenceDictionaryValidation;
        }
    }

    /**
     * doesn't provide a configuration argument, and always returns false, useful for tools that do not want to perform
     * sequence dictionary validation, like aligners
     */
    class NoValidationCollection implements SequenceDictionaryValidationArgumentCollection {
        @Override
        public boolean performSequenceDictionaryValidation() {
            return false;
        }
    }
}
