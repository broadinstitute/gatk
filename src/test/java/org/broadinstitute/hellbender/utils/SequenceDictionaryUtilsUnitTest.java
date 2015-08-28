package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;


import static org.broadinstitute.hellbender.utils.SequenceDictionaryUtils.*;
import static org.broadinstitute.hellbender.utils.SequenceDictionaryUtils.SequenceDictionaryCompatibility.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SequenceDictionaryUtilsUnitTest extends BaseTest {

    private static Logger logger = Logger.getLogger(SequenceDictionaryUtilsUnitTest.class);


    @DataProvider( name = "SequenceDictionaryDataProvider" )
    public Object[][] generateSequenceDictionaryTestData() {
        final SAMSequenceRecord CHRM_HG19 = new SAMSequenceRecord("chrM", 16571);
        final SAMSequenceRecord CHR_NONSTANDARD1 = new SAMSequenceRecord("NonStandard1", 8675309);
        final SAMSequenceRecord CHR_NONSTANDARD2 = new SAMSequenceRecord("NonStandard2", 8675308);

        final Class<UserException.IncompatibleSequenceDictionaries> NO_COMMON_CONTIGS_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;
        final Class<UserException.IncompatibleSequenceDictionaries> UNEQUAL_COMMON_CONTIGS_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;
        final Class<UserException.LexicographicallySortedSequenceDictionary> NON_CANONICAL_HUMAN_ORDER_EXCEPTION = UserException.LexicographicallySortedSequenceDictionary.class;
        final Class<UserException.IncompatibleSequenceDictionaries> OUT_OF_ORDER_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;
        final Class<UserException.IncompatibleSequenceDictionaries> DIFFERENT_INDICES_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;

        final List<SimpleInterval> hg19AllContigsIntervalSet = Arrays.asList(
                new SimpleInterval("chrM", 1, 1),
                new SimpleInterval("chr1", 1, 1),
                new SimpleInterval("chr2", 1, 1),
                new SimpleInterval("chr10", 1, 1));
        final List<SimpleInterval> hg19PartialContigsIntervalSet = Arrays.asList(
                new SimpleInterval("chrM", 1, 1),
                new SimpleInterval("chr1", 1, 1));
        return new Object[][]  {
                // Identical dictionaries:

                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19),                        IDENTICAL, null, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), IDENTICAL, null, false},
                { Arrays.asList(CHR1_B37),                         Arrays.asList(CHR1_B37),                         IDENTICAL, null, false},
                { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),    Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),    IDENTICAL, null, false},

                // Dictionaries with a common subset:
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                                   COMMON_SUBSET, null, false},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19),                    Arrays.asList(CHRM_HG19, CHR1_HG19, CHR10_HG19),                              COMMON_SUBSET, UserException.IncompatibleSequenceDictionaries.class, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                        Arrays.asList(CHR1_HG19, CHR_NONSTANDARD2),                                   COMMON_SUBSET, null, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19),                        Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19),                                   COMMON_SUBSET, null, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19),                        Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHRM_HG19),                        COMMON_SUBSET, null, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD2),            COMMON_SUBSET, null, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHR2_HG19, CHR10_HG19),            COMMON_SUBSET, null, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHR2_HG19, CHR10_HG19, CHRM_HG19), COMMON_SUBSET, null, false},
                { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37, CHR_NONSTANDARD1),    Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37, CHR_NONSTANDARD2),               COMMON_SUBSET, null, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                              COMMON_SUBSET, UserException.IncompatibleSequenceDictionaries.class, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1),            COMMON_SUBSET, null, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                              COMMON_SUBSET, null, false},

                // Dictionaries with no common contigs:
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                     NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_B37),                      NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37), NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),             Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37), NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false},

                // Dictionaries with unequal common contigs:
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG18),                                          UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_B36),                                           Arrays.asList(CHR1_B37),                                           UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18),                   UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),                      Arrays.asList(CHR1_B36, CHR2_B36, CHR10_B36),                      UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18, CHR_NONSTANDARD2), UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG18, CHR2_HG18, CHR10_HG18), UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18),                   UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false},

                // One or both dictionaries in non-canonical human order:
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false},
                { Arrays.asList(CHR1_HG18, CHR10_HG18, CHR2_HG18), Arrays.asList(CHR1_HG18, CHR10_HG18, CHR2_HG18), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false},
                { Arrays.asList(CHR1_B37, CHR10_B37, CHR2_B37),    Arrays.asList(CHR1_B37, CHR10_B37, CHR2_B37),    NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false},
                { Arrays.asList(CHR1_B36, CHR10_B36, CHR2_B36),    Arrays.asList(CHR1_B36, CHR10_B36, CHR2_B36),    NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false},

                // Dictionaries with a common subset, but different relative ordering within that subset
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19, CHRM_HG19),                   OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHRM_HG19, CHR2_HG19, CHR1_HG19),                   OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHRM_HG19, CHR1_HG19),                   OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false},
                { Arrays.asList(CHR1_B37, CHR2_B37),              Arrays.asList(CHR2_B37, CHR1_B37),                                OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false},


                // Dictionaries with a common subset in the same relative order, but with different indices.
                // This will only throw an exception during validation if isReadsToReferenceComparison is true,
                // and there are intervals overlapping the misindexed contigs:

                // These have isReadsToReferenceComparison == true and overlapping intervals, so we expect an exception:
                { Arrays.asList(CHRM_HG19, CHR1_HG19),                                                 Arrays.asList(CHR1_HG19),                                          DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19),                    DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1),  DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19),  DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19 ),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19),                    DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true},

                // We expect exceptions for these because the common contigs' indicies don't match correctly.
                { Arrays.asList(CHR2_HG19, CHR10_HG19),                              Arrays.asList(CHR10_HG19),                       DIFFERENT_INDICES, UserException.IncompatibleSequenceDictionaries.class, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19),             Arrays.asList(CHR1_HG19, CHR2_HG19),             DIFFERENT_INDICES, UserException.IncompatibleSequenceDictionaries.class, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), DIFFERENT_INDICES, UserException.IncompatibleSequenceDictionaries.class, true},

                // These have isReadsToReferenceComparison == false, so we don't expect an exception:
                { Arrays.asList(CHRM_HG19, CHR1_HG19),                              Arrays.asList(CHR1_HG19),                       DIFFERENT_INDICES, null, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHRM_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19), DIFFERENT_INDICES, null, false},

                // tests for SUPERSET

                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                  Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                           Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1),    Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                                SUPERSET, null, false},
        };
    }

    @Test( dataProvider = "SequenceDictionaryDataProvider" )
    public void testSequenceDictionaryValidation( final List<SAMSequenceRecord> firstDictionaryContigs,
                                                  final List<SAMSequenceRecord> secondDictionaryContigs,
                                                  final SequenceDictionaryUtils.SequenceDictionaryCompatibility dictionaryCompatibility, //not needed by this test
                                                  final Class<? extends UserException> expectedExceptionUponValidation,
                                                  final boolean isReadsToReferenceComparison) {
        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);
        final String testDescription = String.format("First dictionary: %s  Second dictionary: %s",
                SequenceDictionaryUtils.getDictionaryAsString(firstDictionary),
                SequenceDictionaryUtils.getDictionaryAsString(secondDictionary));
        Exception exceptionThrown = null;
        try {
            SequenceDictionaryUtils.validateDictionaries(
                    "firstDictionary",
                    firstDictionary,
                    "secondDictionary",
                    secondDictionary,
                    isReadsToReferenceComparison);
        }
        catch ( Exception e ) {
            exceptionThrown = e;
        }
        if ( expectedExceptionUponValidation != null ) {
            Assert.assertTrue(exceptionThrown != null && expectedExceptionUponValidation.isInstance(exceptionThrown),
                    String.format("Expected exception %s but saw %s instead. %s",
                            expectedExceptionUponValidation.getSimpleName(),
                            exceptionThrown == null ? "no exception" : exceptionThrown.getClass().getSimpleName(),
                            testDescription));
        }
        else {
            Assert.assertTrue(exceptionThrown == null,
                    String.format("Expected no exception but saw exception %s instead. %s",
                            exceptionThrown != null ? exceptionThrown.getClass().getSimpleName() : "none",
                            testDescription));
        }
    }

    @Test( dataProvider = "SequenceDictionaryDataProvider" )
    public void testSequenceDictionaryComparison( final List<SAMSequenceRecord> firstDictionaryContigs,
                                                  final List<SAMSequenceRecord> secondDictionaryContigs,
                                                  final SequenceDictionaryUtils.SequenceDictionaryCompatibility dictionaryCompatibility,
                                                  final Class<? extends UserException> expectedExceptionUponValidation,
                                                  final boolean isReadsToReferenceComparison) {

        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);
        final String testDescription = String.format("First dictionary: %s  Second dictionary: %s",
                SequenceDictionaryUtils.getDictionaryAsString(firstDictionary),
                SequenceDictionaryUtils.getDictionaryAsString(secondDictionary));

        final SequenceDictionaryUtils.SequenceDictionaryCompatibility reportedCompatibility =
                SequenceDictionaryUtils.compareDictionaries(firstDictionary, secondDictionary);

        Assert.assertTrue(reportedCompatibility == dictionaryCompatibility,
                String.format("Dictionary comparison should have returned %s but instead returned %s. %s",
                        dictionaryCompatibility, reportedCompatibility, testDescription));
    }

    private SAMSequenceDictionary createSequenceDictionary( final List<SAMSequenceRecord> contigs ) {
        final List<SAMSequenceRecord> clonedContigs = new ArrayList<SAMSequenceRecord>(contigs.size());

        // Clone the individual SAMSequenceRecords to avoid contig-index issues with shared objects
        // across multiple dictionaries in tests
        for ( SAMSequenceRecord contig : contigs ) {
            clonedContigs.add(contig.clone());
        }

        return new SAMSequenceDictionary(clonedContigs);
    }
}