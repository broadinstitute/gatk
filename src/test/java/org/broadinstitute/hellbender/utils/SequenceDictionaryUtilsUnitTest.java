package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.utils.SequenceDictionaryUtils.*;
import static org.broadinstitute.hellbender.utils.SequenceDictionaryUtils.SequenceDictionaryCompatibility.*;

public final class SequenceDictionaryUtilsUnitTest extends GATKBaseTest {

    private static Logger logger = LogManager.getLogger(SequenceDictionaryUtilsUnitTest.class);

    @DataProvider( name = "testSequenceRecordsAreEquivalentDataProvider" )
    public Object[][] testSequenceRecordsAreEquivalentDataProvider() {
        final SAMSequenceRecord CHRM_HG19 = new SAMSequenceRecord("chrM", 16571);
        final SAMSequenceRecord CHR_NONSTANDARD1 = new SAMSequenceRecord("NonStandard1", 8675309);
        final SAMSequenceRecord CHR1_HG19_WITH_UNKNOWN_LENGTH = new SAMSequenceRecord(CHR1_HG19.getSequenceName(), SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH);
        final SAMSequenceRecord CHR1_HG19_WITH_DIFFERENT_LENGTH = new SAMSequenceRecord(CHR1_HG19.getSequenceName(), 123456);
        return new Object[][]{
                {CHR1_HG19, CHR1_HG19, true},
                {CHR1_HG19, CHRM_HG19, false},
                {CHR1_HG19, CHR_NONSTANDARD1, false},
                {null, null, true},
                {CHR1_HG19, null, false},
                {null, CHR1_HG19, false},
                {CHR1_HG19, CHR1_HG19_WITH_UNKNOWN_LENGTH, true},
                {CHR1_HG19, CHR1_HG19_WITH_DIFFERENT_LENGTH, false},
                {CHR1_HG19_WITH_UNKNOWN_LENGTH, CHR1_HG19, true},
                {CHR1_HG19_WITH_DIFFERENT_LENGTH, CHR1_HG19, false},
        };
    }

    @Test(dataProvider = "testSequenceRecordsAreEquivalentDataProvider")
    public void testSequenceRecordsAreEquivalent(final SAMSequenceRecord one, final SAMSequenceRecord two, final boolean expected){
        final boolean actual = SequenceDictionaryUtils.sequenceRecordsAreEquivalent(one, two);
        Assert.assertEquals(actual, expected);
    }

    @DataProvider( name = "SequenceDictionaryDataProvider" )
    public Object[][] generateSequenceDictionaryTestData() {
        final SAMSequenceRecord CHRM_HG19 = new SAMSequenceRecord("chrM", 16571);
        final SAMSequenceRecord CHR_NONSTANDARD1 = new SAMSequenceRecord("NonStandard1", 8675309);
        final SAMSequenceRecord CHR_NONSTANDARD2 = new SAMSequenceRecord("NonStandard2", 8675308);
        final SAMSequenceRecord CHR1_HG19_WITH_UNKNOWN_LENGTH = new SAMSequenceRecord(CHR1_HG19.getSequenceName(), SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH);
        final SAMSequenceRecord CHR1_HG19_WITH_DIFFERENT_LENGTH = new SAMSequenceRecord(CHR1_HG19.getSequenceName(), 123456);

        final SAMSequenceRecord CHR1_HG19_WITH_ATTRIBUTES = new SAMSequenceRecord(CHR1_HG19.getSequenceName(), CHR1_HG19.getSequenceLength());
        CHR1_HG19_WITH_ATTRIBUTES.setAttribute("M5", "0dec9660ec1efaaf33281c0d5ea2560f");
        CHR1_HG19_WITH_ATTRIBUTES.setAttribute("UR", "file:/foo/bar");

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
                {Arrays.asList(CHR1_HG19),                         Arrays.asList(CHR1_HG19),                        IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19),                        IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19),                        IDENTICAL, null, true,  false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19),                        IDENTICAL, null, false, true},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19),                        IDENTICAL, null, true,  true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_B37),                         Arrays.asList(CHR1_B37),                         IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),    Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),    IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19_WITH_UNKNOWN_LENGTH),    IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_HG19_WITH_UNKNOWN_LENGTH),    Arrays.asList(CHR1_HG19),                        IDENTICAL, null, false, false},

                // Dictionaries with a common subset:
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                                   COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                                   COMMON_SUBSET, null, false, true},
                // If requireSuperset == true, we should get an exception upon COMMON_SUBSET:
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19),                    Arrays.asList(CHRM_HG19, CHR1_HG19, CHR10_HG19),                              COMMON_SUBSET, UserException.IncompatibleSequenceDictionaries.class, true, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                        Arrays.asList(CHR1_HG19, CHR_NONSTANDARD2),                                   COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19),                        Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19),                                   COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19),                        Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHRM_HG19),                        COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD2),            COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHR2_HG19, CHR10_HG19),            COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHR2_HG19, CHR10_HG19, CHRM_HG19), COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37, CHR_NONSTANDARD1),    Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37, CHR_NONSTANDARD2),               COMMON_SUBSET, null, false, false},
                // If requireSuperset == true, we should get an exception upon COMMON_SUBSET:
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                              COMMON_SUBSET, UserException.IncompatibleSequenceDictionaries.class, true, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1),            COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                              COMMON_SUBSET, null, false, false},
                // If checkContigOrdering == false, ordering of the common contigs should not matter:
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR2_HG19, CHR1_HG19, CHR10_HG19),                              COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR10_HG19, CHR2_HG19, CHR1_HG19),                              COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR2_HG19, CHR10_HG19, CHR1_HG19),                              COMMON_SUBSET, null, false, false},

                // Dictionaries with no common contigs:
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                     NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                     NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                     NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, true,  false},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                     NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, true,  true},
                { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_B37),                      NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37), NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),             Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37), NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, false},

                // Dictionaries with unequal common contigs:
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG19_WITH_DIFFERENT_LENGTH),                    UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19_WITH_DIFFERENT_LENGTH),                    Arrays.asList(CHR1_HG19),                                          UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG18),                                          UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG18),                                          UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, true,  false},
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG18),                                          UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG18),                                          UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, true,  true},
                { Arrays.asList(CHR1_B36),                                           Arrays.asList(CHR1_B37),                                           UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18),                   UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),                      Arrays.asList(CHR1_B36, CHR2_B36, CHR10_B36),                      UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18, CHR_NONSTANDARD2), UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG18, CHR2_HG18, CHR10_HG18), UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18),                   UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, false},

                // One or both dictionaries in non-canonical human order:
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, true,  true},
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, true,  true},
                { Arrays.asList(CHR1_HG18, CHR10_HG18, CHR2_HG18), Arrays.asList(CHR1_HG18, CHR10_HG18, CHR2_HG18), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_B37, CHR10_B37, CHR2_B37),    Arrays.asList(CHR1_B37, CHR10_B37, CHR2_B37),    NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_B36, CHR10_B36, CHR2_B36),    Arrays.asList(CHR1_B36, CHR10_B36, CHR2_B36),    NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, true},
                // If checkContigOrdering == false, we should not get NON_CANONICAL_HUMAN_ORDER:
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), IDENTICAL, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19, CHR_NONSTANDARD1), COMMON_SUBSET, null, false, false},

                // Dictionaries with a common subset, but different relative ordering within that subset
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, true,  true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, true,  true},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19, CHRM_HG19),                   OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHRM_HG19, CHR2_HG19, CHR1_HG19),                   OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHRM_HG19, CHR1_HG19),                   OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, true},
                { Arrays.asList(CHR1_B37, CHR2_B37),              Arrays.asList(CHR2_B37, CHR1_B37),                                OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, true},
                // If checkContigOrdering == false, we should not get OUT_OF_ORDER:
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              SUPERSET,      null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19, CHR_NONSTANDARD1),            COMMON_SUBSET, null, false, false},

                // Dictionaries with a common subset in the same relative order, but with different indices.
                // This will only throw an exception during validation if checkContigOrdering is true

                // These have checkContigOrdering == true, so we expect DIFFERENT_INDICES and an exception:
                { Arrays.asList(CHRM_HG19, CHR1_HG19),                                                 Arrays.asList(CHR1_HG19),                                          DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                // Setting requireSuperset == true should make no difference here (we should still get DIFFERENT_INDICES and an exception):
                { Arrays.asList(CHRM_HG19, CHR1_HG19),                                                 Arrays.asList(CHR1_HG19),                                          DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true,  true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19),                    DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1),  DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19),  DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19),                               DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19 ),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19),                    DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, false, true},

                // Same test cases as above, but these have checkContigOrdering == false, so we expect SUPERSET or COMMON_SUBSET instead of DIFFERENT_INDICES, and no exception:
                { Arrays.asList(CHRM_HG19, CHR1_HG19),                                                 Arrays.asList(CHR1_HG19),                                          SUPERSET,      null, false, false},
                { Arrays.asList(CHRM_HG19, CHR1_HG19),                                                 Arrays.asList(CHR1_HG19),                                          SUPERSET,      null, true,  false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19),                    COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1),  COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19),  COMMON_SUBSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19),                               SUPERSET,      null, false, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   SUPERSET,      null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19 ),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19),                    SUPERSET,      null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), SUPERSET,      null, false, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), SUPERSET,      null, false, false},

                // tests for SUPERSET
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                  Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                  Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, false, true},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                  Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, true,  false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                  Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, true,  true},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                           Arrays.asList(CHR1_HG19),                                                       SUPERSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1),    Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                                SUPERSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHR10_HG19),    Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                                SUPERSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHR10_HG19),    Arrays.asList(CHR1_HG19, CHR2_HG19),                                            SUPERSET, null, false, false},
                // Extended attributes should be ignored when determining whether a superset exists:
                { Arrays.asList(CHR1_HG19, CHR2_HG19),                                  Arrays.asList(CHR1_HG19_WITH_ATTRIBUTES),                                       SUPERSET, null, false, false},
                { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHRM_HG19),           Arrays.asList(CHR1_HG19_WITH_ATTRIBUTES, CHR10_HG19),                           SUPERSET, null, false, false}
        };
    }

    @Test( dataProvider = "SequenceDictionaryDataProvider" )
    public void testSequenceDictionaryValidation( final List<SAMSequenceRecord> firstDictionaryContigs,
                                                  final List<SAMSequenceRecord> secondDictionaryContigs,
                                                  final SequenceDictionaryUtils.SequenceDictionaryCompatibility dictionaryCompatibility, //not needed by this test
                                                  final Class<? extends UserException> expectedExceptionUponValidation,
                                                  final boolean requireSuperset,
                                                  final boolean checkContigOrdering) {
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
                    requireSuperset,
                    checkContigOrdering);
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
                                                  final boolean requireSuperset,
                                                  final boolean checkContigOrdering) {

        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);
        final String testDescription = String.format("First dictionary: %s  Second dictionary: %s",
                SequenceDictionaryUtils.getDictionaryAsString(firstDictionary),
                SequenceDictionaryUtils.getDictionaryAsString(secondDictionary));

        final SequenceDictionaryUtils.SequenceDictionaryCompatibility reportedCompatibility =
                SequenceDictionaryUtils.compareDictionaries(firstDictionary, secondDictionary, checkContigOrdering);

        Assert.assertTrue(reportedCompatibility == dictionaryCompatibility,
                String.format("Dictionary comparison should have returned %s but instead returned %s. %s",
                        dictionaryCompatibility, reportedCompatibility, testDescription));
    }

    @DataProvider(name = "StandardValidationIgnoresContigOrderData")
    public Object[][] getStandardValidationIgnoresContigOrderData() {
        return new Object[][] {
                { Arrays.asList(CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19) },
                { Arrays.asList(CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19, CHR10_HG19) },
                { Arrays.asList(CHR1_HG19, CHR2_HG19), Arrays.asList(CHR10_HG19, CHR2_HG19, CHR1_HG19) },
                { Arrays.asList(CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR10_HG19, CHR1_HG19) },

        };
    }

    @Test(dataProvider = "StandardValidationIgnoresContigOrderData")
    public void testStandardValidationIgnoresContigOrder( final List<SAMSequenceRecord> firstDictionaryContigs, final List<SAMSequenceRecord> secondDictionaryContigs ) {
        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);

        // Standard validation (the overload of validateDictionaries() that doesn't take any boolean args)
        // should ignore differences in ordering of common contigs, so we shouldn't get an exception here
        SequenceDictionaryUtils.validateDictionaries("first", firstDictionary, "second", secondDictionary);
    }

    @DataProvider(name = "NonSupersetData")
    public Object[][] getNonSupersetData() {
        return new Object[][] {
                { Arrays.asList(CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19, CHR10_HG19) },
                { Arrays.asList(CHR1_HG19),            Arrays.asList(CHR10_HG19, CHR2_HG19, CHR1_HG19) }
        };
    }

    @Test(dataProvider = "NonSupersetData")
    public void testStandardValidationDoesNotRequireSuperset( final List<SAMSequenceRecord> firstDictionaryContigs, final List<SAMSequenceRecord> secondDictionaryContigs ) {
        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);

        // Standard validation (the overload of validateDictionaries() that doesn't take any boolean args)
        // should not require a superset relationship, so we shouldn't get an exception here
        SequenceDictionaryUtils.validateDictionaries("first", firstDictionary, "second", secondDictionary);
    }

    @Test(dataProvider = "NonSupersetData", expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testCRAMValidationDoesRequireSuperset( final List<SAMSequenceRecord> refDictionaryContigs, final List<SAMSequenceRecord> cramDictionaryContigs ) {
        final SAMSequenceDictionary refDictionary = createSequenceDictionary(refDictionaryContigs);
        final SAMSequenceDictionary cramDictionary = createSequenceDictionary(cramDictionaryContigs);

        // CRAM validation against the reference SHOULD require a superset relationship, so we should
        // get an exception here
        SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(refDictionary, cramDictionary);
    }

    @DataProvider(name = "SupersetData")
    public Object[][] getSupersetData() {
        return new Object[][] {
                { Arrays.asList(CHR2_HG19, CHR1_HG19, CHR10_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19, CHR10_HG19)}, //exactly same
                { Arrays.asList(CHR2_HG19, CHR1_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19) },
                { Arrays.asList(CHR10_HG19, CHR2_HG19, CHR1_HG19), Arrays.asList(CHR1_HG19) }
        };
    }

    @Test(dataProvider = "SupersetData")
    public void testCRAMValidationDoesAcceptSuperset( final List<SAMSequenceRecord> refDictionaryContigs, final List<SAMSequenceRecord> cramDictionaryContigs ) {
        final SAMSequenceDictionary refDictionary = createSequenceDictionary(refDictionaryContigs);
        final SAMSequenceDictionary cramDictionary = createSequenceDictionary(cramDictionaryContigs);

        //In these inputs , cram contigs are subsets of ref contigs and so it should be accepted
        SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(refDictionary, cramDictionary);
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

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetContigNamesListExpectingException() {
        getContigNamesList(null);
    }

    @Test
    public void testGetContigNamesList() {

        final SAMSequenceDictionary samSequenceDictionary = new SAMSequenceDictionary(Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37));

        Assert.assertEquals(getContigNamesList(samSequenceDictionary), Arrays.asList("1", "2", "10"));
    }
}