package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryUtils;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.*;
import java.util.stream.Collectors;

/**
 *
 * A series of utility functions that enable the GATK to compare two sequence dictionaries -- from the reference,
 * from BAMs, or from feature sources -- for consistency.  The system supports two basic modes: get an enum state that
 * describes at a high level the consistency between two dictionaries, or a validateDictionaries that will
 * blow up with a UserException if the dicts are too incompatible.
 *
 * Dictionaries are tested for contig name overlaps, consistency in ordering in these overlap set, and length,
 * if available.
 */
public final class SequenceDictionaryUtils {

    private SequenceDictionaryUtils(){}

    /**
     * Tests for compatibility between two sequence dictionaries, using standard validation settings appropriate
     * for the GATK. If the dictionaries are incompatible, then UserExceptions are thrown with detailed error messages.
     *
     * The standard validation settings used by this method are:
     *
     * -Require the dictionaries to share a common subset of equivalent contigs
     *
     * -Do not require dict1 to be a superset of dict2.
     *
     * -Do not perform checks related to contig ordering: don't throw if the common contigs are in
     *  different orders with respect to each other, occur at different absolute indices, or are
     *  lexicographically sorted human dictionaries. GATK uses contig names rather than contig
     *  indices, and so should not be sensitive to contig ordering issues.
     *
     * For comparing a CRAM dictionary against a reference dictionary, call
     * {@link #validateCRAMDictionaryAgainstReference(SAMSequenceDictionary, SAMSequenceDictionary)} instead.
     *
     * @param name1 name associated with dict1
     * @param dict1 the sequence dictionary dict1
     * @param name2 name associated with dict2
     * @param dict2 the sequence dictionary dict2
     */
    public static void validateDictionaries( final String name1,
                                             final SAMSequenceDictionary dict1,
                                             final String name2,
                                             final SAMSequenceDictionary dict2) {
        final boolean requireSuperset = false;
        final boolean checkContigOrdering = false;

        validateDictionaries(name1, dict1, name2, dict2, requireSuperset, checkContigOrdering);
    }

    /**
     * Tests for compatibility between a reference dictionary and a CRAM dictionary, using appropriate
     * validation settings. If the dictionaries are incompatible, then UserExceptions are thrown with
     * detailed error messages.
     *
     * The standard validation settings used by this method are:
     *
     * -Require the reference dictionary to be a superset of the cram dictionary
     *
     * -Do not perform checks related to contig ordering: don't throw if the common contigs are in
     *  different orders with respect to each other, occur at different absolute indices, or are
     *  lexicographically sorted human dictionaries. GATK uses contig names rather than contig
     *  indices, and so should not be sensitive to contig ordering issues.
     *
     * @param referenceDictionary the sequence dictionary for the reference
     * @param cramDictionary sequence dictionary from a CRAM file
     */
    public static void validateCRAMDictionaryAgainstReference( final SAMSequenceDictionary referenceDictionary,
                                                               final SAMSequenceDictionary cramDictionary ) {
        // For CRAM, we require the reference dictionary to be a superset of the reads dictionary
        final boolean requireSuperset = true;
        final boolean checkContigOrdering = false;

        validateDictionaries("reference", referenceDictionary, "reads", cramDictionary, requireSuperset, checkContigOrdering);
    }


    /**
     * Tests for compatibility between two sequence dictionaries.  If the dictionaries are incompatible, then
     * UserExceptions are thrown with detailed error messages.
     *
     * Two sequence dictionaries are compatible if they share a common subset of equivalent contigs,
     * where equivalent contigs are defined as having the same name and length.
     *
     * @param name1 name associated with dict1
     * @param dict1 the sequence dictionary dict1
     * @param name2 name associated with dict2
     * @param dict2 the sequence dictionary dict2
     * @param requireSuperset if true, require that dict1 be a superset of dict2, rather than dict1 and dict2 sharing a common subset
     * @param checkContigOrdering if true, require common contigs to be in the same relative order with respect to each other
     *                            and occur at the same absolute indices, and forbid lexicographically-sorted human dictionaries
     */
    public static void validateDictionaries( final String name1,
                                             final SAMSequenceDictionary dict1,
                                             final String name2,
                                             final SAMSequenceDictionary dict2,
                                             final boolean requireSuperset,
                                             final boolean checkContigOrdering ) {
        Utils.nonNull(dict1, "Something went wrong with sequence dictionary detection, check that "+name1+" has a valid sequence dictionary");
        Utils.nonNull(dict2, "Something went wrong with sequence dictionary detection, check that "+name2+" has a valid sequence dictionary");

        final SAMSequenceDictionaryUtils.SequenceDictionaryCompatibility type = SAMSequenceDictionaryUtils.compareDictionaries(dict1, dict2, checkContigOrdering);

        switch ( type ) {
            case IDENTICAL:
                return;
            case SUPERSET:
                return;
            case COMMON_SUBSET:
                if ( requireSuperset ) {
                    final Set<String> contigs1 = dict1.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toSet());
                    final List<String> missingContigs = dict2.getSequences().stream()
                            .map(SAMSequenceRecord::getSequenceName)
                            .filter(contig -> !contigs1.contains(contig))
                            .collect(Collectors.toList());
                    throw new UserException.IncompatibleSequenceDictionaries(String.format("Dictionary %s is missing contigs found in dictionary %s.  Missing contigs: \n %s \n", name1, name2, String.join(", ", missingContigs)), name1, dict1, name2, dict2);
                }
                return;
            case NO_COMMON_CONTIGS:
                throw new UserException.IncompatibleSequenceDictionaries("No overlapping contigs found", name1, dict1, name2, dict2);

            case UNEQUAL_COMMON_CONTIGS: {
                final List<SAMSequenceRecord> x = SAMSequenceDictionaryUtils.findDisequalCommonContigs(
                        SAMSequenceDictionaryUtils.getCommonContigsByName(dict1, dict2), dict1, dict2);
                final SAMSequenceRecord elt1 = x.get(0);
                final SAMSequenceRecord elt2 = x.get(1);
                throw new UserException.IncompatibleSequenceDictionaries(
                        String.format("Found contigs with the same name but different lengths:\n  contig %s = %s / %d\n  contig %s = %s / %d",
                        name1, elt1.getSequenceName(), elt1.getSequenceLength(),
                        name2, elt2.getSequenceName(), elt2.getSequenceLength()),
                        name1, dict1, name2, dict2
                );
            }

            case NON_CANONICAL_HUMAN_ORDER: {
                // We only get NON_CANONICAL_HUMAN_ORDER if the caller explicitly requested that we check contig ordering,
                // so we should always throw when we see it.
                final UserException ex;
                if ( SAMSequenceDictionaryUtils.nonCanonicalHumanContigOrder(dict1) ) {
                    ex = new UserException.LexicographicallySortedSequenceDictionary(name1, dict1);
                }
                else {
                    ex = new UserException.LexicographicallySortedSequenceDictionary(name2, dict2);
                }

                throw ex;
            }

            case OUT_OF_ORDER: {
                // We only get OUT_OF_ORDER if the caller explicitly requested that we check contig ordering,
                // so we should always throw when we see it.
                throw new UserException.IncompatibleSequenceDictionaries(
                                "The relative ordering of the common contigs in " + name1 + " and " + name2 +
                                " is not the same; to fix this please see: "
                                + "(https://www.broadinstitute.org/gatk/guide/article?id=1328), "
                                + " which describes reordering contigs in BAM and VCF files.",
                                name1, dict1, name2, dict2);
            }

            case DIFFERENT_INDICES: {
                // We only get DIFFERENT_INDICES if the caller explicitly requested that we check contig ordering,
                // so we should always throw when we see it.
                final String msg = "One or more contigs common to both dictionaries have " +
                        "different indices (ie., absolute positions) in each dictionary. Code " +
                        "that is sensitive to contig ordering can fail when this is the case. " +
                        "You should fix the sequence dictionaries so that all shared contigs " +
                        "occur at the same absolute positions in both dictionaries.";
                throw new UserException.IncompatibleSequenceDictionaries(msg, name1, dict1, name2, dict2);
            }
            default:
                throw new GATKException("Unexpected SequenceDictionaryComparison type: " + type);
        }
    }

}
