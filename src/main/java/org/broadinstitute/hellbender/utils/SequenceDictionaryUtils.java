package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.*;

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
public class SequenceDictionaryUtils {
    //
    // for detecting lexicographically sorted human references
    //
    private static final boolean ENABLE_LEXICOGRAPHIC_REQUIREMENT_FOR_HUMAN = true;

    // The following sets of contig records are used to perform the non-canonical human ordering check.
    // This check ensures that the order is 1,2,3... instead of 1, 10, 11, 12...2, 20, 21...

    // hg18
    protected static final SAMSequenceRecord CHR1_HG18 = new SAMSequenceRecord("chr1", 247249719);
    protected static final SAMSequenceRecord CHR2_HG18 = new SAMSequenceRecord("chr2", 242951149);
    protected static final SAMSequenceRecord CHR10_HG18 = new SAMSequenceRecord("chr10", 135374737);

    // hg19
    protected static final SAMSequenceRecord CHR1_HG19 = new SAMSequenceRecord("chr1", 249250621);
    protected static final SAMSequenceRecord CHR2_HG19 = new SAMSequenceRecord("chr2", 243199373);
    protected static final SAMSequenceRecord CHR10_HG19 = new SAMSequenceRecord("chr10", 135534747);

    // b36
    protected static final SAMSequenceRecord CHR1_B36 = new SAMSequenceRecord("1", 247249719);
    protected static final SAMSequenceRecord CHR2_B36 = new SAMSequenceRecord("2", 242951149);
    protected static final SAMSequenceRecord CHR10_B36 = new SAMSequenceRecord("10", 135374737);

    // b37
    protected static final SAMSequenceRecord CHR1_B37 = new SAMSequenceRecord("1", 249250621);
    protected static final SAMSequenceRecord CHR2_B37 = new SAMSequenceRecord("2", 243199373);
    protected static final SAMSequenceRecord CHR10_B37 = new SAMSequenceRecord("10", 135534747);


    public enum SequenceDictionaryCompatibility {
        IDENTICAL,                      // the dictionaries are identical
        COMMON_SUBSET,                  // there exists a common subset of equivalent contigs
        NO_COMMON_CONTIGS,              // no overlap between dictionaries
        UNEQUAL_COMMON_CONTIGS,         // common subset has contigs that have the same name but different lengths
        NON_CANONICAL_HUMAN_ORDER,      // human reference detected but the order of the contigs is non-standard (lexicographic, for examine)
        OUT_OF_ORDER,                   // the two dictionaries overlap but the overlapping contigs occur in different
        // orders with respect to each other
        DIFFERENT_INDICES               // the two dictionaries overlap and the overlapping contigs occur in the same
        // order with respect to each other, but one or more of them have different
        // indices in the two dictionaries. Eg., { chrM, chr1, chr2 } vs. { chr1, chr2 }
    }

    /**
     * Tests for compatibility between two sequence dictionaries.  If the dictionaries are incompatible, then
     * UserExceptions are thrown with detailed error messages.
     * Two sequence dictionaries are compatible if they have the same names and lengths for contigs, or share some contigs.
     *
     * @param name1 name associated with dict1
     * @param dict1 the sequence dictionary dict1
     * @param name2 name associated with dict2
     * @param dict2 the sequence dictionary dict2
     * @param isReadsToReferenceComparison true if one of the dictionaries comes from a reads data source (eg., a BAM),
     *                                     and the other from a reference data source
     * @param intervals the user-specified genomic intervals: only required when isReadsToReferenceComparison is true,
     *                  otherwise can be null
     */
    public static void validateDictionaries( final String name1,
                                             final SAMSequenceDictionary dict1,
                                             final String name2,
                                             final SAMSequenceDictionary dict2,
                                             final boolean isReadsToReferenceComparison,
                                             final List<SimpleInterval> intervals ) {

        final SequenceDictionaryCompatibility type = compareDictionaries(dict1, dict2);

        switch ( type ) {
            case IDENTICAL:
                return;
            case COMMON_SUBSET:
                return;
            case NO_COMMON_CONTIGS:
                throw new UserException.IncompatibleSequenceDictionaries("No overlapping contigs found", name1, dict1, name2, dict2);

            case UNEQUAL_COMMON_CONTIGS: {
                List<SAMSequenceRecord> x = findDisequalCommonContigs(getCommonContigsByName(dict1, dict2), dict1, dict2);
                SAMSequenceRecord elt1 = x.get(0);
                SAMSequenceRecord elt2 = x.get(1);
                UserException ex = new UserException.IncompatibleSequenceDictionaries(String.format("Found contigs with the same name but different lengths:\n  contig %s = %s / %d\n  contig %s = %s / %d",
                        name1, elt1.getSequenceName(), elt1.getSequenceLength(),
                        name2, elt2.getSequenceName(), elt2.getSequenceLength()),
                        name1, dict1, name2, dict2);

                throw ex;
            }

            case NON_CANONICAL_HUMAN_ORDER: {
                UserException ex;
                if ( nonCanonicalHumanContigOrder(dict1) ) {
                    ex = new UserException.LexicographicallySortedSequenceDictionary(name1, dict1);
                }
                else {
                    ex = new UserException.LexicographicallySortedSequenceDictionary(name2, dict2);
                }

                throw ex;
            }

            case OUT_OF_ORDER: {
                UserException ex = new UserException.IncompatibleSequenceDictionaries(
                        "The relative ordering of the common contigs in " + name1 + " and " + name2 +
                                " is not the same; to fix this please see: "
                                + "(https://www.broadinstitute.org/gatk/guide/article?id=1328), "
                                + " which describes reordering contigs in BAM and VCF files.",
                        name1, dict1, name2, dict2);
                throw ex;
            }

            case DIFFERENT_INDICES: {
                // This is currently only known to be problematic when the index mismatch is between a bam and the
                // reference AND when the user's intervals actually include one or more of the contigs that are
                // indexed differently from the reference. In this case, the engine will fail to correctly serve
                // up the reads from those contigs, so throw an exception.
                if ( isReadsToReferenceComparison && intervals != null ) {

                    final Set<String> misindexedContigs = findMisindexedContigsInIntervals(intervals, dict1, dict2);

                    if ( ! misindexedContigs.isEmpty() ) {
                        final String msg = "The following contigs included in the intervals to process have " +
                                        "different indices in the sequence dictionaries for the reads vs. " +
                                        "the reference: %s.  As a result, the GATK engine will not correctly " +
                                        "process reads from these contigs. You should either fix the sequence " +
                                        "dictionaries for your reads so that these contigs have the same indices " +
                                        "as in the sequence dictionary for your reference, or exclude these contigs " +
                                        "from your intervals.";
                        final UserException ex = new UserException.IncompatibleSequenceDictionaries(msg, name1, dict1, name2, dict2);
                        throw ex;
                    }
                }
                break;
            }
            default:
                throw new GATKException("Unexpected SequenceDictionaryComparison type: " + type);
        }
    }

    /**
     * Workhorse routine that takes two dictionaries and returns their compatibility.
     *
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return A SequenceDictionaryCompatibility enum value describing the compatibility of the two dictionaries
     */
    public static SequenceDictionaryCompatibility compareDictionaries( final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2) {
        if ( nonCanonicalHumanContigOrder(dict1) || nonCanonicalHumanContigOrder(dict2) )
            return SequenceDictionaryCompatibility.NON_CANONICAL_HUMAN_ORDER;

        final Set<String> commonContigs = getCommonContigsByName(dict1, dict2);

        if (commonContigs.size() == 0)
            return SequenceDictionaryCompatibility.NO_COMMON_CONTIGS;
        else if ( ! commonContigsHaveSameLengths(commonContigs, dict1, dict2) )
            return SequenceDictionaryCompatibility.UNEQUAL_COMMON_CONTIGS;
        else if ( ! commonContigsAreInSameRelativeOrder(commonContigs, dict1, dict2) )
            return SequenceDictionaryCompatibility.OUT_OF_ORDER;
        else if ( commonContigs.size() == dict1.size() && commonContigs.size() == dict2.size() )
            return SequenceDictionaryCompatibility.IDENTICAL;
        else if ( ! commonContigsAreAtSameIndices(commonContigs, dict1, dict2) )
            return SequenceDictionaryCompatibility.DIFFERENT_INDICES;
        else {
            return SequenceDictionaryCompatibility.COMMON_SUBSET;
        }
    }

    /**
     * Utility function that tests whether the commonContigs in both dicts are equivalent.  Equivalence means
     * that the seq records have the same length, if both are non-zero.
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return true if all of the common contigs are equivalent
     */
    private static boolean commonContigsHaveSameLengths(Set<String> commonContigs, SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        return findDisequalCommonContigs(commonContigs, dict1, dict2) == null;
    }

    /**
     * Returns a List(x,y) that contains two disequal sequence records among the common contigs in both dicts.  Returns
     * null if all common contigs are equivalent
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return
     */
    private static List<SAMSequenceRecord> findDisequalCommonContigs(Set<String> commonContigs, SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        for ( String name : commonContigs ) {
            SAMSequenceRecord elt1 = dict1.getSequence(name);
            SAMSequenceRecord elt2 = dict2.getSequence(name);
            if ( ! sequenceRecordsAreEquivalent(elt1, elt2) )
                return Arrays.asList(elt1,elt2);
        }

        return null;
    }

    /**
     * Helper routine that returns two sequence records are equivalent, defined as having the same name and
     * lengths, if both are non-zero
     *
     * @param me
     * @param that
     * @return
     */
    private static boolean sequenceRecordsAreEquivalent(final SAMSequenceRecord me, final SAMSequenceRecord that) {
        if (me == that) return true;
        if (that == null) return false;

        if (me.getSequenceLength() != 0 && that.getSequenceLength() != 0 && me.getSequenceLength() != that.getSequenceLength())
            return false;

        // todo -- reenable if we want to be really strict here
//        if (me.getExtendedAttribute(SAMSequenceRecord.MD5_TAG) != null && that.getExtendedAttribute(SAMSequenceRecord.MD5_TAG) != null) {
//            final BigInteger thisMd5 = new BigInteger((String)me.getExtendedAttribute(SAMSequenceRecord.MD5_TAG), 16);
//            final BigInteger thatMd5 = new BigInteger((String)that.getExtendedAttribute(SAMSequenceRecord.MD5_TAG), 16);
//            if (!thisMd5.equals(thatMd5)) {
//                return false;
//            }
//        }
//        else {
        if (me.getSequenceName() != that.getSequenceName())
            return false; // Compare using == since we intern() the Strings
//        }

        return true;
    }

    /**
     * A very simple (and naive) algorithm to determine (1) if the dict is a human reference (hg18, hg19, b36, or b37) and if it's
     * lexicographically sorted.  Works by matching lengths of the static chr1, chr10, and chr2, and then if these
     * are all matched, requiring that the order be chr1, chr2, chr10.
     *
     * @param dict
     * @return
     */
    private static boolean nonCanonicalHumanContigOrder(SAMSequenceDictionary dict) {
        if ( ! ENABLE_LEXICOGRAPHIC_REQUIREMENT_FOR_HUMAN ) // if we don't want to enable this test, just return false
            return false;

        SAMSequenceRecord chr1 = null, chr2 = null, chr10 = null;
        for ( SAMSequenceRecord elt : dict.getSequences() ) {
            if ( isHumanSeqRecord(elt, CHR1_HG18, CHR1_HG19, CHR1_B36, CHR1_B37) ) chr1 = elt;
            if ( isHumanSeqRecord(elt, CHR2_HG18, CHR2_HG19, CHR2_B36, CHR2_B37) ) chr2 = elt;
            if ( isHumanSeqRecord(elt, CHR10_HG18, CHR10_HG19, CHR10_B36, CHR10_B37) ) chr10 = elt;
        }
        if ( chr1 != null  && chr2 != null && chr10 != null) {
            return ! ( chr1.getSequenceIndex() < chr2.getSequenceIndex() && chr2.getSequenceIndex() < chr10.getSequenceIndex() );
        }

        return false;
    }

    /**
     * Trivial helper that returns true if elt has the same name and length as rec1 or rec2
     * @param elt record to test
     * @param recs the list of records to check for name and length equivalence
     * @return true if elt has the same name and length as any of the recs
     */
    private static boolean isHumanSeqRecord(SAMSequenceRecord elt, SAMSequenceRecord... recs) {
        for (SAMSequenceRecord rec : recs) {
            if (elt.getSequenceLength() == rec.getSequenceLength() && elt.getSequenceName().equals(rec.getSequenceName())) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns true if the common contigs in dict1 and dict2 are in the same relative order, without regard to
     * absolute index position. This is accomplished by getting the common contigs in both dictionaries, sorting
     * these according to their indices, and then walking through the sorted list to ensure that each ordered contig
     * is equivalent
     *
     * @param commonContigs names of the contigs common to both dictionaries
     * @param dict1 first SAMSequenceDictionary
     * @param dict2 second SAMSequenceDictionary
     * @return true if the common contigs occur in the same relative order in both dict1 and dict2, otherwise false
     */
    private static boolean commonContigsAreInSameRelativeOrder(Set<String> commonContigs, SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        List<SAMSequenceRecord> list1 = sortSequenceListByIndex(getSequencesOfName(commonContigs, dict1));
        List<SAMSequenceRecord> list2 = sortSequenceListByIndex(getSequencesOfName(commonContigs, dict2));

        for ( int i = 0; i < list1.size(); i++ ) {
            SAMSequenceRecord elt1 = list1.get(i);
            SAMSequenceRecord elt2 = list2.get(i);
            if ( ! elt1.getSequenceName().equals(elt2.getSequenceName()) )
                return false;
        }

        return true;
    }

    /**
     * Gets the subset of SAMSequenceRecords in commonContigs in dict
     *
     * @param commonContigs
     * @param dict
     * @return
     */
    private static List<SAMSequenceRecord> getSequencesOfName(Set<String> commonContigs, SAMSequenceDictionary dict) {
        List<SAMSequenceRecord> l = new ArrayList<SAMSequenceRecord>(commonContigs.size());
        for ( String name : commonContigs ) {
            l.add(dict.getSequence(name) );
        }

        return l;
    }

    /**
     * Compares sequence records by their order
     */
    private static class CompareSequenceRecordsByIndex implements Comparator<SAMSequenceRecord> {
        public int compare(SAMSequenceRecord x, SAMSequenceRecord y) {
            return Integer.valueOf(x.getSequenceIndex()).compareTo(y.getSequenceIndex());
        }
    }

    /**
     * Returns a sorted list of SAMSequenceRecords sorted by their indices.  Note that the
     * list is modified in place, so the returned list is == to the unsorted list.
     *
     * @param unsorted
     * @return
     */
    private static List<SAMSequenceRecord> sortSequenceListByIndex(List<SAMSequenceRecord> unsorted) {
        Collections.sort(unsorted, new CompareSequenceRecordsByIndex());
        return unsorted;
    }

    /**
     * Checks whether the common contigs in the given sequence dictionaries occur at the same indices
     * in both dictionaries
     *
     * @param commonContigs Set of names of the contigs that occur in both dictionaries
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return true if the contigs common to dict1 and dict2 occur at the same indices in both dictionaries,
     *         otherwise false
     */
    private static boolean commonContigsAreAtSameIndices( final Set<String> commonContigs, final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2 ) {
        for ( String commonContig : commonContigs ) {
            SAMSequenceRecord dict1Record = dict1.getSequence(commonContig);
            SAMSequenceRecord dict2Record = dict2.getSequence(commonContig);

            // Each common contig must have the same index in both dictionaries
            if ( dict1Record.getSequenceIndex() != dict2Record.getSequenceIndex() ) {
                return false;
            }
        }

        return true;
    }

    /**
     * Gets the set of names of the contigs found in both sequence dictionaries that have different indices
     * in the two dictionaries.
     *
     * @param commonContigs Set of names of the contigs common to both dictionaries
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return a Set containing the names of the common contigs indexed differently in dict1 vs. dict2,
     *         or an empty Set if there are no such contigs
     */
    private static Set<String> getDifferentlyIndexedCommonContigs( final Set<String> commonContigs,
                                                                   final SAMSequenceDictionary dict1,
                                                                   final SAMSequenceDictionary dict2 ) {

        final Set<String> differentlyIndexedCommonContigs = new LinkedHashSet<String>(Utils.optimumHashSize(commonContigs.size()));

        for ( String commonContig : commonContigs ) {
            if ( dict1.getSequence(commonContig).getSequenceIndex() != dict2.getSequence(commonContig).getSequenceIndex() ) {
                differentlyIndexedCommonContigs.add(commonContig);
            }
        }

        return differentlyIndexedCommonContigs;
    }

    /**
     * Finds the names of any contigs indexed differently in the two sequence dictionaries that also
     * occur in the provided set of intervals.
     *
     * @param intervals GenomeLocSortedSet containing the intervals to check
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return a Set of the names of the contigs indexed differently in dict1 vs dict2 that also
     *         occur in the provided intervals, or an empty Set if there are no such contigs
     */
    private static Set<String> findMisindexedContigsInIntervals( final List<SimpleInterval> intervals,
                                                                 final SAMSequenceDictionary dict1,
                                                                 final SAMSequenceDictionary dict2 ) {

        final Set<String> differentlyIndexedCommonContigs = getDifferentlyIndexedCommonContigs(getCommonContigsByName(dict1, dict2), dict1, dict2);
        final Set<String> misindexedContigsInIntervals = new LinkedHashSet<String>(Utils.optimumHashSize(differentlyIndexedCommonContigs.size()));

        // We know differentlyIndexedCommonContigs is a HashSet, so this loop is O(intervals)
        for ( SimpleInterval interval : intervals ) {
            if ( differentlyIndexedCommonContigs.contains(interval.getContig()) ) {
                misindexedContigsInIntervals.add(interval.getContig());
            }
        }

        return misindexedContigsInIntervals;
    }

    /**
     * Returns the set of contig names found in both dicts.
     * @param dict1
     * @param dict2
     * @return
     */
    public static Set<String> getCommonContigsByName(SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        Set<String> intersectingSequenceNames = getContigNames(dict1);
        intersectingSequenceNames.retainAll(getContigNames(dict2));
        return intersectingSequenceNames;
    }

    public static Set<String> getContigNames(SAMSequenceDictionary dict) {
        Set<String> contigNames = new HashSet<String>(Utils.optimumHashSize(dict.size()));
        for (SAMSequenceRecord dictionaryEntry : dict.getSequences())
            contigNames.add(dictionaryEntry.getSequenceName());
        return contigNames;
    }

    /**
     * Returns a compact String representation of the sequence dictionary it's passed
     *
     * The format of the returned String is:
     * [ contig1Name(length: contig1Length) contig2Name(length: contig2Length) ... ]
     *
     * @param dict a non-null SAMSequenceDictionary
     * @return A String containing all of the contig names and lengths from the sequence dictionary it's passed
     */
    public static String getDictionaryAsString( final SAMSequenceDictionary dict ) {
        if ( dict == null ) {
            throw new IllegalArgumentException("Sequence dictionary must be non-null");
        }

        StringBuilder s = new StringBuilder("[ ");

        for ( SAMSequenceRecord dictionaryEntry : dict.getSequences() ) {
            s.append(dictionaryEntry.getSequenceName());
            s.append("(length:");
            s.append(dictionaryEntry.getSequenceLength());
            s.append(") ");
        }

        s.append("]");

        return s.toString();
    }
}