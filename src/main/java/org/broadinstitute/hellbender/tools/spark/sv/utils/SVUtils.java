package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Useful scraps of this and that.
 */
public final class SVUtils {

    public static final String GATKSV_CONTIG_ALIGNMENTS_READ_GROUP_ID = "GATKSVContigAlignments";

    /**
     * Returns the SINGLE sample name contained in the given {@code header}
     * @param header header supposedly holding the sample name
     * @return the sample name contained in {@code header}
     * @throws IllegalArgumentException if the {@code header} contains no or more than one sample names.
     */
    public static String getSampleId(final SAMFileHeader header) {
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        final List<String> sampleNames = readGroups.stream().map(SAMReadGroupRecord::getSample).distinct().collect(Collectors.toList());

        if (sampleNames.isEmpty())
            throw new IllegalArgumentException("Read groups contain no sample name. ");

        if (sampleNames.size() > 1)
            throw new IllegalArgumentException("Read groups must contain reads from one and only one sample, " +
                    "but we are finding the following ones in the given header: \t" + sampleNames.toString());

        return sampleNames.get(0);
    }

    /**
     * todo: this should be fixed in a new major version of htsjdk
     * this exist because for whatever reason,
     * VC.getAttributeAsStringList() sometimes returns a giant single string, while using
     * VC.getAttributeAsString() gives back an array.....
     */
    public static List<String> getAttributeAsStringList(final VariantContext vc, final String attributeKey) {
        return getAttributeAsStringStream(vc, attributeKey)
                .collect(Collectors.toList());
    }

    public static Stream<String> getAttributeAsStringStream(final VariantContext vc, final String attributeKey) {
        if ( ! vc.hasAttribute(attributeKey) )
            return Stream.empty();
        return vc.getAttributeAsStringList(attributeKey, "").stream()
                .flatMap(s -> {
                    if ( s.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR) ) {
                        return Arrays.stream( s.split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR) );
                    } else {
                        return Stream.of(s);
                    }
                });
    }

    // TODO: 9/4/18 should they be in SVInterval class?
    /**
     * Convert from 1-based closed locatable [1, e]
     *         to 1-based semi-open interval [1, e + 1).
     */
    public static SVInterval convertLocatable(@Nonnull final Locatable locatable, @Nonnull final SAMSequenceDictionary sequenceDictionary) {
        final int sequenceIndex = sequenceDictionary.getSequenceIndex(locatable.getContig());
        if (sequenceIndex < 0)
            throw new IllegalArgumentException("Provided locatable: " + locatable.toString() +
                    " doesn't seem to live on the chromosomes held in the provided dictionary: " + sequenceDictionary.toString());
        return new SVInterval(sequenceIndex,
                locatable.getStart(), locatable.getEnd() + 1);
    }

    /**
     * Opposite of {@link #convertLocatable(Locatable, SAMSequenceDictionary)}
     * Note that this method throws {@link IllegalArgumentException} if the provided SVInterval is length-0.
     */
    public static SimpleInterval convertSVInterval(@Nonnull final SVInterval svInterval,
                                                   @Nonnull final SAMSequenceDictionary sequenceDictionary) {
        final SAMSequenceRecord sequence = sequenceDictionary.getSequence(svInterval.getContig());
        if (sequence == null)
            throw new IllegalArgumentException("Provided SVInterval: " + svInterval.toString() +
                    " doesn't seem to live on the chromosomes held in the provided dictionary: " + sequenceDictionary.toString());
        if (svInterval.getLength() == 0)
            throw new IllegalArgumentException("Provided SVInterval: " + svInterval.toString() +
                    " has zero length");
        return new SimpleInterval(sequence.getSequenceName(),
                                  svInterval.getStart(),
                                  svInterval.getEnd() - 1);
    }

    /**
     * Canonical chromosomes are defined, for homo sapiens, as chromosomes 1-22, chromosome X and Y.
     * @param nonCanonicalContigNamesFile   path to file holding non-canonical chromosomes; when {@code null}, all chromosomes defined in {@code dictionary} are treated as canonical
     * @param dictionary                    dictionary assumed to hold all chromosomes of the organism, i.e. a super set of those included in {@code nonCanonicalContigNamesFile}
     */
    public static Set<String> getCanonicalChromosomes(@Nullable final String nonCanonicalContigNamesFile,
                                                      @Nonnull final SAMSequenceDictionary dictionary) {
        final Set<String> results = Utils.nonNull(dictionary).getSequences().stream().map(SAMSequenceRecord::getSequenceName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
        if (nonCanonicalContigNamesFile == null)
            return results;

        try (final Stream<String> nonCanonical = Files.lines(IOUtils.getPath(( Utils.nonNull(nonCanonicalContigNamesFile) )))) {
            nonCanonical.forEach(results::remove);
            return results;
        } catch ( final IOException ioe ) {
            throw new UserException.CouldNotReadInputFile("Can't read nonCanonicalContigNamesFile file " + nonCanonicalContigNamesFile, ioe);
        }
    }

    // =================================================================================================================

    /** return a good initialCapacity for a HashMap that will hold a given number of elements */
    public static int hashMapCapacity( final int nElements )
    {
        return (int)((nElements*4L)/3) + 1;
    }

    /** count the number of items available from an iterator */
    public static <T> int iteratorSize( final Iterator<T> itr ) {
        int result = 0;
        while ( itr.hasNext() ) { result += 1; itr.next(); }
        return result;
    }

    public static int iteratorSize( final LongIterator itr ) {
        int result = 0;
        while ( itr.hasNext() ) { result += 1; itr.next(); }
        return result;
    }

    public static <T> Iterator<T> singletonIterator( final T t ) {
        return Collections.singletonList(t).iterator();
    }

    public static Collection<SVKmer> uniquify(final Collection<SVKmer> coll1, final Collection<SVKmer> coll2) {
        Utils.nonNull(coll1, "first collection of kmers is null");
        Utils.nonNull(coll2, "second collection of kmers is null");

        final HopscotchSet<SVKmer> kmers = new HopscotchSet<>(coll1.size() + coll2.size());
        kmers.addAll(coll1);
        kmers.addAll(coll2);
        return kmers;
    }

    /** Concatenate two lists. */
    public static <T extends Object> List<T> concatenateLists(final List<T> list1, final List<T> list2) {
        Utils.validateArg(list1.getClass().equals(list2.getClass()), "Lists to be concatenated are of different classes");

        final List<T> result = new ArrayList<T>(list1.size() + list2.size());
        result.addAll(list1);
        result.addAll(list2);
        return result;
    }

    public static class IteratorFilter<T> implements Iterator<T> {
        private final Iterator<T> itr;
        private final Predicate<T> predicate;
        private T obj;

        public IteratorFilter( final Iterator<T> itr, final Predicate<T> predicate ) {
            this.itr = itr;
            this.predicate = predicate;
            advance();
        }

        @Override public boolean hasNext() { return obj != null; }

        @Override
        public T next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException("IteratorFilter is exhausted.");
            }
            final T result = obj;
            advance();
            return result;
        }

        private void advance() {
            obj = null;
            while ( itr.hasNext() ) {
                final T next = itr.next();
                if ( predicate.test(next) ) {
                    obj = next;
                    break;
                }
            }
        }
    }

    /**
     * Provides a stream collector that will collect items into an array list with a given initial capacity.
     */
    public static <T> Collector<T, ?, ArrayList<T>> arrayListCollector(final int size) {
        return Collectors.toCollection( () -> new ArrayList<>(size));
    }

    public static <T extends Enum<T>> EnumMap<T, Long> getZeroInitializedEnumMap(final Class<T> clazz) {
        EnumMap<T, Long> map = new EnumMap<>(clazz);
        for (final T t : clazz.getEnumConstants()) {
            map.put(t, 0L);
        }
        return map;
    }
    
    // =================================================================================================================

    //Workaround for seed 14695981039346656037 that doesn't fit in a signed long
    private static final long FNV64_DEFAULT_SEED = new BigInteger("14695981039346656037").longValue();

    /**
     * 64-bit FNV-1a hash for long's
     */

    public static long fnvLong64( final long toHash ) {
        return fnvLong64(FNV64_DEFAULT_SEED, toHash);
    }

    public static long fnvLong64( long start, final long toHash ) {
        final long mult = 1099511628211L;
        start ^= (toHash >> 56) & 0xffL;
        start *= mult;
        start ^= (toHash >> 48) & 0xffL;
        start *= mult;
        start ^= (toHash >> 40) & 0xffL;
        start *= mult;
        start ^= (toHash >> 32) & 0xffL;
        start *= mult;
        start ^= (toHash >> 24) & 0xffL;
        start *= mult;
        start ^= (toHash >> 16) & 0xffL;
        start *= mult;
        start ^= (toHash >> 8) & 0xffL;
        start *= mult;
        start ^= toHash & 0xffL;
        start *= mult;
        return start;
    }

    /**
     * 64-bit FNV-1a hash for byte arrays
     */
    public static long fnvByteArray64(final byte[] toHash) {
        // TODO: this is a mistake:  the constant should be the FNV64_DEFAULT_SEED, but it's the multiplier instead.
        return fnvByteArray64(1099511628211L, toHash);
    }

    public static long fnvByteArray64(long start, final byte[] toHash) {
        for (int i = 0; i < toHash.length; i += 8) {
            long val = 0;
            for (int j = 0; j < 8 && i + j < toHash.length; j++) {
                val = (val << 8) | toHash[i + j];
            }
            start = fnvLong64(start, val);
        }
        return start;
    }

}
