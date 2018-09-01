package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import javax.annotation.Nonnull;
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

    public static String getSampleId(final SAMFileHeader header) {
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        final Set<String> sampleSet = readGroups.stream().map(SAMReadGroupRecord::getSample).collect(Collectors.toSet());

        Utils.validate(sampleSet.size() == 1,
                "Read groups must contain reads from one and only one sample, " +
                        "but we are finding the following ones in the given header: \t" + sampleSet.toString());

        final String sample = sampleSet.iterator().next();
        return sample;
    }

    /**
     * Given {@code sortOrder}, provide appropriate comparator.
     * Currently only support coordinate or query-name order,
     * and throws UserException if other values are specified.
     */
    public static SAMRecordComparator getSamRecordComparator(final SAMFileHeader.SortOrder sortOrder) {
        final SAMRecordComparator samRecordComparator;
        switch (sortOrder) {
            case coordinate:
                samRecordComparator = new SAMRecordCoordinateComparator();
                break;
            case queryname:
                samRecordComparator = new SAMRecordQueryNameComparator();
                break;
            default:
                throw new UserException("Unsupported assembly alignment sort order specified");
        }
        return samRecordComparator;
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

    /**
     * Given chromosome name and requested position (1-based coordinate system),
     * return a 1 bp long interval that spans only the requested base.
     */
    public static SimpleInterval makeOneBpInterval(final String chr, final int pos) {
        return new SimpleInterval(chr, pos, pos);
    }

    /**
     * Canonical chromosomes are defined, for homo sapiens, as chromosomes 1-22, chromosome X and Y.
     */
    public static Set<String> getCanonicalChromosomes(final String nonCanonicalContigNamesFile,
                                                      @Nonnull final SAMSequenceDictionary dictionary) {
        final LinkedHashSet<String> allContigs = Utils.nonNull(dictionary).getSequences().stream().map(SAMSequenceRecord::getSequenceName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
        if (nonCanonicalContigNamesFile == null)
            return allContigs;

        try (final Stream<String> nonCanonical = Files.lines(IOUtils.getPath(( Utils.nonNull(nonCanonicalContigNamesFile) )))) {
            nonCanonical.forEach(allContigs::remove);
            return allContigs;
        } catch ( final IOException ioe ) {
            throw new UserException("Can't read nonCanonicalContigNamesFile file "+nonCanonicalContigNamesFile, ioe);
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
