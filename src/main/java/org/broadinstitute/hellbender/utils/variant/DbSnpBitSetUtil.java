package org.broadinstitute.hellbender.utils.variant;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;

/**
 * Utility class to use with DbSnp files to determine is a locus is
 * a dbSnp site.
 */
public final class DbSnpBitSetUtil {

    private final Map<String, BitSet> sequenceToBitSet = new LinkedHashMap<>();

    /** Little tuple class to contain one bitset for SNPs and another for Indels. */
    public static class DbSnpBitSets {
        public DbSnpBitSetUtil snps;
        public DbSnpBitSetUtil indels;
    }

    /** Private empty contructor for use by factory methods only. */
    private DbSnpBitSetUtil() { }

    /** Constructor that creates a bit set with bits set to true for all variant types. */
    public DbSnpBitSetUtil(final File dbSnpFile, final SAMSequenceDictionary sequenceDictionary) {
        this(dbSnpFile, sequenceDictionary, EnumSet.noneOf(DbSnpVariantType.class));
    }

    /**
     * Constructor.
     *
     * For each sequence, creates  a BitSet that denotes whether a dbSNP entry
     * is present at each base in the reference sequence.  The set is
     * reference.length() + 1 so that it can be indexed by 1-based reference base.
     * True means dbSNP present, false means no dbSNP present.
     *
     * @param dbSnpFile in VCF format.
     * @param sequenceDictionary Optionally, a sequence dictionary corresponding to the dbSnp file, else null.
     * If present, BitSets will be allocated more efficiently because the maximum size will be known.
     * @param variantsToMatch what types of variants to load.
     */
    public DbSnpBitSetUtil(final File dbSnpFile,
                           final SAMSequenceDictionary sequenceDictionary,
                           final Collection<DbSnpVariantType> variantsToMatch) {
        Utils.nonNull(dbSnpFile);
        final Map<DbSnpBitSetUtil, Set<DbSnpVariantType>> tmp = new LinkedHashMap<>();
        tmp.put(this, EnumSet.copyOf(variantsToMatch));
        loadVcf(dbSnpFile, sequenceDictionary, tmp);
    }

    /** Private helper method to read through the VCF and create one or more bit sets. */
    private static void loadVcf(final File dbSnpFile,
                                final SAMSequenceDictionary sequenceDictionary,
                                final Map<DbSnpBitSetUtil, Set<DbSnpVariantType>> bitSetsToVariantTypes) {

        final VCFFileReader variantReader = new VCFFileReader(dbSnpFile);
        final CloseableIterator<VariantContext> variantIterator = variantReader.iterator();

        while (variantIterator.hasNext()) {
            final VariantContext kv = variantIterator.next();

            for (final Map.Entry<DbSnpBitSetUtil, Set<DbSnpVariantType>> tuple : bitSetsToVariantTypes.entrySet()) {
                final DbSnpBitSetUtil bitset            = tuple.getKey();
                final Set<DbSnpVariantType> variantsToMatch  = tuple.getValue();

                BitSet bits = bitset.sequenceToBitSet.get(kv.getContig());
                if (bits == null) {
                    final int nBits;
                    if (sequenceDictionary == null) nBits = kv.getEnd() + 1;
                    else nBits = sequenceDictionary.getSequence(kv.getContig()).getSequenceLength() + 1;
                    bits = new BitSet(nBits);
                    bitset.sequenceToBitSet.put(kv.getContig(), bits);
                }
                if (variantsToMatch.isEmpty() ||
                        (kv.isSNP() && variantsToMatch.contains(DbSnpVariantType.SNP)) ||
                        (kv.isIndel() && variantsToMatch.contains(DbSnpVariantType.insertion)) ||
                        (kv.isIndel() && variantsToMatch.contains(DbSnpVariantType.deletion))) {

                    for (int i = kv.getStart(); i <= kv.getEnd(); i++)  bits.set(i, true);
                }
            }
        }

        CloserUtil.close(variantIterator);
        CloserUtil.close(variantReader);
    }

    /**
     * Returns true if there is a dbSnp entry at pos in sequenceName, otherwise false
     */
    public boolean isDbSnpSite(final String sequenceName, final int pos) {
        // When we have a dbSnpFile with no sequence dictionary, this line will be necessary
        if (sequenceToBitSet.get(sequenceName) == null) {
            return false;
        }
        if (pos > sequenceToBitSet.get(sequenceName).length()) {
            return false;
        }
        return sequenceToBitSet.get(sequenceName).get(pos);
    }

}
