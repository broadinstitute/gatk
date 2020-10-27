package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;

/**
 * Helper class for calculating a somatic LOD band in the SomaticGVCF writer
 *
 * A band contains LOD and DP values for a contiguous stretch of hom-ref genotypes,
 * and provides summary information about the entire block of genotypes.
 *
 * Genotypes within the TLODBlock are restricted to hom-ref genotypes within a band of LOD scores
 *
 * LODs are stored as ints to facilitate GVCF class polymorphism and to avoid double precision issues
 */
final class TLODBlock extends GVCFBlock {

    private double minBlockLOD = Double.POSITIVE_INFINITY;

    //effectively the number of decimal points to use for equality comparisons
    private int partitionPrecision;

    /**
     * Create a new HomRefBlock
     *
     * @param startingVC the VariantContext that starts this band (for starting position information)
     * @param lowerLODBoundAsBinnedInt the lower LOD Bound (inclusive) to use in this band, represented as a binned integer, i.e. multiplied by 10^partitionPrecision
     * @param upperLODBoundAsBinnedInt the upper LOD Bound (exclusive) to use in this band
     */
    TLODBlock(final VariantContext startingVC, final int lowerLODBoundAsBinnedInt, final int upperLODBoundAsBinnedInt, final int partitionPrecision) {
        super(startingVC, (int)Math.floor(lowerLODBoundAsBinnedInt), (int)Math.floor(upperLODBoundAsBinnedInt));
        this.partitionPrecision = partitionPrecision;
    }

    private int convertLODtoInt(final double LOD) {
        return (int)Math.floor(LOD * Math.pow(10, partitionPrecision));
    }

    private double convertIntToLOD(final int binnedValue) {
        return (double)binnedValue / Math.pow(10, partitionPrecision);
    }

    /** Get the min TLOD observed within this band, can be null if no TLODs have yet been observed */
    public double getMinBlockLOD() {
        return minBlockLOD;
    }

    public double getLODLowerBound() {
        return convertIntToLOD(getGQLowerBound());
    }

    public double getLODUpperBound() {
        return convertIntToLOD(getGQUpperBound());
    }

    boolean withinBounds(final double lod) {
        return withinBounds(convertLODtoInt(lod));
    }

    // create a single Genotype with GQ and DP annotations
    @Override
    Genotype createHomRefGenotype(final String sampleName, final boolean floorBlock) {
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, Collections.nCopies(2, getRef()));  //FIXME: for somatic stuff we output the genotype as diploid because that's familiar for human
        gb.noAD().noPL().noAttributes(); // clear all attributes

        gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, minBlockLOD);
        gb.DP(getMedianDP());
        gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());

        return gb.make();
    }

    /**
     * Add information from this Genotype to this band.
     *
     * @param pos Current genomic position. Must be 1 base after the previous position
     * @param genotype A non-null Genotype with TLOD and DP attributes
     */
    @Override
    public void add(final int pos, final int newEnd, final Genotype genotype) {
        Utils.nonNull(genotype, "genotype cannot be null");
        if ( pos > end + 1 ) { throw new IllegalArgumentException("adding genotype at pos " + pos + " isn't contiguous with previous end " + end); }
        // Make sure the LOD is within the bounds of this band
        final double currentLOD = Double.parseDouble(genotype.getExtendedAttribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY).toString());
        if ( !withinBounds(currentLOD)) {
            throw new IllegalArgumentException("cannot add a genotype with LOD=" + currentLOD + " because it's not within bounds ["
                    + this.getLODLowerBound() + ',' + this.getLODUpperBound() + ')');
        }

        if( minBlockLOD == Double.POSITIVE_INFINITY || currentLOD < minBlockLOD) {
            minBlockLOD = currentLOD;
        }

        end = newEnd;
        DPs.add(Math.max(genotype.getDP(), 0)); // DP must be >= 0
    }

    @Override
    public String toString() {
        return "TLODBlock{" +
                "minLOD=" + getLODLowerBound() +
                ", maxLOD=" + getLODUpperBound() +
                '}';
    }
}

