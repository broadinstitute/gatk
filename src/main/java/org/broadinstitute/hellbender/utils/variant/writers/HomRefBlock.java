package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * Helper class for calculating a GQ band in the GVCF writer
 *
 * A band contains GQ and DP values for a contiguous stretch of hom-ref genotypes,
 * and provides summary information about the entire block of genotypes.
 *
 * Genotypes within the HomRefBlock are restricted to hom-ref genotypes within a band of GQ scores
 */
final class HomRefBlock implements Locatable {

    private static final int HOM_REF_PL_POSITION = 0;  //the first value in the minPL[] is always the HomRef

    private final VariantContext startingVC;
    private final int minGQ, maxGQ;
    private final List<Integer> DPs = new ArrayList<>();
    private final Allele ref;
    private final int ploidy;

    private int end;
    private int[] minPLs = null;

    /**
     * Create a new HomRefBlock
     *
     * @param startingVC the VariantContext that starts this band (for starting position information)
     * @param lowerGQBound the lowerGQBound (inclusive) to use in this band
     * @param upperGQBound the upperGQBound (exclusive) to use in this band
     */
    public HomRefBlock(final VariantContext startingVC, final int lowerGQBound, final int upperGQBound, final int defaultPloidy) {
        Utils.nonNull(startingVC, "startingVC cannot be null");
        if ( lowerGQBound > upperGQBound ) { throw new IllegalArgumentException("bad lowerGQBound " + lowerGQBound + " as it's >= upperGQBound " + upperGQBound); }

        this.startingVC = startingVC;
        this.end = getStart() - 1;
        this.ref = startingVC.getReference();
        this.minGQ = lowerGQBound;
        this.maxGQ = upperGQBound;
        this.ploidy = startingVC.getMaxPloidy(defaultPloidy);
    }

    /**
     * Calculate the genotype Quality by subtracting the first smallest pl from the second smallest
     *
     * @param minPLs list of genotype likelihoods
     * @return the genotype quality based on the pls
     */
    @VisibleForTesting
    static int genotypeQualityFromPLs(final int[] minPLs){
        if (minPLs == null || minPLs.length < 3){
            throw new GATKException("minPLs must be at least size 3");
        }

        final int[] sortedPls = Arrays.copyOf(minPLs, minPLs.length);
        Arrays.sort(sortedPls);

        if (sortedPls[0] != minPLs[HOM_REF_PL_POSITION]) {
            throw new GATKException("This should be a home ref block, but the lowest pl was not for hom ref");
        }

        final int rawQuality = sortedPls[1] - sortedPls[0];

        // cap the quality to the highest quality that will be emitted by the VCFEncoder
        // this isn't strictly necessary since it will be capped when written anyway
        // it should help avoid confusion by preventing the quality from changing when written and loaded from disk
        return Math.min(rawQuality, VCFConstants.MAX_GENOTYPE_QUAL);
    }

    /**
     * Convert a HomRefBlock into a VariantContext
     *
     * @param sampleName sample name to give this variant context
     * @return a VariantContext representing the gVCF encoding for this block.
     * It will return {@code null} if input {@code block} is {@code null}, indicating that there
     * is no variant-context to be output into the VCF.
     */
    public VariantContext toVariantContext(String sampleName) {
        final VariantContextBuilder vcb = new VariantContextBuilder(getStartingVC());
        vcb.attributes(new HashMap<>(2)); // clear the attributes
        vcb.stop(getEnd());
        vcb.attribute(VCFConstants.END_KEY, getEnd());
        final Genotype genotype = createHomRefGenotype(sampleName);

        return vcb.genotypes(genotype).make();
    }

    // create a single Genotype with GQ and DP annotations
    private Genotype createHomRefGenotype(String sampleName) {
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, Collections.nCopies(getPloidy(), getRef()));
        gb.noAD().noPL().noAttributes(); // clear all attributes

        final int[] minPLs = getMinPLs();
        gb.PL(minPLs);
        gb.GQ(genotypeQualityFromPLs(minPLs));
        gb.DP(getMedianDP());
        gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());

        return gb.make();
    }

    /**
     * Add information from this Genotype to this band
     * @param genotype a non-null Genotype with GQ and DP attributes
     */
    public void add(final int pos, final Genotype genotype) {
        Utils.nonNull(genotype, "genotype cannot be null");
        if ( ! genotype.hasPL() ) { throw new IllegalArgumentException("genotype must have PL field");}
        if ( pos != end + 1 ) { throw new IllegalArgumentException("adding genotype at pos " + pos + " isn't contiguous with previous end " + end); }
        if ( genotype.getPloidy() != ploidy) { throw new IllegalArgumentException("cannot add a genotype with a different ploidy: " + genotype.getPloidy() + " != " + ploidy); }
        if ( !withinBounds(genotype.getGQ())) {
            throw new IllegalArgumentException("cannot add a genotype with GQ=" + genotype.getGQ() + " because it's not within bounds ["
                    + this.getGQLowerBound() + ',' + this.getGQUpperBound() + ')');
        }

        if( minPLs == null ) {
            minPLs = genotype.getPL();
        } else { // otherwise take the min with the provided genotype's PLs
            final int[] pls = genotype.getPL();
            if (pls.length != minPLs.length) {
                throw new GATKException("trying to merge different PL array sizes: " + pls.length + " != " + minPLs.length);
            }
            for (int i = 0; i < pls.length; i++) {
                minPLs[i] = Math.min(minPLs[i], pls[i]);
            }
        }
        end = pos;
        DPs.add(Math.max(genotype.getDP(), 0)); // DP must be >= 0
    }

    /**
     * Is the GQ value within the bounds of this GQ (GQ >= minGQ && GQ < maxGQ)
     * @param GQ the GQ value to test
     * @return true if within bounds, false otherwise
     */
    public boolean withinBounds(final int GQ) {
        return GQ >= minGQ && GQ < maxGQ;
    }

    /** Get the min DP observed within this band */
    public int getMinDP() {
        return Collections.min(DPs);
    }

    /** Get the median DP observed within this band
     * If there are an even number of DPs recorded in this band the median is the mean of the two middle values */
    public int getMedianDP() {
        return (int) Math.round(MathUtils.median(DPs));
    }

    /** Get the min PLs observed within this band, can be null if no PLs have yet been observed */
    public int[] getMinPLs() {
        return minPLs;
    }

    int getGQUpperBound() {
        return maxGQ;
    }
    int getGQLowerBound() {
        return minGQ;
    }

    public boolean isContiguous(final VariantContext vc) {
        return (vc.getEnd() == getEnd() + 1) && startingVC.getContig().equals(vc.getContig());
    }

    public VariantContext getStartingVC() {
        return startingVC;
    }

    @Override
    public String getContig() {
        return startingVC.getContig();
    }

    @Override
    public int getStart() {
        return startingVC.getStart();
    }

    @Override
    public int getEnd() {
        return end;
    }

    public Allele getRef() {
        return ref;
    }

    public int getSize() {
        return getEnd() - getStart() + 1;
    }

    @Override
    public String toString() {
        return "HomRefBlock{" +
                "minGQ=" + minGQ +
                ", maxGQ=" + maxGQ +
                '}';
    }

    /**
    * @return the ploidy of this hom-ref block.
     */
    public int getPloidy() {
        return ploidy;
    }
}
