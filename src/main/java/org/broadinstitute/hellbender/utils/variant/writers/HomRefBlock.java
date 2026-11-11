package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.variantutils.PosteriorProbabilitiesUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Helper class for calculating a GQ band in the GVCF writer
 *
 * A band contains GQ and DP values for a contiguous stretch of hom-ref genotypes,
 * and provides summary information about the entire block of genotypes.
 *
 * Genotypes within the HomRefBlock are restricted to hom-ref genotypes within a band of GQ scores
 */
final class HomRefBlock extends GVCFBlock {

    private static final int HOM_REF_PL_POSITION = 0;  //the first value in the minPL[] is always the HomRef

    private final int ploidy;

    // The block lowers these in place as genotypes are added, so it must own them: it copies every array it adopts. A
    // genotype's PL array is often shared with other genotypes (GenotypeBuilder(Genotype) keeps the same array), and
    // PosteriorProbabilitiesUtils.parsePosteriorsIntoPhredSpace returns the genotype's own PL array when the posteriors
    // are missing. Otherwise it returns a new array, and the copy costs one small array per block.
    private int[] minPLs = null;
    private int[] minPPs = null;
    private int minGQ = -1;

    /**
     * Create a new HomRefBlock
     *
     * @param startingVC the VariantContext that starts this band (for starting position information)
     * @param lowerGQBound the lowerGQBound (inclusive) to use in this band
     * @param upperGQBound the upperGQBound (exclusive) to use in this band
     */
    public HomRefBlock(final VariantContext startingVC, final int lowerGQBound, final int upperGQBound, final int defaultPloidy) {
        super(startingVC, lowerGQBound, upperGQBound);
        Utils.nonNull(startingVC, "startingVC cannot be null");
        Utils.validateArg(upperGQBound <= VCFConstants.MAX_GENOTYPE_QUAL + 1, "upperGQBound must be <= " + (VCFConstants.MAX_GENOTYPE_QUAL + 1));
        if ( lowerGQBound > upperGQBound ) { throw new IllegalArgumentException("bad lowerGQBound " + lowerGQBound + " as it's >= upperGQBound " + upperGQBound); }

        this.ploidy = startingVC.getMaxPloidy(defaultPloidy);
        final Genotype g = startingVC.getGenotype(0);
        if (g.hasPL()) {
            this.minPLs = ArrayUtils.clone(g.getPL());
        }
        if (g.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY)) {
            this.minPPs = ArrayUtils.clone(PosteriorProbabilitiesUtils.parsePosteriorsIntoPhredSpace(g));
        }
        if (minPPs != null) {
            minGQ = GATKVariantContextUtils.calculateGQFromPLs(minPPs);
        } else if (minPLs != null) {
            minGQ = GATKVariantContextUtils.calculateGQFromPLs(minPLs);
        } else if (g.hasGQ()) {
            minGQ = g.getGQ();
        }
    }



    // create a single Genotype with GQ and DP annotations
    @Override
    Genotype createHomRefGenotype(final String sampleName, final boolean floorBlocks) {
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, Collections.nCopies(getPloidy(), getRef()));
        gb.noAD().noPL().noAttributes(); // clear all attributes

        final int[] minPLs = getMinPLs();
        final int[] minPPs = getMinPPs();
        if (!floorBlocks) {
            if (minPLs != null) {
                gb.PL(minPLs);
            }
            gb.GQ(getMinGQ());
            if (!DPs.isEmpty()) {
                gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());
            }
        }
        else {
            gb.GQ(getGQLowerBound());
        }
        if (!DPs.isEmpty()) {
            gb.DP(getMedianDP());
        }
        if (minPPs != null) {
            gb.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY, Utils.listFromPrimitives(minPPs));
        }

        return gb.make();
    }

    /**
     * Add a homRef block to the current block
     *
     * @param pos current genomic position, must be directly following this block (i.e. contiguous and non-overlapping)
     * @param newEnd new calculated block end position
     * @param genotype A non-null Genotype with GQ and DP attributes
     */
    @Override
    public void add(final int pos, final int newEnd, final Genotype genotype) {
        Utils.nonNull(genotype, "genotype cannot be null");
        if ( pos > end + 1 ) { throw new IllegalArgumentException("adding genotype at pos " + pos + " isn't contiguous with previous end " + end); }
        if ( pos < end + 1 ) { throw new IllegalArgumentException("adding genotype at pos " + pos + " overlaps previous end " + end); }
        if ( genotype.getPloidy() != ploidy) { throw new IllegalArgumentException("cannot add a genotype with a different ploidy: " + genotype.getPloidy() + " != " + ploidy); }
        // Make sure the GQ is within the bounds of this band. Treat GQs > 99 as 99.
        if ( !withinBounds(Math.min(genotype.getGQ(), VCFConstants.MAX_GENOTYPE_QUAL))) {
            throw new IllegalArgumentException("cannot add a genotype with GQ=" + genotype.getGQ() + " because it's not within bounds ["
                    + this.getGQLowerBound() + ',' + this.getGQUpperBound() + ')');
        }

        // A genotype without posteriors counts toward them by its PLs, as parsePosteriorsIntoPhredSpace reads missing
        // posteriors. When the first posteriors arrive, the genotypes before them count by their PLs: the minimum PLs
        // from before this genotype's are merged in. Genotypes with neither PLs nor posteriors do not count.
        final boolean posteriorsStartHere = minPPs == null && genotype.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
        final int[] earlierMinPLs = posteriorsStartHere ? ArrayUtils.clone(minPLs) : null;

        if( minPLs == null ) {
            minPLs = ArrayUtils.clone(genotype.getPL());
        }
        else if (genotype.hasPL()) { // otherwise take the min with the provided genotype's PLs
            lowerToMin(minPLs, genotype.getPL(), "PL");
        }

        if (minPPs != null || posteriorsStartHere) {
            final int[] pps = PosteriorProbabilitiesUtils.parsePosteriorsIntoPhredSpace(genotype);
            if (posteriorsStartHere) {
                minPPs = ArrayUtils.clone(pps);
                if (minPPs != null && earlierMinPLs != null) {
                    lowerToMin(minPPs, earlierMinPLs, "PL and PP");
                }
            } else if (pps != null) { // otherwise take the min with the provided genotype's posteriors
                lowerToMin(minPPs, pps, "PP");
            }
        }
        if (minPPs != null) {
            minGQ = GATKVariantContextUtils.calculateGQFromPLs(minPPs);
        } else if (minPLs != null) {
            minGQ = GATKVariantContextUtils.calculateGQFromPLs(minPLs);
        } else {
            minGQ = minGQ == -1 ? genotype.getGQ() : Math.min(minGQ, genotype.getGQ());
        }

        end = newEnd;
        if (genotype.hasDP()) {
            DPs.add(Math.max(genotype.getDP(), 0)); // DP must be >= 0
        }
    }

    /**
     * Lowers each value of {@code mins} to the corresponding value of {@code values}.
     *
     * @param what the name of the values, for the error message
     * @throws GATKException if the arrays differ in length
     */
    private static void lowerToMin(final int[] mins, final int[] values, final String what) {
        if (values.length != mins.length) {
            throw new GATKException("trying to merge different " + what + " array sizes: " + values.length + " != " + mins.length);
        }
        for (int i = 0; i < values.length; i++) {
            mins[i] = Math.min(mins[i], values[i]);
        }
    }

    /** Get the min PLs observed within this band, can be null if no PLs have yet been observed */
    public int[] getMinPLs() {
        return minPLs;
    }

    /** Get the min PPs observed within this band, can be null if no PPs have yet been observed */
    public int[] getMinPPs() {
        return minPPs;
    }

    /**
    * @return the ploidy of this hom-ref block.
     */
    public int getPloidy() {
        return ploidy;
    }

    public int getMinGQ() {
        return minGQ;
    }
}
