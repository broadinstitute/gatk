package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
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
final class HomRefBlock implements Locatable {

    private static final int HOM_REF_PL_POSITION = 0;  //the first value in the minPL[] is always the HomRef

    private final VariantContext startingVC;
    private final int minGQ, maxGQ;
    private final List<Integer> DPs = new ArrayList<>();
    private final Allele ref;
    private final int ploidy;

    private int end;
    private int[] minPLs = null;
    private int[] minPPs = null;

    /**
     * Create a new HomRefBlock
     *
     * @param startingVC the VariantContext that starts this band (for starting position information)
     * @param lowerGQBound the lowerGQBound (inclusive) to use in this band
     * @param upperGQBound the upperGQBound (exclusive) to use in this band
     */
    public HomRefBlock(final VariantContext startingVC, final int lowerGQBound, final int upperGQBound, final int defaultPloidy) {
        Utils.nonNull(startingVC, "startingVC cannot be null");
        Utils.validateArg(upperGQBound <= VCFConstants.MAX_GENOTYPE_QUAL + 1, "upperGQBound must be <= " + (VCFConstants.MAX_GENOTYPE_QUAL + 1));
        if ( lowerGQBound > upperGQBound ) { throw new IllegalArgumentException("bad lowerGQBound " + lowerGQBound + " as it's >= upperGQBound " + upperGQBound); }

        this.startingVC = startingVC;
        this.end = getStart() - 1;
        this.ref = startingVC.getReference();
        this.minGQ = lowerGQBound;
        this.maxGQ = upperGQBound;
        this.ploidy = startingVC.getMaxPloidy(defaultPloidy);
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
        vcb.attributes(new LinkedHashMap<>(2)); // clear the attributes
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
        final int[] minPPs = getMinPPs();
        gb.PL(minPLs);
        gb.GQ(GATKVariantContextUtils.calculateGQFromPLs(minPPs != null? minPPs : minPLs));
        gb.DP(getMedianDP());
        gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());
        if (minPPs != null) {
            gb.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY, Utils.listFromPrimitives(minPPs));
        }

        return gb.make();
    }

    /**
     * Add information from this Genotype to this band.
     *
     * Treats GQ values > 99 as 99.
     *
     * @param pos Current genomic position. Must be 1 base after the previous position
     * @param genotype A non-null Genotype with GQ and DP attributes
     */
    public void add(final int pos, final Genotype genotype) {
        add(pos, pos, genotype);
    }

    /**
     * Add a homRef block to the current block
     *
     * @param pos current genomic position
     * @param newEnd new calculated block end position
     * @param genotype A non-null Genotype with GQ and DP attributes
     */
    public void add(final int pos, final int newEnd, final Genotype genotype) {
        Utils.nonNull(genotype, "genotype cannot be null");
        if ( ! genotype.hasPL() ) { throw new IllegalArgumentException("genotype must have PL field");}
        if ( pos != end + 1 ) { throw new IllegalArgumentException("adding genotype at pos " + pos + " isn't contiguous with previous end " + end); }
        if ( genotype.getPloidy() != ploidy) { throw new IllegalArgumentException("cannot add a genotype with a different ploidy: " + genotype.getPloidy() + " != " + ploidy); }
        // Make sure the GQ is within the bounds of this band. Treat GQs > 99 as 99.
        if ( !withinBounds(Math.min(genotype.getGQ(), VCFConstants.MAX_GENOTYPE_QUAL))) {
            throw new IllegalArgumentException("cannot add a genotype with GQ=" + genotype.getGQ() + " because it's not within bounds ["
                    + this.getGQLowerBound() + ',' + this.getGQUpperBound() + ')');
        }

        if( minPLs == null ) {
            minPLs = genotype.getPL();
        }
        else { // otherwise take the min with the provided genotype's PLs
            final int[] pls = genotype.getPL();
            if (pls.length != minPLs.length) {
                throw new GATKException("trying to merge different PL array sizes: " + pls.length + " != " + minPLs.length);
            }
            for (int i = 0; i < pls.length; i++) {
                minPLs[i] = Math.min(minPLs[i], pls[i]);
            }
        }

        if( genotype.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY)) {
            if (minPPs == null ) {
                minPPs = PosteriorProbabilitiesUtils.parsePosteriorsIntoPhredSpace(genotype);
            }
            else { // otherwise take the min with the provided genotype's PLs
                final int[] pps = PosteriorProbabilitiesUtils.parsePosteriorsIntoPhredSpace(genotype);
                if (pps.length != minPPs.length) {
                    throw new GATKException("trying to merge different PP array sizes: " + pps.length + " != " + minPPs.length);
                }
                for (int i = 0; i < pps.length; i++) {
                    minPPs[i] = Math.min(minPPs[i], pps[i]);
                }
            }
        }

        end = newEnd;
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

    /** Get the min PPs observed within this band, can be null if no PPs have yet been observed */
    public int[] getMinPPs() {
        return minPPs;
    }

    int getGQUpperBound() {
        return maxGQ;
    }
    int getGQLowerBound() {
        return minGQ;
    }

    public boolean isContiguous(final VariantContext vc) {
        return (vc.getStart() == getEnd() + 1) && startingVC.getContig().equals(vc.getContig());
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

    public boolean canMergeWith(final Genotype g) {
        Utils.nonNull(g);
        return withinBounds(g.getGQ())
                && getPloidy() == g.getPloidy()
                && (getMinPLs() == null || !g.hasPL() || getMinPLs().length == g.getPL().length);
    }
}
