package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

public abstract class GVCFBlock implements Locatable {
    protected final VariantContext startingVC;
    protected final int minGQ;
    protected final int maxGQ;
    protected final Allele ref;
    protected final List<Integer> DPs = new ArrayList<>();
    protected int end;

    public GVCFBlock(final VariantContext startingVC, final int lowerGQBound, final int upperGQBound) {
        Utils.nonNull(startingVC, "startingVC cannot be null");
        this.startingVC = startingVC;
        this.minGQ = lowerGQBound;
        this.maxGQ = upperGQBound;
        this.ref = startingVC.getReference();
        this.end = startingVC.getAttributeAsInt(VCFConstants.END_KEY, getStart());

        final Genotype g = startingVC.getGenotype(0);
        if (g.hasDP()) {
            DPs.add(Math.max(g.getDP(), 0)); // DP must be >= 0
        }
    }

    public void add(int pos, Genotype genotype) {add(pos, pos, genotype);}

    /**
     * Convert a HomRefBlock into a VariantContext
     *
     * @param sampleName sample name to give this variant context
     * @return a VariantContext representing the gVCF encoding for this block.
     * It will return {@code null} if input {@code block} is {@code null}, indicating that there
     * is no variant-context to be output into the VCF.
     */
    public VariantContext toVariantContext(String sampleName, final boolean floorBlocks) {
        final VariantContextBuilder vcb = new VariantContextBuilder(getStartingVC());
        vcb.attributes(new LinkedHashMap<>(2)); // clear the attributes
        vcb.stop(getEnd());
        vcb.attribute(VCFConstants.END_KEY, getEnd());
        final Genotype genotype = createHomRefGenotype(sampleName, floorBlocks);

        return vcb.genotypes(genotype).make();
    }

    abstract Genotype createHomRefGenotype(String sampleName, boolean floorBlocks);

    public abstract void add(int pos, int newEnd, Genotype genotype);

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

    int getGQUpperBound() {
        return maxGQ;
    }

    int getGQLowerBound() {
        return minGQ;
    }

    /**
     * Allow overlapping input blocks, as in the case where shards are split with duplicate boundary events
     * @param vc
     * @return false if there is a gap between this block and input vc
     */
    public boolean isContiguous(final VariantContext vc) {
        return vc.withinDistanceOf(this, 1);
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
}
