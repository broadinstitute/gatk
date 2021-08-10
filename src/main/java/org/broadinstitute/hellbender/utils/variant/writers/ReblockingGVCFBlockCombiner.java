package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ReblockGVCF;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.PushPullTransformer;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Combines variants into GVCF blocks.
 */
public class ReblockingGVCFBlockCombiner extends GVCFBlockCombiner implements PushPullTransformer<VariantContext> {
    private static final Logger logger = LogManager.getLogger(ReblockingGVCFBlockCombiner.class);

    private final List<VariantContextBuilder> homRefBlockBuffer = new ArrayList<>(10);  //10 is a generous estimate for the number of overlapping deletions
    private static final Comparator<? super VariantContextBuilder> VCB_COMPARATOR = Comparator.comparingLong(VariantContextBuilder::getStart);

    final private CachingIndexedFastaSequenceFile referenceReader;

    /**
     * fields updated on the fly during GVCFWriter operation
     */
    private int vcfOutputEnd = 1;
    private int bufferEnd = 0;
    final private boolean dropLowQuals;
    final private boolean allowMissingHomRefData;
    final private double rgqThreshold;
    private String currentContig = null;

    ReblockingGVCFBlockCombiner(final List<? extends Number> gqPartitions, final boolean floorBlocks,
                                       final CachingIndexedFastaSequenceFile referenceReader, final ReblockingOptions options) {
        super(gqPartitions, floorBlocks);
        this.referenceReader = referenceReader;
        this.dropLowQuals = options.getDropLowQualsOpt();
        this.allowMissingHomRefData = options.getAllowMissingHomRefDataOpt();
        this.rgqThreshold = options.getRgqThresholdOpt();
    }

    /**
     * Add hom-ref site from vc to this gVCF hom-ref state tracking, emitting any pending states if appropriate
     *
     * @param vc a non-null VariantContext
     * @param g  a non-null genotype from VariantContext
     * @return a VariantContext to be emitted, or null if non is appropriate
     */
    protected VariantContext addHomRefSite(final VariantContext vc, final Genotype g) {
        final Genotype genotype = vc.getGenotype(0);
        final VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);

        if (dropLowQuals && (!genotype.hasGQ() || genotype.getGQ() < rgqThreshold || genotype.getGQ() == 0)) {
            return null;
        } else if (genotype.isHomRef()) {
            if (!genotype.hasPL()) {
                if (genotype.hasGQ()) {
                    logger.warn("PL is missing for hom ref genotype at at least one position for sample " + genotype.getSampleName() + ": " + vc.getContig() + ":" + vc.getStart() +
                            ".  Using GQ to determine quality.");
                    final int gq = genotype.getGQ();
                    final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                    vcBuilder.genotypes(gBuilder.GQ(gq).make());
                    return super.addHomRefSite(vcBuilder.make(), gBuilder.make());
                } else {
                    final String message = "Homozygous reference genotypes must contain GQ or PL. Both are missing for hom ref genotype at "
                            + vc.getContig() + ":" + vc.getStart() + " for sample " + genotype.getSampleName() + ".";
                    if (allowMissingHomRefData) {
                        logger.warn(message);
                        final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                        vcBuilder.genotypes(gBuilder.GQ(0).PL(new int[]{0,0,0}).make());
                        return super.addHomRefSite(vcBuilder.make(), gBuilder.make());
                    } else {
                        throw new UserException.BadInput(message);
                    }
                }
            }
            return super.addHomRefSite(vcBuilder.make(), genotype);
            //some external data has no-called genotypes with good likelihoods
        } else if (!genotype.isCalled() && genotype.hasPL() && genotype.getPL()[0] == 0) {
            return super.addHomRefSite(vcBuilder.make(), genotype);
        }
        else {
            return null;
        }
    }

    /**
     * Add the input as a new reference block or write and remove ref blocks that end before the variantContext if it is variant
     * Trim ref block if the variant occurs in the middle of a block
     * @param variantContextToOutput may be variant or reference, can overlap existing ref blocks in buffer,
     *                              but should never start before vcfOutputEnd
     */
    @Override
    public void submit(final VariantContext variantContextToOutput) {
        if (variantContextToOutput == null) {
            return;
        }
        Utils.validate(variantContextToOutput.getStart() <= variantContextToOutput.getEnd(),
                () -> "Input variant context at position " + currentContig + ":" + variantContextToOutput.getStart() + " has negative length: start=" + variantContextToOutput.getStart() + " end=" + variantContextToOutput.getEnd());

        if (currentContig == null) {
            currentContig = variantContextToOutput.getContig();
        } else if (!variantContextToOutput.getContig().equals(currentContig)) {
            flushRefBlockBuffer();
            currentContig = variantContextToOutput.getContig();
            vcfOutputEnd = 1;  //must be one to be a valid SimpleInterval
        }
        final VariantContextBuilder newHomRefBlock = new VariantContextBuilder(variantContextToOutput);

        final Genotype g = variantContextToOutput.getGenotype(0);
        if (isHomRef(g)) {
            if (variantContextToOutput.getStart() <= vcfOutputEnd) {
                if (variantContextToOutput.getEnd() <= vcfOutputEnd) {
                    return;
                }
                moveBuilderStart(newHomRefBlock, vcfOutputEnd + 1, referenceReader);
            }
        }

        final List<VariantContextBuilder> completedBlocks = new ArrayList<>();
        final List<VariantContextBuilder> tailBuffer = new ArrayList<>();
        for (final VariantContextBuilder builder : homRefBlockBuffer) {
            final int blockStart = (int)builder.getStart();
            final int variantEnd = variantContextToOutput.getEnd();
            if (blockStart > variantEnd) {
                if ((!g.isHomRef() || (g.hasPL() && g.getPL()[0] != 0))) {
                    super.submit(variantContextToOutput);
                    vcfOutputEnd = Math.max(vcfOutputEnd, variantEnd);
                } else {
                    homRefBlockBuffer.add(newHomRefBlock);
                }
                homRefBlockBuffer.removeAll(completedBlocks);
                homRefBlockBuffer.addAll(tailBuffer);
                homRefBlockBuffer.sort(VCB_COMPARATOR);
                return;
            }
            int blockEnd = (int)builder.getStop();
            final int variantStart = variantContextToOutput.getStart();
            if (blockEnd >= variantStart) {  //then trim out overlap
                blockEnd = trimBlockToVariant(variantContextToOutput, completedBlocks, tailBuffer, builder);
            }
            //only flush ref blocks if we're outputting a variant, otherwise ref blocks can be out of order
            if (blockStart < variantStart && !g.isHomRef()) {
                super.submit(builder.make());
                vcfOutputEnd = blockEnd;
                completedBlocks.add(builder);
            }
        }
        homRefBlockBuffer.removeAll(completedBlocks);
        if (isHomRef(g)) {
            if (ReblockGVCF.isHomRefBlock(variantContextToOutput) && newHomRefBlock.getStart() < vcfOutputEnd) {
                throw new IllegalStateException("Reference positions added to buffer should not overlap positions already output to VCF. "
                        + variantContextToOutput.getStart() + " overlaps position " + currentContig + ":" + vcfOutputEnd + " already emitted.");
            }
            homRefBlockBuffer.add(newHomRefBlock);
            bufferEnd = Math.max(bufferEnd, (int)newHomRefBlock.getStop());
        } else if (bufferEnd >= variantContextToOutput.getEnd() || homRefBlockBuffer.isEmpty()) {  //
            super.submit(variantContextToOutput);
            vcfOutputEnd = Math.max(vcfOutputEnd, variantContextToOutput.getEnd());  //there may have been a previous, larger deletion
        }
        homRefBlockBuffer.addAll(tailBuffer);
        homRefBlockBuffer.sort(VCB_COMPARATOR);  //this may seem lazy, but it's more robust to assumptions about overlap being violated
    }

    //funky logic for DRAGEN GVCFs where call may not match PLs
    private boolean isHomRef(final Genotype g) {
        return (g.isHomRef() && !g.hasPL()) || (g.hasPL() && g.getPL()[0] == 0);
    }

    /**
     * Remove overlap between a variant and a ref block by trimming the ref block
     * @param variantContextToOutput    a variant
     * @param completedBlocks   the list of blocks that are finalized and will be removed from the buffer
     * @param tailBuffer    the list of ref blocks that will stay in the hom ref block buffer without getting output
     * @param builder   the current ref block, overlapping the variant
     * @return  the new end of the current ref block
     */
    private int trimBlockToVariant(final VariantContext variantContextToOutput, final List<VariantContextBuilder> completedBlocks,
                                   final List<VariantContextBuilder> tailBuffer, final VariantContextBuilder builder) {
        final int blockStart = (int)builder.getStart();
        final int variantEnd = variantContextToOutput.getEnd();
        int blockEnd = (int)builder.getStop();
        final int variantStart = variantContextToOutput.getStart();
        if (blockEnd > variantEnd && blockStart < variantStart) {  //then this block will be split -- add a post-variant block
            final VariantContextBuilder blockTailBuilder = new VariantContextBuilder(builder);
            moveBuilderStart(blockTailBuilder, variantEnd + 1, referenceReader);
            tailBuffer.add(blockTailBuilder);
            builder.stop(variantStart-1);
            builder.attribute(VCFConstants.END_KEY, variantStart-1);
            blockEnd = variantStart-1;
        }
        if (blockStart < variantStart) { //right trim
            builder.attribute(VCFConstants.END_KEY, variantStart - 1);
            builder.stop(variantStart - 1);
        } else {  //left trim
            if (variantContextToOutput.contains(new SimpleInterval(currentContig, blockStart, blockEnd))) {
                //if block is entirely overlapped by VC to output, then remove it from the buffer, but don't output
                completedBlocks.add(builder);
            } else {
                if (blockEnd < variantEnd + 1) {
                    throw new GATKException.ShouldNeverReachHereException("ref block end overlaps variant end; current builder: " + builder.getStart() + " to " + builder.getStop());
                }
                moveBuilderStart(builder, variantEnd + 1, referenceReader);
            }
        }
        if (builder.getStart() > builder.getStop()) {
            throw new GATKException.ShouldNeverReachHereException("builder start follows stop; current builder: " + builder.getStart() + " to " + builder.getStop());
        }
        return blockEnd;
    }

    /**
     * Write all the reference blocks in the buffer to the output VCF
     */
    private void flushRefBlockBuffer() {
        for (final VariantContextBuilder builder : homRefBlockBuffer) {
            //toOutput.add(builder.make());
            try {
                super.submit(builder.make());
            } catch (Exception e) {
                throw new GATKException("builder threw an exception at " + builder.getContig() + ":" + builder.getStart(), e);
            }
            vcfOutputEnd = (int)builder.getStop();
        }
        homRefBlockBuffer.clear();
        bufferEnd = 0;
    }

    /**
     * Modifies ref block builder to change start position and update ref allele accordingly in VC and genotypes
     * @param builder   a builder for a reference block, contains only NON_REF, no other ALTs
     * @param newStart  the new position for the reference block
     */
    public static void moveBuilderStart(final VariantContextBuilder builder, final int newStart, final CachingIndexedFastaSequenceFile referenceReader) {
        final byte[] newRef = ReferenceUtils.getRefBaseAtPosition(referenceReader, builder.getContig(), newStart);
        final Allele newRefAllele = Allele.create(newRef, true);
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        for (final Genotype g : builder.getGenotypes()) {
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            final List<Allele> newGTAlleles = g.getAlleles().stream().map(a -> a.isReference() ? newRefAllele : a).collect(Collectors.toList());
            gb.alleles(newGTAlleles);
            genotypesArray.add(gb.make());
        }
        final List<Allele> newVCAlleles = builder.getAlleles().stream().map(a -> a.isReference() ? newRefAllele : a).collect(Collectors.toList());
        builder.start(newStart).alleles(newVCAlleles).genotypes(genotypesArray);
    }

    public int getVcfOutputEnd() {
        return vcfOutputEnd;
    }

    public String getCurrentContig() { return currentContig; }

    public int getBufferEnd() { return bufferEnd; }

    public int getBufferStart() {
        return (int)homRefBlockBuffer.get(0).getStart();
    }

    @Override
    public void signalEndOfInput() {
        flushRefBlockBuffer();
        super.signalEndOfInput();
    }
}