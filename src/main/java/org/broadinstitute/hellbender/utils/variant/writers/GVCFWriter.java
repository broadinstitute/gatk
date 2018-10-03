package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.List;

import static htsjdk.variant.vcf.VCFConstants.MAX_GENOTYPE_QUAL;

/**
 * Genome-wide VCF writer
 */
public final class GVCFWriter implements VariantContextWriter {

    public final static String GVCF_BLOCK = "GVCFBlock";

    /** Where we'll ultimately write our VCF records */
    private final VariantContextWriter underlyingWriter;

    private final RangeMap<Integer, Range<Integer>> gqPartitions;
    private final int defaultPloidy;

    /** fields updated on the fly during GVCFWriter operation */
    private int nextAvailableStart = -1;
    private String contigOfNextAvailableStart = null;
    private String sampleName = null;
    private HomRefBlock currentBlock = null;

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param gqPartitions     a list of GQ partitions, this list must be non-empty and every element must be larger than previous element
     * @param defaultPloidy    the assumed ploidy for input variant context without one.
     */
    public GVCFWriter(final VariantContextWriter underlyingWriter, final List<Integer> gqPartitions, final int defaultPloidy) {
        this.underlyingWriter = Utils.nonNull(underlyingWriter);
        this.gqPartitions = parsePartitions(gqPartitions);
        this.defaultPloidy = defaultPloidy;
    }

    /**
     * Create {@link HomRefBlock}s which will collectively accept variants of any genotype quality
     *
     * Each individual block covers a band of genotype qualities with the splits between bands occurring at values in {@code gqPartitions}.
     * There will be {@code gqPartitions.size() +1} bands produced covering the entire possible range of genotype qualities from 0 to {@link VCFConstants#MAX_GENOTYPE_QUAL}.
     *
     * @param gqPartitions proposed GQ partitions
     * @return a list of HomRefBlocks accepting bands of genotypes qualities split at the points specified in gqPartitions
     */
    @VisibleForTesting
    static RangeMap<Integer,Range<Integer>> parsePartitions(final List<Integer> gqPartitions) {
        Utils.nonEmpty(gqPartitions);
        Utils.containsNoNull(gqPartitions, "The list of GQ partitions contains a null integer");
        final RangeMap<Integer, Range<Integer>> result = TreeRangeMap.create();
        int lastThreshold = 0;
        for (final Integer value : gqPartitions) {
            if (value < 0) {
                throw new IllegalArgumentException("The list of GQ partitions contains a non-positive integer.");
            } else if (value > MAX_GENOTYPE_QUAL + 1) {
                throw new IllegalArgumentException(String.format("The value %d in the list of GQ partitions is greater than VCFConstants.MAX_GENOTYPE_QUAL + 1 = %d.", value, MAX_GENOTYPE_QUAL + 1));
            } else if (value < lastThreshold) {
                throw new IllegalArgumentException(String.format("The list of GQ partitions is out of order. Previous value is %d but the next is %d.", lastThreshold, value));
            } else if (value == lastThreshold) {
                throw new IllegalArgumentException(String.format("The value %d appears more than once in the list of GQ partitions.", value));
            }

            result.put(Range.closedOpen(lastThreshold, value), Range.closedOpen(lastThreshold, value));
            lastThreshold = value;
        }

        if (lastThreshold <= MAX_GENOTYPE_QUAL) {
            result.put(Range.closedOpen(lastThreshold, MAX_GENOTYPE_QUAL + 1), Range.closedOpen(lastThreshold,MAX_GENOTYPE_QUAL + 1));
        }

        return result;
    }

    /**
     * Write the VCF header
     *
     * Adds standard GVCF fields to the header
     *
     * @param header a non-null header
     */
    @Override
    public void writeHeader(VCFHeader header) {
        Utils.nonNull(header, "header cannot be null");

        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MIN_DP_FORMAT_KEY));

        for (final Range<Integer> partition : gqPartitions.asMapOfRanges().keySet()) {
            header.addMetaDataLine(rangeToVCFHeaderLine(partition));
        }

        underlyingWriter.writeHeader(header);
    }

    static VCFHeaderLine rangeToVCFHeaderLine(Range<Integer> genotypeQualityBand) {
        // Need to uniquify the key for the header line using the min/max GQ, since
        // VCFHeader does not allow lines with duplicate keys.
        final String key = String.format(GVCF_BLOCK+"%d-%d", genotypeQualityBand.lowerEndpoint(), genotypeQualityBand.upperEndpoint());
        return new VCFHeaderLine(key, "minGQ=" + genotypeQualityBand.lowerEndpoint() + "(inclusive),maxGQ=" + genotypeQualityBand.upperEndpoint() + "(exclusive)");
    }


    /**
     * Close this GVCF writer.  Finalizes any pending hom-ref blocks and emits those to the underlyingWriter as well
     */
    @Override
    public void close() {
        try {
            emitCurrentBlock();
        } finally {
            underlyingWriter.close();
        }
    }

    @Override
    public boolean checkError() {
        return underlyingWriter.checkError();
    }

    /**
     * Add hom-ref site from vc to this gVCF hom-ref state tracking, emitting any pending states if appropriate
     *
     * @param vc a non-null VariantContext
     * @param g  a non-null genotype from VariantContext
     * @return a VariantContext to be emitted, or null if non is appropriate
     */
    protected VariantContext addHomRefSite(final VariantContext vc, final Genotype g) {

        if (nextAvailableStart != -1) {
            //there's a use case here related to ReblockGVCFs for overlapping deletions on different haplotypes
            if ( vc.getStart() <= nextAvailableStart && vc.getContig().equals(contigOfNextAvailableStart) ) {
                if (vc.getEnd() <= nextAvailableStart)
                    return null;
            }
            // otherwise, reset to non-relevant
            nextAvailableStart = -1;
            contigOfNextAvailableStart = null;
        }

        final VariantContext result;
        if (genotypeCanBeMergedInCurrentBlock(g)) {
            currentBlock.add(vc.getStart(), vc.getAttributeAsInt(VCFConstants.END_KEY, vc.getStart()), g);
            result = null;
        } else {
            result = currentBlock != null ? currentBlock.toVariantContext(sampleName): null;
            currentBlock = createNewBlock(vc, g);
        }
        return result;
    }

    private boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        return currentBlock != null
                && currentBlock.withinBounds(Math.min(g.getGQ(), MAX_GENOTYPE_QUAL))
                && currentBlock.getPloidy() == g.getPloidy()
                && (currentBlock.getMinPLs() == null || !g.hasPL() || (currentBlock.getMinPLs().length == g.getPL().length));
    }

    /**
     * Flush the current hom-ref block, if necessary, to the underlying writer, and reset the currentBlock to null
     */
    private void emitCurrentBlock() {
        if (currentBlock != null) {
            underlyingWriter.add(currentBlock.toVariantContext(sampleName));
            currentBlock = null;
        }
    }


    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g  the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    private HomRefBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the GQ of g
        final int gq = Math.min(g.getGQ(), MAX_GENOTYPE_QUAL);
        final Range<Integer> partition = gqPartitions.get(gq);

        if( partition == null) {
            throw new GATKException("GQ " + g + " from " + vc + " didn't fit into any partition");
        }

        // create the block, add g to it, and return it for use
        final HomRefBlock block = new HomRefBlock(vc, partition.lowerEndpoint(), partition.upperEndpoint(), defaultPloidy);
        block.add(vc.getStart(), vc.getAttributeAsInt(VCFConstants.END_KEY, vc.getStart()), g);
        return block;
    }

    /**
     * Add a VariantContext to this writer for emission
     *
     * Requires that the VC have exactly one genotype
     *
     * @param vc a non-null VariantContext
     */
    @Override
    public void add(VariantContext vc) {
        Utils.nonNull(vc);
        Utils.validateArg(vc.hasGenotypes(), "GVCF assumes that the VariantContext has genotypes");
        Utils.validateArg(vc.getGenotypes().size() == 1, () -> "GVCF assumes that the VariantContext has exactly one genotype but saw " + vc.getGenotypes().size());

        if (sampleName == null) {
            sampleName = vc.getGenotype(0).getSampleName();
        }

        if (currentBlock != null && !currentBlock.isContiguous(vc)) {
            // we've made a non-contiguous step (across interval, onto another chr), so finalize
            emitCurrentBlock();
        }

        final Genotype g = vc.getGenotype(0);
        if (g.isHomRef() && vc.hasAlternateAllele(Allele.NON_REF_ALLELE) && vc.isBiallelic()) {
            // create bands
            final VariantContext maybeCompletedBand = addHomRefSite(vc, g);
            if (maybeCompletedBand != null) {
                underlyingWriter.add(maybeCompletedBand);
            }
        } else {
            // g is variant, so flush the bands and emit vc
            emitCurrentBlock();
            nextAvailableStart = vc.getEnd();
            contigOfNextAvailableStart = vc.getContig();
            underlyingWriter.add(vc);
        }

    }

    @Override
    public void setHeader(VCFHeader header) {
        underlyingWriter.setHeader(header);
    }
}
