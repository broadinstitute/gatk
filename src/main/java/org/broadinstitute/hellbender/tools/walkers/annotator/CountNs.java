package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import java.util.*;


/**
 * Apply a read-based annotation that reports the number of Ns seen at a given site. This is intended for use on consensus called data.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of Ns at the pileup")
public class CountNs extends GenotypeAnnotation {
    /**
     * Calculate annotations for each allele based on given VariantContext and likelihoods for a given genotype's sample
     * and add the annotations to the GenotypeBuilder.  By default annotations are only calculated for alt alleles but
     * implementations may override the {@code includeRefAllele()} method.  See parent class docs in {@link GenotypeAnnotation}.
     */

    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        if ( g == null || likelihoods == null ) {
            return;
        }
        long Count = likelihoods.sampleReads(likelihoods.indexOfSample(g.getSampleName())).stream().filter(read -> doesReadHaveN(read, vc)).count();

        gb.attribute(GATKVCFConstants.N_COUNT_KEY, Count);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(GATKVCFConstants.N_COUNT_KEY, 1, VCFHeaderLineType.Integer, "Counts Ns at site"));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.N_COUNT_KEY); }

    private Boolean doesReadHaveN(final GATKRead read, final VariantContext vc) {
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        return ( offset != ReadUtils.CLIPPING_GOAL_NOT_REACHED && !AlignmentUtils.isInsideDeletion(read.getCigar(), offset) && read.getBase(offset) == 'N');
    }
}
