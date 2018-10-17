package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Original Alignment annotation counts the number of alt reads where the original alignment contig doesn't match the current alignment contig
 *
 * <p>If reads were realigned to multiple references (for example the full human reference followed by just the
 * mitochondrial contig) and the original alignment tag was recorded before realigning, then we can count the number of alt
 * reads that have been realigned from other contigs to this one. This can be useful for downstream filtering if an alt
 * allele has all or most of its support originally from a different contig. In the mitochondrial case this can be useful
 * for filtering known NuMTs that are present in other contigs in the reference.  </p>
 *
 * <h3>Caveat</h3>
 * <p>This annotation can only be calculated if an OA tag is present in the bam and can only be run by Mutect2.</p>
 *
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of alt reads with an OA tag that doesn't match the current alignment contig.")
public class OriginalAlignment extends GenotypeAnnotation implements StandardMutectAnnotation {
    protected final OneShotLogger warning = new OneShotLogger(this.getClass());
    public static final String KEY = GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY;

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);

        final double[] lods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        if (lods==null) {
            warning.warn(String.format("One or more variant contexts is missing the 'TLOD' annotation, %s will not be computed for these VariantContexts", GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY));
            return;
        }
        final int indexOfMaxLod = MathUtils.maxElementIndex(lods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxLod);
        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        final String currentContig = ref.getInterval().getContig();

        final long nonChrMAlt = bestAlleles.stream()
                .filter(ba -> ba.read.hasAttribute(AddOriginalAlignmentTags.OA_TAG_NAME) && ba.isInformative() && ba.allele.equals(altAlelle) &&
                        !AddOriginalAlignmentTags.getOAContig(ba.read).equals(currentContig))
                .count();
        gb.attribute(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, nonChrMAlt);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(KEY));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(KEY);
    }
}
