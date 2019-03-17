package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.sam.markduplicates.MarkDuplicates;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Finds a lower bound on the number of unique reads at a locus that support a non-reference allele.
 *
 * <p>Multiple reads with the same start position and fragment length are grouped and counted only once as they are
 * likely duplicates.  In most cases such reads should be filtered using a tool such as {@link MarkDuplicates}.  This annotation
 * is designed for use with unique molecular identifiers (UMIs), in which case reads with the same start and fragment length but different
 * UMIs would appear to be independent.  This is not a default annotation of any GATK tool but can be enabled on the command line
 * with --annotation UniqueAltReadCount.</p>
 *
 * <p>Although these reads have different UMIs, sometimes they really are PCR duplicates.
 * We now believe that these duplicates are the result of a false-priming event that occurs during PCR amplification
 * in which excess adapter remains after the ligation step and fails to be completely
 * cleaned up during SPRI. This excess adapter is thought to act as a PCR primer during amplification, which leads to
 * the synthesis of a molecule with the wrong UMI.</p>
 *
 * <p>This annotation does not require or use any BAM file duplicate flags or UMI information, just the read alignments.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of non-duplicate-insert ALT reads (UNIQ_ALT_READ_COUNT)")
public class UniqueAltReadCount extends InfoFieldAnnotation {
    public static final String KEY = GATKVCFConstants.UNIQUE_ALT_READ_SET_COUNT_KEY;

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(KEY));
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        final Allele altAllele = vc.getAlternateAllele(0); // assume single-allelic

        // Build a map from the (Start Position, Fragment Size) tuple to the count of reads with that
        // start position and fragment size
        Map<ImmutablePair<Integer, Integer>, Long> duplicateReadMap = likelihoods.bestAllelesBreakingTies().stream()
                .filter(ba -> ba.allele.equals(altAllele) && ba.isInformative())
                .map(ba -> new ImmutablePair<>(ba.evidence.getStart(), ba.evidence.getFragmentLength()))
                .collect(Collectors.groupingBy(x -> x, Collectors.counting()));

        return ImmutableMap.of(KEY, duplicateReadMap.size());
    }
}
