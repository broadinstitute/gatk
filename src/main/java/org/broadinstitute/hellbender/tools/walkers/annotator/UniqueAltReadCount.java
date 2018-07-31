package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
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
public class UniqueAltReadCount extends GenotypeAnnotation {
    public static final String UNIQUE_ALT_READ_SET_COUNT_KEY = "UNIQ_ALT_READ_COUNT";

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(UNIQUE_ALT_READ_SET_COUNT_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(UNIQUE_ALT_READ_SET_COUNT_KEY, 1, VCFHeaderLineType.Integer,
                "Number of ALT reads with unique start and mate end positions at a variant site"));
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        if (g.isHomRef()) {
            // skip the normal sample
            return;
        }

        final Allele altAllele = vc.getAlternateAllele(0); // assume single-allelic
        final String tumorSampleName = g.getSampleName();
        Collection<ReadLikelihoods<Allele>.BestAllele> tumorBestAlleles = likelihoods.bestAllelesBreakingTies(tumorSampleName);

        // Build a map from the (Start Position, Fragment Size) tuple to the count of reads with that
        // start position and fragment size
        Map<ImmutablePair<Integer, Integer>, Long> duplicateReadMap = tumorBestAlleles.stream()
                .filter(ba -> ba.allele.equals(altAllele) && ba.isInformative())
                .map(ba -> new ImmutablePair<>(ba.read.getStart(), ba.read.getFragmentLength()))
                .collect(Collectors.groupingBy(x -> x, Collectors.counting()));

        gb.attribute(UNIQUE_ALT_READ_SET_COUNT_KEY, duplicateReadMap.size());
    }
}
