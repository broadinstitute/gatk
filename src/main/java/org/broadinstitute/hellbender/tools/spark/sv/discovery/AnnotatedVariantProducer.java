package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Given identified pair of breakpoints for a simple SV and its supportive evidence, i.e. chimeric alignments,
 * produce an annotated {@link VariantContext}.
 */
public class AnnotatedVariantProducer implements Serializable {
    private static final long serialVersionUID = 1L;


    public static VariantContext produceAnnotatedVcFromEvidenceTargetLink(final EvidenceTargetLink evidenceTargetLink,
                                                                          final SvType svType) {
        final PairedStrandedIntervals pairedStrandedIntervals = evidenceTargetLink.getPairedStrandedIntervals();
        final StrandedInterval strandedIntervalLeft = pairedStrandedIntervals.getLeft();
        final StrandedInterval strandedIntervalRight = pairedStrandedIntervals.getRight();
        final int start = strandedIntervalLeft.getInterval().midpoint();
        final int end = strandedIntervalRight.getInterval().midpoint();
        final VariantContextBuilder builder = svType
                .getBasicInformation()
                .attribute(GATKSVVCFConstants.CIPOS, produceCIInterval(start, strandedIntervalLeft.getInterval()))
                .attribute(GATKSVVCFConstants.CIEND, produceCIInterval(end, strandedIntervalRight.getInterval()))
                .attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, evidenceTargetLink.getReadPairs())
                .attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, evidenceTargetLink.getSplitReads());
        return builder.make();
    }

    public static List<VariantContext> annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(final List<VariantContext> assemblyDiscoveredVariants,
                                                                                              final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                                              final ReadMetadata metadata,
                                                                                              final SAMSequenceDictionary refDict,
                                                                                              final DiscoverVariantsFromContigAlignmentsArgumentCollection parameters,
                                                                                              final Logger localLogger) {

        final int originalEvidenceLinkSize = evidenceTargetLinks.size();
        final List<VariantContext> result = assemblyDiscoveredVariants
                .stream()
                .map(variant ->
                        annotateWithImpreciseEvidenceLinks(variant, evidenceTargetLinks, refDict, metadata,
                                                        parameters.assemblyImpreciseEvidenceOverlapUncertainty))
                .collect(Collectors.toList());
        localLogger.info("Used " + (originalEvidenceLinkSize - evidenceTargetLinks.size()) + " evidence target links to annotate assembled breakpoints");
        return result;
    }

    //==================================================================================================================

    private static VariantContext annotateWithImpreciseEvidenceLinks(final VariantContext variant,
                                                                     final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                     final SAMSequenceDictionary referenceSequenceDictionary,
                                                                     final ReadMetadata metadata,
                                                                     final int defaultUncertainty) {
        if (variant.getStructuralVariantType() == StructuralVariantType.DEL) {
            SVContext svc = SVContext.of(variant);
            final int padding = (metadata == null) ? defaultUncertainty : (metadata.getMaxMedianFragmentSize() / 2);
            PairedStrandedIntervals svcIntervals = svc.getPairedStrandedIntervals(metadata, referenceSequenceDictionary, padding);

            final Iterator<Tuple2<PairedStrandedIntervals, EvidenceTargetLink>> overlappers = evidenceTargetLinks.overlappers(svcIntervals);
            int readPairs = 0;
            int splitReads = 0;
            while (overlappers.hasNext()) {
                final Tuple2<PairedStrandedIntervals, EvidenceTargetLink> next = overlappers.next();
                readPairs += next._2.getReadPairs();
                splitReads += next._2.getSplitReads();
                overlappers.remove();
            }
            final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
            if (readPairs > 0) {
                variantContextBuilder.attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, readPairs);
            }
            if (splitReads > 0) {
                variantContextBuilder.attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, splitReads);
            }

            return variantContextBuilder.make();
        } else {
            return variant;
        }
    }

    /**
     * Produces the string representation of a VCF 4.2-style SV CI interval centered around 'point'.
     */
    @VisibleForTesting
    static String produceCIInterval(final int point, final SVInterval ciInterval) {
        Utils.validate(ciInterval.getStart() <= point && ciInterval.getEnd() >= point, "Interval must contain point");
        return String.join(",",
                String.valueOf(ciInterval.getStart() - point),
                String.valueOf(ciInterval.getEnd() - point));
    }

    /**
     * Apply filters (that implements {@link StructuralVariantFilter}) given list of variants,
     * and write the variants to a single VCF file.
     * @param variants variants to which filters are to be applied and written to file
     */
    public static List<VariantContext> filterMergedVariantList(
                                final List<VariantContext> variants,
                                final DiscoverVariantsFromContigAlignmentsArgumentCollection discoveryArgs ) {
        final List<VariantContext> variantsWithFilterApplied = new ArrayList<>(variants.size());
        final List<StructuralVariantFilter> filters = Arrays.asList(
                new SVMappingQualityFilter(discoveryArgs.minMQ),
                new SVAlignmentLengthFilter(discoveryArgs.minAlignLength));
        for (final VariantContext variant : variants) {
            String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_INS) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP)) {
                if (Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)) < StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND )
                    continue;
            }
            variantsWithFilterApplied.add(applyFilters(variant, filters));
        }
        return variantsWithFilterApplied;
    }

    /**
     * Filters out variants by testing against provided
     * filter key, threshold.
     *
     * Variants with value below specified threshold (or null value)
     * are filtered out citing given reason.
     *
     * @throws ClassCastException if the value corresponding to provided key cannot be casted as a {@link Double}
     */
    private static VariantContext applyFilters(final VariantContext variantContext,
                                               final List<StructuralVariantFilter> filters) {

        final Set<String> appliedFilters = new HashSet<>();
        for (final StructuralVariantFilter filter : filters) {
            if ( !filter.test(variantContext) )
                appliedFilters.add(filter.getName());
        }

        if (appliedFilters.isEmpty())
            return variantContext;
        else {
            return new VariantContextBuilder(variantContext).filters(appliedFilters).make();
        }
    }
}
