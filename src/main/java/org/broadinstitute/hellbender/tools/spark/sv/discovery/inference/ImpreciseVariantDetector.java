package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

public class ImpreciseVariantDetector {


    public static List<VariantContext> callImpreciseDeletionFromEvidenceLinks(final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                              final ReadMetadata metadata,
                                                                              final ReferenceMultiSource reference,
                                                                              final int impreciseEvidenceVariantCallingThreshold,
                                                                              final Logger localLogger) {
        final List<VariantContext> impreciseVariants =
                Utils.stream(evidenceTargetLinks)
                        .map(p -> p._2)
                        .filter(EvidenceTargetLink::isImpreciseDeletion)
                        .filter(e -> e.getReadPairs() + e.getSplitReads() > impreciseEvidenceVariantCallingThreshold)
                        .map(e -> createImpreciseDeletionVariant(e, metadata, reference))
                        .collect(Collectors.toList());

        localLogger.info("Called " + impreciseVariants.size() + " imprecise deletion variants");
        return impreciseVariants;
    }

    private static VariantContext createImpreciseDeletionVariant(final EvidenceTargetLink e,
                                                                 final ReadMetadata metadata,
                                                                 final ReferenceMultiSource reference) {
        final SvType svType = new SimpleSVType.ImpreciseDeletion(e, metadata);
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromEvidenceTargetLink(e, svType, metadata, reference);
    }
}
