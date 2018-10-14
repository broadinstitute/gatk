package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
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
                                                                              final ReferenceMultiSparkSource reference,
                                                                              final int impreciseVariantEvidenceThreshold,
                                                                              final int maxCallableImpreciseVariantDeletionSize,
                                                                              final Logger localLogger) {
        final List<VariantContext> impreciseVariants =
                Utils.stream(evidenceTargetLinks)
                        .map(p -> p._2)
                        .filter(EvidenceTargetLink::isImpreciseDeletion)
                        // todo: for now, we set a limit on the size of deletion we can call solely based on imprecise evidence
                        .filter(e -> e.getDistance() < maxCallableImpreciseVariantDeletionSize)
                        .filter(e -> e.getReadPairs() + e.getSplitReads() > impreciseVariantEvidenceThreshold)
                        .map(e -> createImpreciseDeletionVariant(e, metadata, reference))
                        .collect(Collectors.toList());

        localLogger.info("Called " + impreciseVariants.size() + " imprecise deletion variants");
        return impreciseVariants;
    }

    private static VariantContext createImpreciseDeletionVariant(final EvidenceTargetLink e,
                                                                 final ReadMetadata metadata,
                                                                 final ReferenceMultiSparkSource reference) {
        final int svLength = e.getPairedStrandedIntervals().getLeft().getInterval().midpoint() -
                             e.getPairedStrandedIntervals().getRight().getInterval().midpoint();
        final SvType svType = new SimpleSVType.ImpreciseDeletion(e, svLength, metadata, reference);
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromEvidenceTargetLink(e, svType);
    }
}
