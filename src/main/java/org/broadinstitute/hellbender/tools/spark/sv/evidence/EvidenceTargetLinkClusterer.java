package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;

/**
 * This class is responsible for iterating over a collection of BreakpointEvidence to find clusters of evidence with
 * distal targets (discordant read pairs or split reads) that agree in their location and target intervals and strands.
 * Paired intervals that agree are intersected with one another to produce what should be the smallest possible pair
 * of intervals that might hold the breakpoints. The evidence target links also keep track of the number of discordant
 * read pairs and split reads that went into their creation.
 */
public class EvidenceTargetLinkClusterer {

    private final ReadMetadata readMetadata;
    private final int minEvidenceMapq;

    public EvidenceTargetLinkClusterer(final ReadMetadata readMetadata, final int minEvidenceMapq) {
        this.readMetadata = readMetadata;
        this.minEvidenceMapq = minEvidenceMapq;
    }

    public Iterator<EvidenceTargetLink> cluster(final Iterator<BreakpointEvidence> breakpointEvidenceIterator) throws Exception {
        final List<EvidenceTargetLink> links = new ArrayList<>();
        final PairedStrandedIntervalTree<EvidenceTargetLink> currentLinks = new PairedStrandedIntervalTree<>();
        while (breakpointEvidenceIterator.hasNext()) {
            final BreakpointEvidence nextEvidence = breakpointEvidenceIterator.next();
            if (nextEvidence.hasDistalTargets(readMetadata, minEvidenceMapq)) {
                Utils.validate(nextEvidence instanceof BreakpointEvidence.SplitRead || nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence,
                        "Unknown evidence type with distal target: " + nextEvidence);
                EvidenceTargetLink updatedLink = null;
                final String templateName = ((BreakpointEvidence.ReadEvidence) nextEvidence).getTemplateName();
                // todo: what to do if there are more than one distal targets -- for now just taking the first one
                // todo: this would only be an issue with split reads with more than one SA mapping, and we aren't
                // todo: computing distal targets for those yet

                EvidenceTargetLink nextEvidenceLink = new EvidenceTargetLink(
                        new StrandedInterval(nextEvidence.getLocation(),
                            nextEvidence.isEvidenceUpstreamOfBreakpoint()),
                        nextEvidence.getDistalTargets(readMetadata, minEvidenceMapq).get(0),
                        nextEvidence instanceof BreakpointEvidence.SplitRead
                                ? 1 : 0,
                        nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                ? 1 : 0,
                        nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence ? Collections.singleton(templateName) : new HashSet<>(),
                        nextEvidence instanceof BreakpointEvidence.SplitRead ? Collections.singleton(templateName) : new HashSet<>()
                );

                // if the new link overlaps with an existing link, combine the two links together, adding their evidence counts
                // and intersecting their intervals
                final Iterator<Tuple2<PairedStrandedIntervals, EvidenceTargetLink>> it = currentLinks.overlappers(nextEvidenceLink.getPairedStrandedIntervals());
                if (it.hasNext()) {
                    final Tuple2<PairedStrandedIntervals, EvidenceTargetLink> matchingETL = it.next();
                    final EvidenceTargetLink oldLink = matchingETL._2();
                    if (nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence && oldLink.getReadPairTemplateNames().contains(templateName)) continue;
                    if (nextEvidence instanceof BreakpointEvidence.SplitRead && oldLink.getSplitReadTemplateNames().contains(templateName)) continue;
                    it.remove();
                    final SVInterval newSource = oldLink.source.getInterval().intersect(nextEvidenceLink.source.getInterval());
                    final SVInterval newTarget = oldLink.target.getInterval().intersect(nextEvidenceLink.target.getInterval());
                    int newSplitReadCount = nextEvidence instanceof BreakpointEvidence.SplitRead
                            ? oldLink.splitReads + 1 : oldLink.splitReads;
                    int newReadPairCount = nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                            ? oldLink.readPairs + 1 : oldLink.readPairs;
                    final Set<String> readPairTemplateNames = new HashSet<>(oldLink.getReadPairTemplateNames());
                    if (nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence) readPairTemplateNames.add(templateName);
                    final Set<String> splitReadTemplateNames = new HashSet<>(oldLink.getSplitReadTemplateNames());
                    if (nextEvidence instanceof BreakpointEvidence.SplitRead) splitReadTemplateNames.add(templateName);
                    updatedLink = new EvidenceTargetLink(
                            new StrandedInterval(newSource, oldLink.source.getStrand()),
                            new StrandedInterval(newTarget, oldLink.target.getStrand()),
                            newSplitReadCount,
                            newReadPairCount,
                            readPairTemplateNames,
                            splitReadTemplateNames);
                }
                if (updatedLink == null) {
                    updatedLink = nextEvidenceLink;
                }
                currentLinks.put(updatedLink.getPairedStrandedIntervals(), updatedLink);
            }
        }

        currentLinks.forEach(entry -> links.add(entry._2()));

        return links.iterator();
    }

}
