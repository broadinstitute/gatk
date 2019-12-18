package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;
import java.util.Set;

public class SVCallRecordWithEvidence extends SVCallRecord {

    private final List<SplitReadSite> startSplitReadSites;
    private final List<SplitReadSite> endSplitReadSites;
    private final List<DiscordantPairEvidence> discordantPairs;

    public SVCallRecordWithEvidence(final SVCallRecord record) {
        super(record.getContig(), record.getStart(), record.getStartStrand(), record.getEndContig(), record.getEnd(),
                record.getEndStrand(), record.getType(), record.getLength(), record.getAlgorithms(), record.getSamples());
        this.startSplitReadSites = Collections.emptyList();
        this.endSplitReadSites = Collections.emptyList();
        this.discordantPairs = Collections.emptyList();
    }

    public SVCallRecordWithEvidence(final SVCallRecord record,
                                    List<SplitReadSite> startSplitReadSites,
                                    List<SplitReadSite> endSplitReadSites,
                                    List<DiscordantPairEvidence> discordantPairs) {
        super(record.getContig(), record.getStart(), record.getStartStrand(), record.getEndContig(), record.getEnd(),
                record.getEndStrand(), record.getType(), record.getLength(), record.getAlgorithms(), record.getSamples());
        Utils.nonNull(startSplitReadSites);
        Utils.nonNull(endSplitReadSites);
        Utils.nonNull(discordantPairs);
        Utils.containsNoNull(startSplitReadSites, "Encountered null in start split reads");
        Utils.containsNoNull(endSplitReadSites, "Encountered null in end split reads");
        Utils.containsNoNull(discordantPairs, "Encountered null in discordant pairs");
        this.startSplitReadSites = startSplitReadSites;
        this.endSplitReadSites = endSplitReadSites;
        this.discordantPairs = discordantPairs;
    }

    public SVCallRecordWithEvidence(final String startContig,
                                    final int start,
                                    final boolean startStrand,
                                    final String endContig,
                                    final int end,
                                    final boolean endStrand,
                                    final StructuralVariantType type,
                                    final int length,
                                    final List<String> algorithms,
                                    final Set<String> samples,
                                    final List<SplitReadSite> startSplitReadSites,
                                    final List<SplitReadSite> endSplitReadSites,
                                    final List<DiscordantPairEvidence> discordantPairs) {
        super(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, samples);
        Utils.nonNull(startSplitReadSites);
        Utils.nonNull(endSplitReadSites);
        Utils.nonNull(discordantPairs);
        Utils.containsNoNull(startSplitReadSites, "Encountered null in start split reads");
        Utils.containsNoNull(endSplitReadSites, "Encountered null in end split reads");
        Utils.containsNoNull(discordantPairs, "Encountered null in discordant pairs");
        this.startSplitReadSites = startSplitReadSites;
        this.endSplitReadSites = endSplitReadSites;
        this.discordantPairs = discordantPairs;
    }

    public List<DiscordantPairEvidence> getDiscordantPairs() {
        return discordantPairs;
    }

    public List<SplitReadSite> getStartSplitReadSites() {
        return startSplitReadSites;
    }

    public List<SplitReadSite> getEndSplitReadSites() {
        return endSplitReadSites;
    }
}
