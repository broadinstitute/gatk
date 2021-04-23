package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class SVCallRecordWithEvidence extends SVCallRecord {

    private final List<SplitReadSite> startSplitReadSites;
    private final List<SplitReadSite> endSplitReadSites;
    private final List<DiscordantPairEvidence> discordantPairs;

    public SVCallRecordWithEvidence(final SVCallRecord record) {
        super(record.getContig(), record.getStart(), record.getStartStrand(), record.getEndContig(), record.getEnd(),
                record.getEndStrand(), record.getType(), record.getLength(), record.getAlgorithms(), record.getGenotypes());
        this.startSplitReadSites = Collections.emptyList();
        this.endSplitReadSites = Collections.emptyList();
        this.discordantPairs = Collections.emptyList();
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
                                    final List<Genotype> genotypes,
                                    final List<SplitReadSite> startSplitReadSites,
                                    final List<SplitReadSite> endSplitReadSites,
                                    final List<DiscordantPairEvidence> discordantPairs) {
        super(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, genotypes);
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

    @Override
    public boolean equals(final Object obj) {
        if (!super.equals(obj)) {
            return false;
        }
        if (this.getClass() != obj.getClass()) {
            return false;
        }
        final SVCallRecordWithEvidence b = (SVCallRecordWithEvidence)obj;
        boolean areEqual = this.getDiscordantPairs().containsAll(b.getDiscordantPairs());
        areEqual &= b.getDiscordantPairs().containsAll(this.getDiscordantPairs());
        areEqual &= this.getEndSplitReadSites().containsAll(b.getEndSplitReadSites());
        areEqual &= b.getEndSplitReadSites().containsAll(this.getEndSplitReadSites());
        areEqual &= this.getStartSplitReadSites().containsAll(b.getStartSplitReadSites());
        areEqual &= b.getStartSplitReadSites().containsAll(this.getStartSplitReadSites());
        return areEqual;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), discordantPairs, endSplitReadSites, startSplitReadSites);
    }

    String prettyPrint() {
        return getContig() + ":" + getStart() + "-" + getEnd();
    }
}
