package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public class SVCallRecordWithEvidence extends SVCallRecord {

    private final List<SplitReadSite> startSplitReadSites;
    private final List<SplitReadSite> endSplitReadSites;
    private final List<DiscordantPairEvidence> discordantPairs;
    private final CopyNumberPosteriorDistribution copyNumberDistribution;

    public SVCallRecordWithEvidence(final SVCallRecord record) {
        this(record, Collections.emptyList(), Collections.emptyList(), Collections.emptyList(), null);
    }

    public SVCallRecordWithEvidence(final SVCallRecord record,
                                    final List<SplitReadSite> startSplitReadSites,
                                    final List<SplitReadSite> endSplitReadSites,
                                    final List<DiscordantPairEvidence> discordantPairs,
                                    final CopyNumberPosteriorDistribution copyNumberDistribution) {
        this(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(),
                record.getAlleles(), record.getGenotypes(), record.getAttributes(), startSplitReadSites, endSplitReadSites,
                discordantPairs, copyNumberDistribution);
    }

    public SVCallRecordWithEvidence(final String id,
                                    final String startContig,
                                    final int position1,
                                    final boolean strand1,
                                    final String contig2,
                                    final int position2,
                                    final boolean strand2,
                                    final StructuralVariantType type,
                                    final int length,
                                    final List<String> algorithms,
                                    final List<Allele> alleles,
                                    final List<Genotype> genotypes,
                                    final Map<String,Object> attributes,
                                    final List<SplitReadSite> startSplitReadSites,
                                    final List<SplitReadSite> endSplitReadSites,
                                    final List<DiscordantPairEvidence> discordantPairs,
                                    final CopyNumberPosteriorDistribution copyNumberDistribution) {
        super(id, startContig, position1, strand1, contig2, position2, strand2, type, length, algorithms, alleles, genotypes, attributes);
        Utils.nonNull(startSplitReadSites);
        Utils.nonNull(endSplitReadSites);
        Utils.nonNull(discordantPairs);
        this.startSplitReadSites = Collections.unmodifiableList(startSplitReadSites);
        this.endSplitReadSites = Collections.unmodifiableList(endSplitReadSites);
        this.discordantPairs = Collections.unmodifiableList(discordantPairs);
        this.copyNumberDistribution = copyNumberDistribution;
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

    public CopyNumberPosteriorDistribution getCopyNumberDistribution() { return copyNumberDistribution; }

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
        areEqual &= b.getCopyNumberDistribution() == this.getCopyNumberDistribution()
                || (this.getCopyNumberDistribution() != null && this.getCopyNumberDistribution().equals(b.getCopyNumberDistribution()));
        return areEqual;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), discordantPairs, endSplitReadSites, startSplitReadSites, copyNumberDistribution);
    }

}
