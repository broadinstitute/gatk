package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import java.util.*;

/**
 * Filters out a record if the allele balance for heterozygotes is out of a defined range across all samples.
 * The threshold is set as the minimum fraction of the data drawn from the less-represented allele - e.g. 0.3 would
 * set that whichever allele has lower representation across all heterozygous individuals must account for at least 30% of the
 * total observations.
 */
public final class AlleleBalanceFilter implements VariantFilter {
    /** The filter string used for sites that fail the allele balance filter. */
    public static final String AB_FILTER = "AlleleBalance";

    private final double hetAlleleBalance;

    public AlleleBalanceFilter(final double hetAlleleBalance) {
        this.hetAlleleBalance = hetAlleleBalance;
    }

    private static class Counts { int allele1; int allele2; }

    @Override
    public List<VCFFilterHeaderLine> headerLines() {
        return CollectionUtil.makeList(new VCFFilterHeaderLine(AB_FILTER, "Heterozygote allele balance below required threshold."));
    }

    @Override
    public String filter(final VariantContext ctx) {
        if (ctx.getHetCount() == 0) return null;
        final Map<List<Allele>, Counts> countsMap = new HashMap<>();

        for (final Genotype gt : ctx.getGenotypesOrderedByName()) {
            if (gt.isNoCall() || !gt.isHet()) continue;

            final List<Allele> alleles = gt.getAlleles();
            Counts counts = countsMap.get(alleles);
            if (counts == null) {
                counts = new Counts();
                countsMap.put(alleles, counts);
            }

            counts.allele1 += gt.getAD()[ctx.getAlleleIndex(alleles.get(0))];
            counts.allele2 += gt.getAD()[ctx.getAlleleIndex(alleles.get(1))];
        }

        for (final Counts counts : countsMap.values()) {
            final int total = counts.allele1 + counts.allele2;
            if (total > 0 && Math.min(counts.allele1, counts.allele2) / (double) total < this.hetAlleleBalance) return AB_FILTER;
        }

        return null;
    }
}

