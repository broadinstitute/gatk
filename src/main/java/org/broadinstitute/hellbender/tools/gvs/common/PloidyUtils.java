package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Lazy;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class PloidyUtils {
    /**
     * This copies the PAR definition used in the tool FindMendelianViolations with the addition of the Y regions
     * https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/vcf/MendelianViolations/FindMendelianViolations.java#L203
     * Wikipedia data: https://en.wikipedia.org/wiki/Pseudoautosomal_region
     */
    private static final Set<String> PSEUDO_AUTOSOMAL_REGIONS_DEFINITION = CollectionUtil .makeSet("X:60001-2699520", "X:154931044-155260560", "chrX:10001-2781479", "chrX:155701383-156030895", "Y:10001-2649520", "Y:59034050-59363566", "chrY:10001-2781479", "chrY:56887903-57217415");

    public static Lazy<Set<Interval>> ParIntervals = new Lazy<>(() -> Collections.unmodifiableSet(parseIntervalLists(PSEUDO_AUTOSOMAL_REGIONS_DEFINITION)));

    public static boolean doesVariantOverlapPAR(VariantContext variant) {
        return doesIntervalOverlapPAR(variant.getContig(), variant.getStart(), variant.getEnd());
    }

    public static boolean isLocationInPAR(long location) {
        // Note: SchemaUtils.decodeContig makes implicit hg38 assumptions. If we support different references, we'll
        // need to factor that into how we do our mapping from pure location numbers to contigs and positions
        return doesIntervalOverlapPAR(SchemaUtils.decodeContig(location), SchemaUtils.decodePosition(location), SchemaUtils.decodePosition(location) + 1);
    }

    public static boolean doesIntervalOverlapPAR(String chrom, int positionStart, int positionEnd) {
        // Note: SchemaUtils.decodeContig makes implicit hg38 assumptions. If we support different references, we'll
        // need to factor that into how we do our mapping from pure location numbers to contigs and positions
        Interval locationAsInterval = new Interval(chrom, positionStart, positionEnd);
        for (Interval par : ParIntervals.get()) {
            if (par.intersects(locationAsInterval)) {
                return true;
            }
        }
        return false;
    }

    private static Set<Interval> parseIntervalLists(final Set<String> intervalLists){
        // Parse the PARs
        final Set<Interval> intervals = new HashSet<>(intervalLists.size());
        for (final String par : intervalLists) {
            final String[] splits1 = par.split(":");
            final String[] splits2 = splits1[1].split("-");
            intervals.add(new Interval(splits1[0], Integer.parseInt(splits2[0]), Integer.parseInt(splits2[1])));
        }
        return intervals;
    }

}
