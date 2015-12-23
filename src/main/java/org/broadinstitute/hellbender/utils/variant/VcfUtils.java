package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.vcf.VCFHeader;

import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * Utils for dealing with VCF files.
 */
public final class VcfUtils {

    private VcfUtils(){}

    public static SortedSet<String> getSortedSampleSet(Map<String, VCFHeader> headers, GATKVariantContextUtils.GenotypeMergeType mergeOption) {
        final SortedSet<String> samples = new TreeSet<>();
        for (final Map.Entry<String, VCFHeader> val : headers.entrySet()) {
            VCFHeader header = val.getValue();
            samples.addAll(header.getGenotypeSamples().stream().map(sample -> GATKVariantContextUtils.mergedSampleName(val.getKey(), sample,
                    mergeOption == GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY)).collect(Collectors.toList()));
        }

        return samples;
    }


}
