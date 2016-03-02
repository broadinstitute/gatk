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

    // as determined experimentally Nov-Dec 2013
    public final static GATKVCFIndexType DEFAULT_GVCF_INDEX_TYPE = GATKVCFIndexType.LINEAR;
    public final static Integer DEFAULT_GVCF_INDEX_PARAMETER = 128000;

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


    /**
     * Check if using the GCVF indexing arguments' values
     *
     * @param variantIndexType variant indexing strategy
     * @param variantIndexParameter variant indexing parameter
     * @return true if the index type and parameter are the default GVCF values, false otherwise
     */
    public static boolean usingGVCFIndexingArguments(final GATKVCFIndexType variantIndexType, final int variantIndexParameter) {
        return variantIndexType == DEFAULT_GVCF_INDEX_TYPE && variantIndexParameter == DEFAULT_GVCF_INDEX_PARAMETER;
    }
}

