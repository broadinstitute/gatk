package org.broadinstitute.hellbender.tools.funcotator.metadata;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Utility class to determine the tumor normal pairs that makeup a VCF Header
 */
public class SamplePairExtractor {

    public static List<String> TUMOR_SAMPLE_NAME_LIST = Arrays.asList("TUMOR", "CASE", "MET");
    public static List<String> NORMAL_SAMPLE_NAME_LIST = Arrays.asList("NORMAL", "CONTROL");
    public static String NO_NORMAL = "";

    private SamplePairExtractor(){}

    /**
     * This method should not be needed in germline use cases, but should produce an empty list in such cases.
     *
     * Pairing is done as follows:
     *
     * - Using the M2 convention:  "tumor_sample" in the header and "normal_sample" if applicable
     * - If there is only one sample, the pair will be the sample name, "".
     * - The samples have names such as "TUMOR" and "NORMAL".  See {@link SamplePairExtractor#TUMOR_SAMPLE_NAME_LIST} and
     *     {@link SamplePairExtractor#NORMAL_SAMPLE_NAME_LIST}
     *
     * If none of these are fulfilled, the returned list will be empty.
     *
     * @param vcfHeader Never {@code null}.  Use blank constructor in {@link VCFHeader} rather than null.
     * @return Tumor-normal sample names.  In tumor-only cases, the normal sample will have an empty string.
     *  Never {@code null}
     */
    public static List<TumorNormalPair> extractPossibleTumorNormalPairs(final VCFHeader vcfHeader) {
        Utils.nonNull(vcfHeader);
        // First attempt the M2 nomenclature
        if (vcfHeader.getOtherHeaderLineUnique(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER) != null) {
            final VCFHeaderLine normalMetaDataLine = vcfHeader.getOtherHeaderLineUnique(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER);
            return Collections.singletonList(new TumorNormalPair(vcfHeader.getOtherHeaderLineUnique(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER).getValue(),
                    (normalMetaDataLine == null ?
                            NO_NORMAL : normalMetaDataLine.getValue())));
        }

        // Then try sample names (e.g. "TUMOR", "NORMAL")
        final List<String> sampleNames = vcfHeader.getSampleNamesInOrder();
        if (sampleNames.size() == 1) {
            return Collections.singletonList(new TumorNormalPair(sampleNames.get(0), NO_NORMAL));
        }
        if (sampleNames.size() > 0) {
            final List<String> tumorSamples = sampleNames.stream()
                    .filter(s -> TUMOR_SAMPLE_NAME_LIST.contains(s)).collect(Collectors.toList());
            final List<String> normalSamples = sampleNames.stream()
                    .filter(s -> NORMAL_SAMPLE_NAME_LIST.contains(s)).collect(Collectors.toList());

            return tumorSamples.stream().flatMap(t ->
                    normalSamples.stream().map(n -> new TumorNormalPair(t,n))).collect(Collectors.toList());
        }

        return Collections.emptyList();
    }
}
