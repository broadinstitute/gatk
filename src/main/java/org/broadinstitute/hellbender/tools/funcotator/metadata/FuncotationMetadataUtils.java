package org.broadinstitute.hellbender.tools.funcotator.metadata;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.List;
import java.util.stream.Collectors;

public class FuncotationMetadataUtils {
    private FuncotationMetadataUtils() {}
    public final static String UNKNOWN_DESCRIPTION = "Unknown";

    /**
     * @param fieldNames Never {@code null}
     * @return
     */
    public static FuncotationMetadata createWithUnknownAttributes(final List<String> fieldNames) {
        return VcfFuncotationMetadata.create(
                fieldNames.stream().map(f -> new VCFInfoHeaderLine(f, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, UNKNOWN_DESCRIPTION))
                        .collect(Collectors.toList())
        );
    }
}
