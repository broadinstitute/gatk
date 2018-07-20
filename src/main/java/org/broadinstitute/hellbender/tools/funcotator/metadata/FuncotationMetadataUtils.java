package org.broadinstitute.hellbender.tools.funcotator.metadata;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;

public class FuncotationMetadataUtils {
    private FuncotationMetadataUtils() {}
    public final static String UNKNOWN_DESCRIPTION = "Unknown";

    /**
     * @param fieldNames Never {@code null}
     * @return {@link FuncotationMetadata} with values populated indicating that we do not really know the metadata.
     * And type is a String.  Never {@code null}
     */
    public static FuncotationMetadata createWithUnknownAttributes(final List<String> fieldNames) {
        Utils.nonNull(fieldNames);
        return VcfFuncotationMetadata.create(
                fieldNames.stream().map(f -> new VCFInfoHeaderLine(f, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, UNKNOWN_DESCRIPTION))
                        .collect(Collectors.toList())
        );
    }

    //TODO: Docs
    //TODO: tests
    public static FuncotationMetadata merge(final FuncotationMetadata funcotationMetadata1, final FuncotationMetadata funcotationMetadata2) {
        final LinkedHashSet<VCFInfoHeaderLine> rawMetadata = new LinkedHashSet<>(funcotationMetadata1.retrieveAllHeaderInfo());
        rawMetadata.addAll(funcotationMetadata2.retrieveAllHeaderInfo());
        return VcfFuncotationMetadata.create(new ArrayList<>(rawMetadata));
    }
}
