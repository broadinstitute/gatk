package org.broadinstitute.hellbender.tools.funcotator.metadata;

import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.utils.Utils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * A concrete class for {@link FuncotationMetadata} that can be easily built from a VCF Header.
 */
public class VcfFuncotationMetadata implements FuncotationMetadata {

    private LinkedHashMap<String, VCFInfoHeaderLine> fieldNameToHeaderLineMap;

    private VcfFuncotationMetadata(final LinkedHashMap<String, VCFInfoHeaderLine> fieldNameToHeaderLineMap) {
        this.fieldNameToHeaderLineMap = fieldNameToHeaderLineMap;
    }

    /**
     * @param vcfInfoHeaderLines Never {@code null}
     * @return Metadata corresponding to VCF info fields.  Never {@code null}
     */
    public static VcfFuncotationMetadata create(final List<VCFInfoHeaderLine> vcfInfoHeaderLines) {
        Utils.nonNull(vcfInfoHeaderLines);
        return new VcfFuncotationMetadata(
            vcfInfoHeaderLines.stream().collect(Collectors.toMap(v -> v.getID(), Function.identity(), (x1,x2) -> x2, LinkedHashMap::new ))
        );
    }

    /**
     * See {@link FuncotationMetadata#retrieveHeaderInfo(String)}
     *
     * @param fieldName field to search.  Never {@code null}
     * @return Never {@code null}
     */
    @Override
    public VCFInfoHeaderLine retrieveHeaderInfo(final String fieldName) {
        return fieldNameToHeaderLineMap.get(fieldName);
    }

    /**
     * See {@link FuncotationMetadata#retrieveAllHeaderInfo()}
     *
     * @return Never {@code null}
     */
    @Override
    public List<VCFInfoHeaderLine> retrieveAllHeaderInfo() {
        return new ArrayList<>(fieldNameToHeaderLineMap.values());
    }

    @Override
    public String toString() {
        return "VcfFuncotationMetadata{" +
                "fieldNameToHeaderLineMap=" + fieldNameToHeaderLineMap +
                '}';
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final VcfFuncotationMetadata that = (VcfFuncotationMetadata) o;

        return fieldNameToHeaderLineMap != null ? fieldNameToHeaderLineMap.equals(that.fieldNameToHeaderLineMap) : that.fieldNameToHeaderLineMap == null;
    }

    @Override
    public int hashCode() {
        return fieldNameToHeaderLineMap != null ? fieldNameToHeaderLineMap.hashCode() : 0;
    }
}
