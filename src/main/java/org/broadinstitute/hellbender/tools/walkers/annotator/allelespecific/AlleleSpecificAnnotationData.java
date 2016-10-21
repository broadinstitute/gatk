package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class to encapsulate the raw data for allele-specific classes compatible with the ReducibleAnnotation interface
 * @param <T> the type of raw data to be stored for later annotation calculation
 */
public final class AlleleSpecificAnnotationData<T> extends ReducibleAnnotationData<T>{
    private final List<Allele> alleleList;
    private final Allele refAllele;

    public AlleleSpecificAnnotationData(final List<Allele> inputAlleles, final String inputData) {
        super(inputData);
        this.attributeMap = new HashMap<>();
        inputAlleles.forEach(a -> {attributeMap.put(a, null);});
        alleleList = inputAlleles;
        refAllele = alleleList.stream().filter(Allele::isReference).findAny().orElse(null);
        checkRefAlleles();
    }

    @Override
    public List<Allele> getAlleles() {return Collections.unmodifiableList(alleleList);}

    /**
     * Get the reference allele for this allele-specific data.
     * (Used in cases where annotations compare some attribute of the alt alleles to that of the reference.)
     * @return  the reference allele for this data
     */
    public Allele getRefAllele() {return refAllele;}

    @Override
    public void setAttributeMap(final Map<Allele, T> inputMap) {
        super.setAttributeMap(inputMap);
    }

    private void checkRefAlleles() {
        final long refCount = alleleList.stream().filter(Allele::isReference).count();
        if (refCount > 1) {
            throw new IllegalArgumentException("ERROR: multiple reference alleles found in annotation data\n");
        }
        if (refCount == 0) {
            throw new IllegalArgumentException("ERROR: no reference alleles found in annotation data\n");
        }
    }
}
