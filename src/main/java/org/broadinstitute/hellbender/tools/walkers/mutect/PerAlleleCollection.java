package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;

/**
 * A container for allele to value mapping.
 *
 * Each PerAlleleCollection may hold a value for each ALT allele and, optionally, a value for the REF allele.
 * For example,
 *
 *   PerAlleleCollection<Double> alleleFractions = PerAlleleCollection.createPerAltAlleleCollection()
 *
 * may be a container for allele fractions for ALT alleles in a variant context. While
 *
 *   PerAlleleCollection<Double> alleleCount = PerAlleleCollection.createPerRefAndAltAlleleCollection()
 *
 * may hold the allele counts for the REF allele and all ALT alleles in a variant context.
 *
 *
 **/
public class PerAlleleCollection<X> {
    // TODO: consider using Optional for ref allele
    private Optional<Allele> refAllele;
    private Optional<X> refValue;
    private Map<Allele, X> altAlleleValueMap;
    private Type type;

    public enum Type {ALT_ONLY, REF_AND_ALT}

    public PerAlleleCollection(final Type type){
        this.type = type;
        altAlleleValueMap = new HashMap<>();
        refAllele = Optional.empty();
    }

    /**
     * Take an allele, REF or ALT, and update its value appropriately
     *
     * @param allele : REF or ALT allele
     * @param value :
     */
    public void set(Allele allele, X value){
        Utils.nonNull(allele, "allele is null");
        Utils.nonNull(value, "value is null");
        Utils.validateArg(type == Type.REF_AND_ALT || allele.isNonReference(), "Collection stores values for alternate alleles only");
        if (allele.isReference()){
            setRef(allele, value);
        } else {
            setAlt(allele, value);
        }
    }

    public void set(final Collection<Allele> alleles, final Function<Allele, X> function) {
        alleles.forEach(a -> set(a, function.apply(a)));
    }

    public void setRef(Allele allele, X value){
        Utils.nonNull(allele, "ref allele is null");
        Utils.nonNull(value, "value is null");
        Utils.validateArg(allele.isReference(), "setting non-reference allele as reference");
        Utils.validateArg(!refAllele.isPresent(), "Resetting the reference allele not permitted");
        refAllele = Optional.of(allele);
        refValue = Optional.of(value);
    }

    public void setAlt(Allele allele, X value){
        Utils.nonNull(allele, "ref allele is null");
        Utils.nonNull(value, "value is null");
        Utils.validateArg(allele.isNonReference(), "Setting reference allele as alt");
        altAlleleValueMap.put(allele, value);
    }

    /**
     * Get the value for an allele, REF or ALT
     * @param allele
     */
    public X get(Allele allele){
        Utils.nonNull(allele, "allele is null");
        if (allele.isReference()){
            Utils.validateArg(allele.equals(refAllele.get()), "Requested ref allele does not match the stored ref allele");
            return getRef();
        } else {
            return getAlt(allele);
        }
    }

    public X getRef(){
        if (type == Type.ALT_ONLY) {
            throw new IllegalStateException("Collection does not hold the REF allele");
        }

        if (refAllele.isPresent()){
            return refValue.get();
        } else {
            throw new IllegalStateException("Collection's ref allele has not been set yet");
        }
    }

    public X getAlt(Allele allele){
        Utils.nonNull(allele, "allele is null");
        Utils.validateArg(allele.isNonReference(), "allele is not an alt allele");
        Utils.validateArg(altAlleleValueMap.containsKey(allele), "Requested alt allele is not in the collection");
        return altAlleleValueMap.get(allele);
    }
    
    public Set<Allele> getAltAlleles(){
        return altAlleleValueMap.keySet();
    }
}