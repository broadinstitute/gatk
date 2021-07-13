package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Enum for the attributes belonging to a VariantContext object
 *
 */

public enum VariantContextAttributeEnum {
    CONTIG("contig"),
    START("start"),
    END("end"),
    ID("ID"),
    ALLELES("alleles"),
    FILTERS("filters"),
    TYPE("type"),
    PHRED_SCALED_QUALITY("phredScaledQual"),
    ATTRIBUTES("attributes"),
    GENOTYPES("genotypes"),
    SOURCE("source");

    private String name;

    VariantContextAttributeEnum (String name){
        this.name = name;
    }

    public String getName(){
        return this.name;
    }

    public Object getValue(VariantContext vc){
        switch(this){
            case CONTIG:
                return vc.getContig();
            case START:
                return vc.getStart();
            case END:
                return vc.getEnd();
            case ID:
                return vc.getID();
            case ALLELES:
                return vc.getAlleles();
            case FILTERS:
                return vc.getFilters();
            case TYPE:
                return vc.getType();
            case PHRED_SCALED_QUALITY:
                return vc.getPhredScaledQual();
            case ATTRIBUTES:
                return vc.getAttributes();
            case GENOTYPES:
                return vc.getGenotypes();
            case SOURCE:
                return vc.getSource();
            default:
                return null;
        }
    }
}
