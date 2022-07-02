package org.broadinstitute.hellbender.testutils;


import htsjdk.variant.variantcontext.Genotype;
/**
 * Enum for the attributes belonging to a Genotype object
 */

public enum GenotypeAttributeEnum {
    SAMPLE_NAME("sampleName"),
    ALLELES("alleles"),
    IS_PHASED("isPhased"),
    DP("DP"),
    AD("AD"),
    GQ("GQ"),
    PL("PL"),
    TYPE("type"),
    FILTERS("filters"),
    GENOTYPE_STRING("genotypeString"),
    LIKELIHOODS("likelihoods"),
    PLOIDY("ploidy"),
    EXTENDED_ATTRIBUTES("extendedAttributes");

    private String name;

    GenotypeAttributeEnum (String name){
        this.name = name;
    }

    public String getName(){
        return this.name;
    }

    public Object getValue(Genotype genotype){
        switch(this){
            case SAMPLE_NAME:
                return genotype.getSampleName();
            case ALLELES:
                return genotype.getAlleles();
            case IS_PHASED:
                return genotype.isPhased();
            case DP:
                return genotype.getDP();
            case AD:
                return genotype.getAD();
            case GQ:
                return genotype.getGQ();
            case PL:
                return genotype.getPL();
            case TYPE:
                return genotype.getType();
            case FILTERS:
                return genotype.getFilters();
            case GENOTYPE_STRING:
                return genotype.getGenotypeString();
            case LIKELIHOODS:
                return genotype.getLikelihoods();
            case PLOIDY:
                return genotype.getPloidy();
            case EXTENDED_ATTRIBUTES:
                return genotype.getExtendedAttributes();
            default:
                return null;
        }
    }
}
