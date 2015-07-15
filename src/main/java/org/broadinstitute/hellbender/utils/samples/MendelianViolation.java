package org.broadinstitute.hellbender.utils.samples;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Class for the identification and tracking of mendelian violation. It can be used in 2 distinct ways:
 * - Either using an instance of the MendelianViolation class to track mendelian violations for each of the families while
 * walking over the variants
 * - Or using the static methods to directly get information about mendelian violation in a family at a given locus
 *
 */
public final class MendelianViolation {

    //Stores occurrences of inheritance
    private EnumMap<GenotypeType, EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> inheritance;

    private int violations_total=0;

    private final boolean allCalledOnly = false;
    private final double minGenotypeQuality;
    private final boolean abortOnSampleNotFound;

    /**
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     */
    public MendelianViolation(final double minGenotypeQualityP) {
        this(minGenotypeQualityP, true);
    }

    /**
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     * @param abortOnSampleNotFound - Whether to stop execution if a family is passed but no relevant genotypes are found. If false, then the family is ignored.
     */
    public MendelianViolation(final double minGenotypeQualityP, final boolean abortOnSampleNotFound) {
        minGenotypeQuality = minGenotypeQualityP;
        this.abortOnSampleNotFound = abortOnSampleNotFound;
        createInheritanceMap();
    }

    //Count of violations of the type HOM_REF/HOM_REF -> HET
    public int getParentsRefRefChildHet() {
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
    }

    public boolean isViolation(final Sample mother, final Sample father, final Sample child, final VariantContext vc){
        violations_total=0;
        clearInheritanceMap();
        updateViolations(mother.getFamilyID(), mother.getID(), father.getID(), child.getID(),vc);
        return violations_total>0;
    }

    private void updateViolations(final String familyId, final String motherId, final String fatherId, final String childId, final VariantContext vc){
        final Genotype gMom = vc.getGenotype(motherId);
        final Genotype gDad = vc.getGenotype(fatherId);
        final Genotype gChild = vc.getGenotype(childId);

        if (gMom == null || gDad == null || gChild == null){
            if(abortOnSampleNotFound) {
                throw new IllegalArgumentException(String.format("Variant %s:%d: Missing genotypes for family %s: mom=%s dad=%s family=%s", vc.getContig(), vc.getStart(), familyId, motherId, fatherId, childId));
            } else
                return;
        }
        //Count No calls
        if(allCalledOnly && (!gMom.isCalled() || !gDad.isCalled() || !gChild.isCalled())){
            //no call
        }
        else if (!gMom.isCalled() && !gDad.isCalled() || !gChild.isCalled()){
            //no call
        }
        //Count lowQual. Note that if min quality is set to 0, even values with no quality associated are returned
        else if (minGenotypeQuality > 0 && (
            gMom.getGQ()   < minGenotypeQuality ||
            gDad.getGQ()   < minGenotypeQuality ||
            gChild.getGQ() < minGenotypeQuality )) {
            //no call
        } else{
            if(isViolation(gMom, gDad, gChild)){
                violations_total++;
            }
            final int count = inheritance.get(gMom.getType()).get(gDad.getType()).get(gChild.getType());
            inheritance.get(gMom.getType()).get(gDad.getType()).put(gChild.getType(),count+1);
        }
    }

    /**
     * Evaluate the genotypes of mom, dad, and child to detect Mendelian violations
     *
     * @param gMom
     * @param gDad
     * @param gChild
     * @return true if the three genotypes represent a Mendelian violation; false otherwise
     */
    public static boolean isViolation(final Genotype gMom, final Genotype gDad, final Genotype gChild) {
        if (gChild.isNoCall()){ //cannot possibly be a violation is child is no call
            return false;
        }
        if(gMom.isHomRef() && gDad.isHomRef() && gChild.isHomRef()) {
            return false;
        }

        //1 parent is no "call
        if(!gMom.isCalled()){
            return (gDad.isHomRef() && gChild.isHomVar()) || (gDad.isHomVar() && gChild.isHomRef());
        }
        else if(!gDad.isCalled()){
            return (gMom.isHomRef() && gChild.isHomVar()) || (gMom.isHomVar() && gChild.isHomRef());
        }
        //Both parents have genotype information
        final Allele childRef = gChild.getAlleles().get(0);
        return !(gMom.getAlleles().contains(childRef) && gDad.getAlleles().contains(gChild.getAlleles().get(1)) ||
            gMom.getAlleles().contains(gChild.getAlleles().get(1)) && gDad.getAlleles().contains(childRef));
    }

    private void createInheritanceMap(){
        inheritance = new EnumMap<>(GenotypeType.class);
        for(GenotypeType mType : GenotypeType.values()){
            inheritance.put(mType, new EnumMap<>(GenotypeType.class));
            for(GenotypeType dType : GenotypeType.values()){
                inheritance.get(mType).put(dType, new EnumMap<>(GenotypeType.class));
                for(GenotypeType cType : GenotypeType.values()){
                    inheritance.get(mType).get(dType).put(cType, 0);
                }
            }
        }
    }

    private void clearInheritanceMap(){
        for(GenotypeType mType : GenotypeType.values()){
            for(GenotypeType dType : GenotypeType.values()){
                for(GenotypeType cType : GenotypeType.values()){
                    inheritance.get(mType).get(dType).put(cType, 0);
                }
            }
        }
    }
}
