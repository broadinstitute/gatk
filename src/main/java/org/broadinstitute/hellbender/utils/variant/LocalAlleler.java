package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

import java.util.*;
import java.util.stream.Collectors;

public class LocalAlleler {

    public static final String LAA = "LAA";
    public static final String LGT = "LGT";
    public static final String LAD = "LAD";

    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc){
        // new LAA
        // GT -> LGT
        // PL -> LPL
        // AD -> LAD
        //todo handle originalGenotype has no GT field
        final Map<String, Object> localAttributes = new LinkedHashMap<>();

        //construct LAA
        final LinkedHashSet<Allele> localAllelesIncludingRef = getLocalAlleles(originalGenotype, vc);
        final List<Integer> localAlleleIndexes = convertToVCFRepresentationOfLocalAlleles(vc, localAllelesIncludingRef);
        localAttributes.put(LAA, localAlleleIndexes);

        //construct LGT
        List<Integer> localGenotypes = getLocalGenotypeList(originalGenotype, localAllelesIncludingRef);
        String localGenotypesString = createLocalGenotypeString(localGenotypes, originalGenotype.isPhased());
        localAttributes.put(LGT, localGenotypesString);

        //construct LAD
        if( originalGenotype.hasAD()){
            ArrayList<Integer> localAlleleDepth = new ArrayList<>(localAlleleIndexes.size());
            // all alleles A*, C,T, AAT, NON_REF  local A*,T,NON_REF  AD 10,11,12,13,14  LAD 10,12,14
            int[] originalAlleleDepth = originalGenotype.getAD();

            //be sure to add the ref AD
            localAlleleDepth.add(originalAlleleDepth[0]);
            for( int index : localAlleleIndexes){
                localAlleleDepth.add(originalAlleleDepth[index]);
            }
            localAttributes.put(LAD, localAlleleDepth);
        }

        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(originalGenotype);
        //add the new attributes to the builder
        localAttributes.forEach(genotypeBuilder::attribute);
        return genotypeBuilder.make();
    }

    private static List<Integer> convertToVCFRepresentationOfLocalAlleles(VariantContext vc, LinkedHashSet<Allele> localAllelesIncludingRef) {
        final List<Integer> localAlleleIndexes = vc.getAlleleIndices(localAllelesIncludingRef);
        //todo this is stupidly inefficient probably
        localAlleleIndexes.remove(0);
        return localAlleleIndexes;
    }

    private static List<Integer> getLocalGenotypeList(Genotype originalGenotype, LinkedHashSet<Allele> localAlleles) {
        List<Integer> localGenotypes = new ArrayList<>(localAlleles.size());
        List<Allele> localAlleleList = new ArrayList<>(localAlleles);
        for(final Allele allele: originalGenotype.getAlleles()){
            localGenotypes.add(localAlleleList.indexOf(allele));
        }
        return localGenotypes;
    }

    private static String createLocalGenotypeString(List<Integer> localGenotypes, boolean phased) {
        String delimiter = phased ? Genotype.PHASED_ALLELE_SEPARATOR : Genotype.UNPHASED_ALLELE_SEPARATOR;
        return localGenotypes.stream().map(String::valueOf).collect(Collectors.joining(delimiter));
    }

    private static LinkedHashSet<Allele> getLocalAlleles(Genotype originalGenotype, VariantContext vc) {
        final LinkedHashSet<Allele> localAlleles = new LinkedHashSet<>();
        //add the reference as the 0th allele always
        localAlleles.add(vc.getReference());
        localAlleles.addAll(originalGenotype.getAlleles());

        // if there is a NON_REF allele, add it
        Allele lastAllele = vc.getAlleles().get(vc.getNAlleles() - 1);
        if(lastAllele.isNonRefAllele()){
            localAlleles.add(Allele.NON_REF_ALLELE);
        }
        return localAlleles;
    }

    private static class LocalAlleleMapper {
        public LocalAlleleMapper(Genotype originalGenotype, VariantContext variantContext){
            LinkedHashSet<Allele> localAlleles = getLocalAlleles(originalGenotype, variantContext);

        }
    }

}
