package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;
import java.util.stream.Collectors;

public class LocalAlleler {

    public static final String LAA = "LAA";
    public static final String LGT = "LGT";

    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc){
        // new LAA
        // GT -> LGT
        // PL -> LPL
        // AD -> LAD

        final Map<String, Object> localAttributes = new LinkedHashMap<>();

        //construct LAA
        final LinkedHashSet<Allele> localAlleles = getLocalAlleles(originalGenotype, vc);
        final List<Integer> localAlleleIndexes = vc.getAlleleIndices(localAlleles);
        //todo this is stupidly inefficient probably
        localAlleleIndexes.remove(0);
        localAttributes.put(LAA, localAlleleIndexes);

        List<Integer> localGenotypes = new ArrayList<>(localAlleles.size());
        List<Allele> localAlleleList = new ArrayList<>(localAlleles);
        for(final Allele allele: originalGenotype.getAlleles()){
            localGenotypes.add(localAlleleList.indexOf(allele));
        }

        String delimiter = originalGenotype.isPhased() ? Genotype.PHASED_ALLELE_SEPARATOR : Genotype.UNPHASED_ALLELE_SEPARATOR;
        String localGenotypesString = localGenotypes.stream().map(String::valueOf).collect(Collectors.joining(delimiter));
        localAttributes.put(LGT, localGenotypesString);

        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(originalGenotype);
        localAttributes.forEach(genotypeBuilder::attribute);
        return genotypeBuilder.make();
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

}
