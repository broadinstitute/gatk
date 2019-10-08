package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

public class LocalAlleler {

    public static final String LAA = "LAA";

    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc){
        // new LAA
        // GT -> LGT
        // PL -> LPL
        // AD -> LAD

        final Map<String, Object> localAttributes = new LinkedHashMap<>();

        final LinkedHashSet<Allele> localAlleles = getLocalAlleles(originalGenotype, vc);
        final List<Integer> localAlleleIndexes = vc.getAlleleIndices(localAlleles);
        localAttributes.put(LAA, localAlleleIndexes);

        List<Integer> localGenotypes = new ArrayList<>(localAlleles.size());
        for(final Allele allele: originalGenotype.getAlleles()){
            List<Allele> localAlleleList = new ArrayList<>(localAlleles);
            localGenotypes.add(localAlleleList.indexOf(allele));
        }
        localAttributes.put("LGT", localGenotypes);
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
