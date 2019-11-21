package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;

import java.util.*;
import java.util.stream.Collectors;

public class LocalAlleler {

    public static final String LAA = "LAA";
    public static final String LGT = "LGT";
    public static final String LAD = "LAD";
    public static final String LPL = "LPL";

    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc) {
        return addLocalFields(originalGenotype, vc, true);
    }

    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc, final boolean removeNonLocalVersions){
        // new LAA
        // GT -> LGT
        // PL -> LPL
        // AD -> LAD
        //todo handle originalGenotype has no GT field
        final Map<String, Object> localAttributes = new LinkedHashMap<>();

        //construct LAA
        final AlleleListPermutation<Allele> localAllelesIncludingRef = getLocalAlleles(originalGenotype, vc);
        final List<Integer> localAlleleIndexes = convertToVCFRepresentationOfLocalAlleles(localAllelesIncludingRef);
        localAttributes.put(LAA, localAlleleIndexes);

        //construct LGT
        List<Integer> localGenotypes = getLocalGenotypeList(originalGenotype, localAllelesIncludingRef);
        String localGenotypesString = createLocalGenotypeString(localGenotypes, originalGenotype.isPhased());
        localAttributes.put(LGT, localGenotypesString);

        //construct LAD
        if( originalGenotype.hasAD()){
            final ArrayList<Integer> localAlleleDepth = new ArrayList<>(localAllelesIncludingRef.numberOfAlleles());
            final int[] originalAlleleDepth = originalGenotype.getAD();

            for( int i = 0; i < localAllelesIncludingRef.numberOfAlleles(); i++) {
                localAlleleDepth.add(originalAlleleDepth[localAllelesIncludingRef.fromIndex(i)]);
            }
            localAttributes.put(LAD, localAlleleDepth);
        }

        //construct LPL
        if( originalGenotype.hasPL()) {
            int[] originalPls = originalGenotype.getPL();
            int[] plIndices = AlleleSubsettingUtils.subsettedPLIndices(originalGenotype.getPloidy(), localAllelesIncludingRef);
            int[] newPls = new int[plIndices.length];
            for (int i = 0; i < plIndices.length; i++) {
                newPls[i] = originalPls[plIndices[i]];
            }
            localAttributes.put(LPL, newPls);
        }

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(originalGenotype);

        if(removeNonLocalVersions){
            genotypeBuilder.noAD()
                    .noPL()
                    .alleles(Collections.emptyList());
        }

        //add the new attributes to the builder
        localAttributes.forEach(genotypeBuilder::attribute);
        return genotypeBuilder.make();
    }

    private static List<Integer> convertToVCFRepresentationOfLocalAlleles(AlleleListPermutation<Allele> localAllelesIncludingRef) {
        final List<Integer> localAlleleIndexes = new ArrayList<>(localAllelesIncludingRef.numberOfAlleles());
        for( int i = 1; i < localAllelesIncludingRef.numberOfAlleles(); i++){
            localAlleleIndexes.add(localAllelesIncludingRef.fromIndex(i));
        }
        return localAlleleIndexes;
    }

    private static List<Integer> getLocalGenotypeList(Genotype originalGenotype, AlleleListPermutation<Allele> localAlleles) {
        List<Integer> localGenotypes = new ArrayList<>(localAlleles.numberOfAlleles());
        for(final Allele allele: originalGenotype.getAlleles()){
            localGenotypes.add(localAlleles.indexOfAllele(allele));
        }
        return localGenotypes;
    }

    private static String createLocalGenotypeString(List<Integer> localGenotypes, boolean phased) {
        String delimiter = phased ? Genotype.PHASED_ALLELE_SEPARATOR : Genotype.UNPHASED_ALLELE_SEPARATOR;
        return localGenotypes.stream().map(String::valueOf).collect(Collectors.joining(delimiter));
    }

    private static AlleleListPermutation<Allele> getLocalAlleles(Genotype originalGenotype, VariantContext vc) {
        final LinkedHashSet<Allele> localAlleles = new LinkedHashSet<>();
        //add the reference as the 0th allele always
        localAlleles.add(vc.getReference());
        localAlleles.addAll(originalGenotype.getAlleles());

        // if there is a NON_REF allele, add it
        if (vc.getAlleles().stream().anyMatch(Allele::isNonRefAllele)){
            localAlleles.add(Allele.NON_REF_ALLELE);
        }
        final IndexedAlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(vc.getAlleles());//.permutation(new IndexedAlleleList<>(allelesToKeep);
        final AlleleListPermutation<Allele> localAllelesMapping = originalAlleleList.permutation(new IndexedAlleleList<>(localAlleles));

        return localAllelesMapping;
    }
}
