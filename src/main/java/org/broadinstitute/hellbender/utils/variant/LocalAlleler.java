package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;

import java.util.*;
import java.util.stream.Collectors;

public class LocalAlleler {

    public static final String LAA = "LAA";
    public static final String LGT = "LGT";
    public static final String LAD = "LAD";
    public static final String LPL = "LPL";

    //todo Handle no-calls differently by merging all evidence into non-ref and reducing to 0/non_ref
    //todo Handle hom-ref in VCF by choosing the second most likely genotype and setting it as a local allele.
    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc) {
        return addLocalFields(originalGenotype, vc, true, false);
    }

    public static Genotype addLocalFields(Genotype originalGenotype, VariantContext vc, final boolean removeNonLocalVersions, final boolean allowMixedOutput){
        Utils.nonNull(originalGenotype);
        Utils.nonNull(vc);
        // new LAA
        // GT -> LGT
        // PL -> LPL
        // AD -> LAD

        //todo handle originalGenotype has no GT field
        final Map<String, Object> localAttributes = new LinkedHashMap<>();

        //construct LAA
        final AlleleListPermutation<Allele> localAllelesIncludingRef = getLocalAlleles(originalGenotype, vc);

        if(allowMixedOutput && localAllelesIncludingRef.isNonPermuted()){
            return originalGenotype;
        }

        final List<Integer> localAlleleIndexes = convertToVCFRepresentationOfLocalAlleles(localAllelesIncludingRef);
        localAttributes.put(LAA, localAlleleIndexes);

        //construct LGT
        final List<Integer> localGenotypes = getLocalGenotypeList(originalGenotype, localAllelesIncludingRef);
        final String localGenotypesString = createLocalGenotypeString(localGenotypes, originalGenotype.isPhased());
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
            final int[] originalPls = originalGenotype.getPL();
            final int[] plIndices = AlleleSubsettingUtils.subsettedPLIndices(originalGenotype.getPloidy(), localAllelesIncludingRef);
            final int[] newPls = new int[plIndices.length];
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
        final List<Integer> localGenotypes = new ArrayList<>(localAlleles.numberOfAlleles());
        for(final Allele allele: originalGenotype.getAlleles()){
            localGenotypes.add(localAlleles.indexOfAllele(allele));
        }
        return localGenotypes;
    }

    private static String createLocalGenotypeString(List<Integer> localGenotypes, boolean phased) {
        final String delimiter = phased ? Genotype.PHASED_ALLELE_SEPARATOR : Genotype.UNPHASED_ALLELE_SEPARATOR;
        return localGenotypes.stream()
                .map(genotypeIndex -> genotypeIndex == -1 ? "." : String.valueOf(genotypeIndex))
                .collect(Collectors.joining(delimiter));
    }

    private static AlleleListPermutation<Allele> getLocalAlleles(Genotype originalGenotype, VariantContext vc) {
        final IndexedAlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(vc.getAlleles());

        //For a no call (./.) we have no information about which alleles might be relevant so we have to include them all.
        //TODO decide if this is the right behavior
        if( originalGenotype.isNoCall()){
            return originalAlleleList.permutation();
        }

        final LinkedHashSet<Allele> localAlleles = new LinkedHashSet<>();
        //add the reference as the 0th allele always
        localAlleles.add(vc.getReference());
        localAlleles.addAll(originalGenotype.getAlleles());

        // if there is a NON_REF allele, add it
        if (vc.getAlleles().stream().anyMatch(Allele::isNonRefAllele)){
            localAlleles.add(Allele.NON_REF_ALLELE);
        }

        try {
            return originalAlleleList.permutation(new IndexedAlleleList<>(localAlleles));
        } catch (final Exception e){
            throw new RuntimeException("Error at " + vc.getStart(), e);
        }


    }
}
