package org.broadinstitute.hellbender.tools.evoquer;

import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.collections4.SetUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

public class AlleleSubsettingUtilsForJointCalling {
    private static final Logger logger = LogManager.getLogger(AlleleSubsettingUtilsForJointCalling.class);

    /**
     *
     * @param unmergedCalls
     * @param maxAlleles
     * @return
     */
    public static List<VariantContext> subsetAlleles(List<VariantContext> unmergedCalls, int maxAlleles, Allele longestRefAllele) {

        Map<Allele, Integer> alleleSpecificQuals = new HashMap<>();
        List<VariantContext> callsWithExpandedAlleles = sumAllQualsAndExpandAlleles(unmergedCalls, alleleSpecificQuals, longestRefAllele);
        logger.info("unique alleles: " + alleleSpecificQuals.size() + alleleSpecificQuals.keySet());

        // choose the target alleles
        List<Allele> targetAlleles = alleleSpecificQuals.entrySet().stream()
                // get the highest values
                .sorted(Map.Entry.<Allele, Integer>comparingByValue().reversed()).limit(maxAlleles)
                .map(entry -> entry.getKey()).collect(Collectors.toList());
        logger.info(targetAlleles);

        Set<Allele> allelesToDrop = new HashSet<>(alleleSpecificQuals.keySet());
        allelesToDrop.removeAll(targetAlleles);

        // how do i do this without hardcoding?
        int[] relevantAlleleIndexes = new int[] {0,2};

        return callsWithExpandedAlleles.stream().map(vc -> {
            Set<Allele> allelesToSubset = SetUtils.intersection(new HashSet<>(vc.getAlleles()), allelesToDrop);
            if (allelesToSubset.isEmpty()) {
                return vc;
            } else {
                logger.info("subsetting alleles for sample: " + vc.getSampleNames());
                Genotype g =  vc.getGenotype(0);
                GenotypeBuilder gb = new GenotypeBuilder(g);

                List<Allele> galleles = g.getAlleles();
                List<Allele> updatedGAlleles = galleles.stream().map(gallele -> {
                    if (allelesToSubset.contains(gallele)) {
                        return Allele.NO_CALL;
                    } else {
                        return gallele;
                    }}).collect(Collectors.toList());
                gb.alleles(updatedGAlleles);
                // TODO use methods to return index mapping?
                GenotypeLikelihoods gll = GenotypeLikelihoods.fromPLs(g.getPL());
                double[] probs = gll.getAsVector();
                double [] modifiedProbs = new double[]{probs[0], MathUtils.log10SumLog10(probs[1], probs[3]), MathUtils.log10SumLog10(new double[]{probs[2], probs[4], probs[5]})};
                gb.PL(modifiedProbs);

                int[] ad = g.getAD();
                int[] modifiedAD = new int[]{ad[0], ad[2]};
                gb.AD(modifiedAD);

                ArrayList<Genotype> newGlist = new ArrayList<>();
                newGlist.add(gb.make());

                Map<String, Object> attrs = vc.getCommonInfo().getAttributes();
                Map<String, Object> modifiedAttrs = attrs.entrySet().stream()
                        .collect(Collectors.toMap(
                                Map.Entry::getKey,
                                entry -> remapAlleleAttributes(entry.getValue().toString(), relevantAlleleIndexes)));
                List<Allele> newAlleles = vc.getAlleles();
                // replace the current ref with the longest ref allele
                newAlleles.remove(0);
                newAlleles.add(0, longestRefAllele);
                newAlleles = newAlleles.stream().filter(allele -> !allelesToDrop.contains(allele)).collect(Collectors.toList());
                return new VariantContextBuilder(vc)
                        .alleles(newAlleles)
                        .attributes(modifiedAttrs)
                        .genotypes(GenotypesContext.create(newGlist)).make();
            }
        }).collect(Collectors.toList());
    }

    private static String remapAlleleAttributes(String values, int[] relevantAlleleIndexes) {
        // how do i make sure these are allele specific annotations only?
        List<String> alleleSpecificValues = AnnotationUtils.getAlleleLengthListOfString(values);
        final List<?> subsetList = alleleSpecificValues.size() > 0 ? ReferenceConfidenceVariantContextMerger.remapRLengthList(alleleSpecificValues, relevantAlleleIndexes)
                : Collections.nCopies(relevantAlleleIndexes.length, "");
        return AnnotationUtils.encodeAnyASList(subsetList);
    }

    /**
     *
     * @param variantContexts
     * @param alleleSpecificQuals
     * @param longestRefAllele
     * @return
     */
    private static List<VariantContext> sumAllQualsAndExpandAlleles(List<VariantContext> variantContexts, Map<Allele, Integer> alleleSpecificQuals, Allele longestRefAllele) {
        return variantContexts.stream().map(vc ->
                sumQualsAndExpandedAllelesForVC(
                        vc,
                        longestRefAllele,
                        alleleSpecificQuals)
        ).collect(Collectors.toList());
    }

    /**
     * We return a variant context that contains expanded alleles if applicable and add the
     * quals for each allele to the allele quals map
     * @param vc
     * @param longestRefAllele
     * @param alleleSpecificQuals in/out param that contains the summed quals for each allele seen
     * @return
     */
    private static VariantContext sumQualsAndExpandedAllelesForVC(final VariantContext vc, Allele longestRefAllele, Map<Allele, Integer> alleleSpecificQuals) {
        Map<Allele, Allele> alleleMappings = new HashMap<>();
        // if you try to create mappings with the same length an exception is thrown
        if (vc.getReference().length() < longestRefAllele.length()) {
            // this is where the magic happens!
            alleleMappings = GATKVariantContextUtils.createAlleleMapping(longestRefAllele, vc, alleleSpecificQuals.keySet());
        }

        List<Allele> alleles = vc.getAlleles();
        VariantContextBuilder vcb = null;
        List<Allele> updatedGenotypeAlleles = new ArrayList<>(vc.getGenotype(0).getAlleles());
        VariantContext result = vc;     // default the result

        String rawDataString = vc.getAttributeAsString("AS_QUALapprox", null);
        String[] rawDataPerAllele = null;
        if (rawDataString != null && !rawDataString.isEmpty()) {
            rawDataPerAllele = rawDataString.split(AnnotationUtils.ALLELE_SPECIFIC_SPLIT_REGEX);
        }

            // TODO separate this out - we still need to subset other alleles
        for (int i = 0; i < alleles.size(); i++) {
            final Allele allele = alleles.get(i);
            if (rawDataPerAllele != null) {
                final String alleleQual = rawDataPerAllele[i];  // should i check that this is actually an int?
                if (!allele.isNonRefAllele() && !alleleQual.isEmpty() && !alleleQual.equals(AnnotationUtils.MISSING_VALUE)) {
                    sumQuals(alleleSpecificQuals, alleleMappings.getOrDefault(allele, allele), Integer.parseInt(alleleQual));
                }
            }
            // now replace the original allele with the expanded allele in the genotype and vc
            if (i == 0) {
                Allele oldRef = vc.getReference();
                vcb = new VariantContextBuilder(vc);
                setVariantContextBuilderAllele(vcb, 0, longestRefAllele);
                if (updatedGenotypeAlleles.contains(oldRef)) {
                    Collections.replaceAll(updatedGenotypeAlleles, oldRef, longestRefAllele);
                }
            }
            else if (alleleMappings.containsKey(allele)) {
                Allele mappedAllele = alleleMappings.get(allele);

                // replace allele in new vc
                setVariantContextBuilderAllele(vcb, i, mappedAllele);

                // replace every instance of allele in the genotype alleles
                if (updatedGenotypeAlleles.contains(allele)) {
                    Collections.replaceAll(updatedGenotypeAlleles, allele, mappedAllele);
                }
            } else {
                logger.warn("no mapping found for allele: " + allele);
            }
        }
        if (vcb != null) {  // this should never be null
            GenotypeBuilder genotypeBuilder = new GenotypeBuilder(vc.getGenotype(0));
            genotypeBuilder.alleles(updatedGenotypeAlleles);
            vcb.genotypes(Collections.singletonList(genotypeBuilder.make()));
            result = vcb.make();
        } else {
            logger.warn("vcb is null! " + vc.getAlleles());
        }
        return result;
    }

    private static void setVariantContextBuilderAllele(VariantContextBuilder vcb, int index, Allele newAllele) {
        List<Allele> modifyableAlleles = vcb.getAlleles();
        modifyableAlleles.set(index, newAllele);
        vcb.alleles(modifyableAlleles);
    }


    private static void sumQuals(Map<Allele, Integer> alleleSpecificQuals, Allele mappedAllele, int alleleQual) {
        int currentQual = 0;
        // refactor to use map update or replace or something
        if (alleleSpecificQuals.containsKey(mappedAllele)) {
            currentQual = alleleSpecificQuals.get(mappedAllele);
        }
        alleleSpecificQuals.put(mappedAllele, currentQual + alleleQual);
    }

}
