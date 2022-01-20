package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/*
 * TODO this whole class needs refactoring. it's cleaned up significantly from VQSR version,
 *  but still has a long way to go. we should decide how strongly to couple to the tranche code and refactor
 *  that at the same time
 */
final class VariantDataManager {
    private static final Logger logger = LogManager.getLogger(VariantDataManager.class);

    private static final int INITIAL_SIZE = 1000000;

    private final List<VariantDatum> data;
    final List<List<Allele>> alternateAlleles;
    private final List<String> annotationKeys;
    private final List<VariantSet> variantSets;
    private final boolean useASannotations; //TODO maybe just pass this and the below as an argument collection
    private final boolean trustAllPolymorphic;

    VariantDataManager(final List<String> annotationKeys,
                       final boolean useASannotations,
                       final boolean trustAllPolymorphic) {
        this.data = new ArrayList<>(INITIAL_SIZE);
        this.alternateAlleles = new ArrayList<>(INITIAL_SIZE);
        final List<String> uniqueAnnotations = annotationKeys.stream().distinct().sorted().collect(Collectors.toList()); // sort alphabetically
        if (annotationKeys.size() != uniqueAnnotations.size()) {
            logger.warn(String.format("Ignoring duplicate annotations for recalibration %s.", Utils.getDuplicatedItems(annotationKeys)));
        }
        this.annotationKeys = new ArrayList<>( uniqueAnnotations );
        variantSets = new ArrayList<>();
        this.useASannotations = useASannotations;
        this.trustAllPolymorphic = trustAllPolymorphic;
    }

    void addDatum(final FeatureContext featureContext, 
                  final VariantContext vc, 
                  final Allele refAllele, 
                  final Allele altAllele,
                  final boolean isExtractTrainingAndTruthOnly) {
        final VariantDatum datum = new VariantDatum();

        // Populate the datum with lots of fields from the VariantContext, unfortunately the VC is too big so we just
        // pull in only the things we absolutely need.
        datum.referenceAllele = refAllele;
        datum.alternateAllele = altAllele;
        decodeAnnotations(datum, vc);

        // non-deterministic because order of calls depends on load of machine TODO SL: not sure what this means?
        datum.loc = new SimpleInterval(vc);

        datum.isBiallelicSNP = vc.isSNP() && vc.isBiallelic();
        datum.isTransition = datum.isBiallelicSNP && GATKVariantContextUtils.isTransition(vc);

        // Loop through the training data sets and if they overlap this locus (and allele, if applicable) then update
        // the training status appropriately. The locus used to find training set variants is retrieved
        // by parseTrainingSets from the FeatureContext argument.
        parseTrainingSets(featureContext, vc, datum, useASannotations, trustAllPolymorphic);

        if (!isExtractTrainingAndTruthOnly || (datum.atTrainingSite || datum.atTruthSite)) {
            data.add(datum);

            // TODO clean this up; need to store alleles to enable passing of training sites-only VCF to score tool for marking of training sites
            if (isExtractTrainingAndTruthOnly) {
                alternateAlleles.add(vc.getAlternateAlleles());
            }
        }
    }

    public List<VariantDatum> getData() {
        return data;
    }

    void addVariantSet(final VariantSet variantSet) {
        variantSets.add(variantSet);
    }

    List<String> getAnnotationKeys() {
        return annotationKeys;
    }

    boolean checkHasTrainingSet() {
        for (final VariantSet variantSet : variantSets) {
            if (variantSet.isTraining) {
                return true;
            }
        }
        return false;
    }

    boolean checkHasTruthSet() {
        for (final VariantSet variantSet : variantSets) {
            if (variantSet.isTruth) {
                return true;
            }
        }
        return false;
    }

    private void decodeAnnotations(final VariantDatum datum,
                                   final VariantContext vc) {
        final double[] annotations = new double[annotationKeys.size()];
        final boolean[] isNull = new boolean[annotationKeys.size()];
        int iii = 0;
        for(final String key : annotationKeys) {
            isNull[iii] = false;
            annotations[iii] = decodeAnnotation(key, vc, useASannotations, datum);
            if( Double.isNaN(annotations[iii]) ) { isNull[iii] = true; }
            iii++;
        }
        datum.annotations = annotations;
        datum.isNull = isNull;
    }

    private static double decodeAnnotation(final String annotationKey,
                                           final VariantContext vc,
                                           final boolean useASannotations,
                                           final VariantDatum datum ) {
        double value;

        try {
            //if we're in allele-specific mode and an allele-specific annotation has been requested, parse the appropriate value from the list
            if (useASannotations && annotationKey.startsWith(GATKVCFConstants.ALLELE_SPECIFIC_PREFIX)) {
                final List<Object> valueList = vc.getAttributeAsList(annotationKey);
                //FIXME: we need to look at the ref allele here too
                if (vc.hasAllele(datum.alternateAllele)) {
                    final int altIndex = vc.getAlleleIndex(datum.alternateAllele) - 1; //- 1 is to convert the index from all alleles (including reference) to just alternate alleles
                    value = Double.parseDouble((String) valueList.get(altIndex));
                } else {
                    //if somehow our alleles got mixed up
                    throw new IllegalStateException("ExtractAnnotationsVariantDatum allele " + datum.alternateAllele + " is not contained in the input VariantContext.");
                }
            } else {
                value = vc.getAttributeAsDouble(annotationKey, Double.NaN);
            }
            if (Double.isInfinite(value)) {
                value = Double.NaN;
            }
        } catch (final NumberFormatException e) {
            value = Double.NaN;
        }
        return value;
    }

    private void parseTrainingSets(final FeatureContext featureContext,
                                   final VariantContext evalVC,
                                   final VariantDatum datum,
                                   final boolean useASannotations,
                                   final boolean trustAllPolymorphic) {
        datum.atTruthSite = false;
        datum.atTrainingSite = false;

        for (final VariantSet variantSet : variantSets) {
            final List<VariantContext> vcs = featureContext.getValues(variantSet.variantSource, featureContext.getInterval().getStart());
            for (final VariantContext trainVC : vcs) {
                if (useASannotations && !doAllelesMatch(trainVC, datum)) {
                    continue;
                }
                if (isValidVariant(evalVC, trainVC, trustAllPolymorphic)) {
                    datum.atTruthSite = datum.atTruthSite || variantSet.isTruth;
                    datum.atTrainingSite = datum.atTrainingSite || variantSet.isTraining;
                }
            }
        }
    }

    private static boolean isValidVariant(final VariantContext evalVC,
                                          final VariantContext trainVC,
                                          final boolean trustAllPolymorphic) {
        return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() && checkVariationClass(evalVC, trainVC) &&
                (trustAllPolymorphic || !trainVC.hasGenotypes() || trainVC.isPolymorphicInSamples());
    }

    private static boolean doAllelesMatch(final VariantContext trainVC,
                                          final VariantDatum datum) {
        //only do this check in the allele-specific case, where each datum represents one allele
        if (datum.alternateAllele == null) {
            return true;
        }
        try {
            return GATKVariantContextUtils.isAlleleInList(datum.referenceAllele, datum.alternateAllele, trainVC.getReference(), trainVC.getAlternateAlleles());
        } catch (final IllegalStateException e) {
            throw new IllegalStateException("Reference allele mismatch at position " + trainVC.getContig() + ':' + trainVC.getStart() + " : ", e);
        }
    }

    private static boolean checkVariationClass(final VariantContext evalVC,
                                               final VariantContext trainVC) {
        switch(trainVC.getType()) {
            case SNP:
            case MNP:
                return checkVariationClass(evalVC, VariantTypeMode.SNP);
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return checkVariationClass(evalVC, VariantTypeMode.INDEL);
            default:
                return false;
        }
    }

    static boolean checkVariationClass(final VariantContext vc,
                                       final VariantTypeMode mode) {
        switch(mode) {
            case SNP:
                return vc.isSNP() || vc.isMNP();
            case INDEL:
                return vc.isStructuralIndel() || vc.isIndel() || vc.isMixed() || vc.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new IllegalStateException("Encountered unknown mode: " + mode);
        }
    }

    static boolean checkVariationClass(final VariantContext vc,
                                       final Allele allele,
                                       final VariantTypeMode mode) {
        switch(mode) {
            case SNP:
                //note that spanning deletions are considered SNPs by this logic
                return vc.getReference().length() == allele.length();
            case INDEL:
                return (vc.getReference().length() != allele.length()) || allele.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new IllegalStateException("Encountered unknown mode: " + mode);
        }
    }

    static void setScores(final List<VariantDatum> data,
                          final double[] scores) {
        IntStream.range(0, data.size()).forEach(i -> data.get(i).score = scores[i]);
    }
}
