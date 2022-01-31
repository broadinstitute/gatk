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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/*
 * TODO this whole class needs refactoring. it's cleaned up significantly from VQSR version,
 *  but still has a long way to go. we should decide how strongly to couple to the tranche code and refactor
 *  that at the same time
 */
final class VariantDataCollection {
    private static final Logger logger = LogManager.getLogger(VariantDataCollection.class);

    private static final int INITIAL_SIZE = 1000000;
    private static final String BIALLELIC_SNP_LABEL = "biallelic_snp";
    private static final String TRANSITION_LABEL = "transition";
    static final List<String> COMPUTED_LABELS = Collections.unmodifiableList(Arrays.asList(BIALLELIC_SNP_LABEL, TRANSITION_LABEL));

    private final List<VariantDatum> data;
    final List<List<Allele>> alternateAlleles;
    private final List<String> annotationKeys;
    private final boolean useASannotations; //TODO maybe just pass this and the below as an argument collection
    private final boolean trustAllPolymorphic;

    VariantDataCollection(final List<String> annotationKeys,
                          final boolean useASannotations,
                          final boolean trustAllPolymorphic) {
        this.data = new ArrayList<>(INITIAL_SIZE);
        this.alternateAlleles = new ArrayList<>(INITIAL_SIZE);
        final List<String> uniqueAnnotations = annotationKeys.stream().distinct().sorted().collect(Collectors.toList()); // sort alphabetically
        if (annotationKeys.size() != uniqueAnnotations.size()) {
            logger.warn(String.format("Ignoring duplicate annotations for recalibration %s.", Utils.getDuplicatedItems(annotationKeys)));
        }
        this.annotationKeys = new ArrayList<>( uniqueAnnotations );
        this.useASannotations = useASannotations;
        this.trustAllPolymorphic = trustAllPolymorphic;
    }

    void addDatum(final VariantContext vc,
                  final Allele refAllele, 
                  final Allele altAllele,
                  final Set<String> labels) {
        final VariantDatum datum = new VariantDatum();

        // Populate the datum with lots of fields from the VariantContext, unfortunately the VC is too big so we just
        // pull in only the things we absolutely need.
        datum.referenceAllele = refAllele;
        datum.alternateAllele = altAllele;
        decodeAnnotations(datum, vc);

        // non-deterministic because order of calls depends on load of machine TODO SL: not sure what this means?
        datum.loc = new SimpleInterval(vc);

        datum.labels = labels;
        if (vc.isSNP() && vc.isBiallelic()) {
            datum.labels.add(BIALLELIC_SNP_LABEL);
            if (GATKVariantContextUtils.isTransition(vc)) {
                datum.labels.add(TRANSITION_LABEL);
            }
        }

        data.add(datum);
    }

    // TODO clean this up; need to store alleles to enable passing of training sites-only VCF to score tool for marking of training sites
    void addAlternateAlleles(final List<Allele> alternateAlleles) {
        this.alternateAlleles.add(alternateAlleles);
    }

    public List<VariantDatum> getData() {
        return data;
    }

    List<String> getAnnotationKeys() {
        return annotationKeys;
    }

    private void decodeAnnotations(final VariantDatum datum,
                                   final VariantContext vc) {
        final double[] annotations = new double[annotationKeys.size()];
        int iii = 0;
        for(final String key : annotationKeys) {
            annotations[iii] = decodeAnnotation(key, vc, useASannotations, datum);
            iii++;
        }
        datum.annotations = annotations;
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

    static void setScores(final List<VariantDatum> data,
                          final double[] scores) {
        IntStream.range(0, data.size()).forEach(i -> data.get(i).score = scores[i]);
    }
}
