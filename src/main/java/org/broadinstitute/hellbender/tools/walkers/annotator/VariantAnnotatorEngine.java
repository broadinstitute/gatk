package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException.BadArgumentValue;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.reflections.ReflectionUtils;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * The class responsible for computing annotations for variants.
 * Annotations are auto-discovered - ie, any class that extends {@link VariantAnnotation} and
 * lives in this package is treated as an annotation and the engine will attempt to create instances of it
 * by calling the non-arg constructor (loading will fail if there is no no-arg constructor).
 */
public final class VariantAnnotatorEngine {
    private final List<InfoFieldAnnotation> infoAnnotations;
    private final List<GenotypeAnnotation> genotypeAnnotations;

    private final VariantOverlapAnnotator variantOverlapAnnotator;

    private VariantAnnotatorEngine(final AnnotationManager annots,
                                   final FeatureInput<VariantContext> dbSNPInput,
                                   final List<FeatureInput<VariantContext>> featureInputs){
        infoAnnotations = annots.createInfoFieldAnnotations();
        genotypeAnnotations = annots.createGenotypeAnnotations();
        variantOverlapAnnotator = initializeOverlapAnnotator(dbSNPInput, featureInputs);
    }

    /**
     * Makes the engine for all known annotation types (minus the excluded ones).
     * @param annotationsToExclude list of annotations to exclude (pass an empty list to indicate that there are no exclusions)
     * @param dbSNPInput input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param comparisonFeatureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     */
    public static VariantAnnotatorEngine ofAllMinusExcluded(final List<String> annotationsToExclude,
                                                            final FeatureInput<VariantContext> dbSNPInput,
                                                            final List<FeatureInput<VariantContext>> comparisonFeatureInputs) {
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        Utils.nonNull(comparisonFeatureInputs, "comparisonFeatureInputs is null");
        return new VariantAnnotatorEngine(AnnotationManager.ofAllMinusExcluded(annotationsToExclude), dbSNPInput, comparisonFeatureInputs);
    }

    /**
     * Makes the engine for given annotation types and groups (minus the excluded ones).
     * @param annotationGroupsToUse list of annotations groups to include
     * @param annotationsToUse     list of of annotations to include
     * @param annotationsToExclude list of annotations to exclude
     * @param dbSNPInput input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param comparisonFeatureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     */
    public static VariantAnnotatorEngine ofSelectedMinusExcluded(final List<String> annotationGroupsToUse,
                                                                 final List<String> annotationsToUse,
                                                                 final List<String> annotationsToExclude,
                                                                 final FeatureInput<VariantContext> dbSNPInput,
                                                                 final List<FeatureInput<VariantContext>> comparisonFeatureInputs) {
        Utils.nonNull(annotationGroupsToUse, "annotationGroupsToUse is null");
        Utils.nonNull(annotationsToUse, "annotationsToUse is null");
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        Utils.nonNull(comparisonFeatureInputs, "comparisonFeatureInputs is null");
        return new VariantAnnotatorEngine(AnnotationManager.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude), dbSNPInput, comparisonFeatureInputs);
    }

    private VariantOverlapAnnotator initializeOverlapAnnotator(final FeatureInput<VariantContext> dbSNPInput, final List<FeatureInput<VariantContext>> featureInputs) {
        final Map<FeatureInput<VariantContext>, String> overlaps = new LinkedHashMap<>();
        for ( final FeatureInput<VariantContext> fi : featureInputs) {
            overlaps.put(fi, fi.getName());
        }
        if (overlaps.values().contains(VCFConstants.DBSNP_KEY)){
            throw new GATKException("The map of overlaps must not contain " + VCFConstants.DBSNP_KEY);
        }
        if (dbSNPInput != null) {
            overlaps.put(dbSNPInput, VCFConstants.DBSNP_KEY); // add overlap detection with DBSNP by default
        }

        return new VariantOverlapAnnotator(dbSNPInput, overlaps);
    }

    /**
     * Returns the list of genotype annotations that will be applied.
     * Note: The returned list is unmodifiable.
     */
    public List<GenotypeAnnotation> getGenotypeAnnotations() {
        return Collections.unmodifiableList(genotypeAnnotations);
    }

    /**
     * Returns the list of info annotations that will be applied.
     * Note: The returned list is unmodifiable.
     */
    public List<InfoFieldAnnotation> getInfoAnnotations() {
        return Collections.unmodifiableList(infoAnnotations);
    }

    /**
     * Returns the set of descriptions to be added to the VCFHeader line (for all annotations in this engine).
     */
    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {
        final Set<VCFHeaderLine> descriptions = new LinkedHashSet<>();

        for ( final InfoFieldAnnotation annotation : infoAnnotations) {
            descriptions.addAll(annotation.getDescriptions());
        }
        for ( final GenotypeAnnotation annotation : genotypeAnnotations) {
            descriptions.addAll(annotation.getDescriptions());
        }
        for ( final String db : variantOverlapAnnotator.getOverlapNames() ) {
            if ( VCFStandardHeaderLines.getInfoLine(db, false) != null ) {
                descriptions.add(VCFStandardHeaderLines.getInfoLine(db));
            } else {
                descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, db + " Membership"));
            }
        }

        if (descriptions.contains(null)){
            throw new IllegalStateException("getVCFAnnotationDescriptions should not contain null. This error is likely due to an incorrect implementation of getDescriptions() in one or more of the annotation classes");
        }
        return descriptions;
    }

    /**
     * Annotates the given variant context - adds all annotations that satisfy the predicate.
     * @param vc the variant context to annotate
     * @param features context containing the features that overlap the given variant
     * @param ref the reference context of the variant to annotate or null if there is none
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample. May be null
     * @param addAnnot function that indicates if the given annotation type should be added to the variant
     */
    public VariantContext annotateContext(final VariantContext vc,
                                          final FeatureContext features,
                                          final ReferenceContext ref,
                                          final ReadLikelihoods<Allele> likelihoods,
                                          final Predicate<VariantAnnotation> addAnnot) {
        Utils.nonNull(vc, "vc cannot be null");
        Utils.nonNull(features, "features cannot be null");

        // annotate genotypes, creating another new VC in the process
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.genotypes(annotateGenotypes(ref, vc, likelihoods, addAnnot));
        final VariantContext newGenotypeAnnotatedVC = builder.make();

        final Map<String, Object> infoAnnotMap = new LinkedHashMap<>(newGenotypeAnnotatedVC.getAttributes());
        for ( final InfoFieldAnnotation annotationType : this.infoAnnotations) {
            if (addAnnot.test(annotationType)){
                final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(ref, newGenotypeAnnotatedVC, likelihoods);
                if ( annotationsFromCurrentType != null ) {
                    infoAnnotMap.putAll(annotationsFromCurrentType);
                }
            }
        }

        // create a new VC with info and genotype annotations
        final VariantContext annotated = builder.attributes(infoAnnotMap).make();

        // annotate db occurrences
        return variantOverlapAnnotator.annotateOverlaps(features, variantOverlapAnnotator.annotateRsID(features, annotated));
    }

    private GenotypesContext annotateGenotypes(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods,
                                               final Predicate<VariantAnnotation> addAnnot) {
        if ( genotypeAnnotations.isEmpty() ) {
            return vc.getGenotypes();
        }

        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            final GenotypeBuilder gb = new GenotypeBuilder(genotype);
            for ( final GenotypeAnnotation annotation : genotypeAnnotations) {
                if (addAnnot.test(annotation)) {
                    annotation.annotate(ref, vc, genotype, gb, likelihoods);
                }
            }
            genotypes.add(gb.make());
        }

        return genotypes;
    }

    private static final class AnnotationManager {

        private final List<String> annotationGroupsToUse;
        private final List<String> annotationsToUse;
        private final List<String> annotationsToExclude;

        private AnnotationManager(final List<String> annotationGroupsToUse, final List<String> annotationsToUse, final List<String> annotationsToExclude){
            this.annotationGroupsToUse = annotationGroupsToUse;
            this.annotationsToUse = annotationsToUse;
            this.annotationsToExclude = annotationsToExclude;

            final Set<String> allAnnotationNames = new LinkedHashSet<>(AnnotationManager.getAllAnnotationNames());
            final Set<String> unknownAnnots = Sets.difference(new LinkedHashSet<>(annotationsToUse), allAnnotationNames);
            assertAnnotationExists(unknownAnnots);

            final Set<String> unknownAnnotsExclude = Sets.difference(new LinkedHashSet<>(annotationsToExclude), allAnnotationNames);
            assertAnnotationExists(unknownAnnotsExclude);

            final Set<String> unknownGroups =  Sets.difference(new LinkedHashSet<>(annotationGroupsToUse), new LinkedHashSet<>(AnnotationManager.getAllAnnotationGroupNames()));
            if (!unknownGroups.isEmpty()){
                throw new BadArgumentValue("group", "Unknown annotation group in " + unknownGroups + ". Known groups are " + AnnotationManager.getAllAnnotationGroupNames());
            }
        }

        private void assertAnnotationExists(final Set<String> missingAnnots){
            if (!missingAnnots.isEmpty()){
                throw new BadArgumentValue("annotation", "Annotation " + missingAnnots + " not found; please check that you have specified the name correctly");
            }
        }

        /**
         * An annotation will be included only when:
         * - it is in one of the annotation groups or
         * - it is listed explicitly
         * - and it is not excluded explicitly.
         */
        static AnnotationManager ofSelectedMinusExcluded(final List<String> annotationGroupsToUse, final List<String> annotationsToUse, final List<String> annotationsToExclude){
            final List<String> groups = new ArrayList<>(annotationGroupsToUse);//make copy
            final List<String> annots = new ArrayList<>(annotationsToUse);//make copy
            final List<String> excludes = new ArrayList<>(annotationsToExclude);//make copy
            return new AnnotationManager(groups, annots, excludes);
        }

        /**
         * An annotation will be included only when it is not excluded explicitly.
         */
        static AnnotationManager ofAllMinusExcluded(final List<String> annotationsToExclude){
            final List<String> groups = getAllAnnotationGroupNames();
            final List<String> annots = getAllAnnotationNames();
            return new AnnotationManager(groups, annots, annotationsToExclude);
        }

        private static List<String> getAllAnnotationNames() {
            final Set<VariantAnnotation> union = Sets.union(new LinkedHashSet<>(makeAllGenotypeAnnotations()), new LinkedHashSet<>(AnnotationManager.makeAllInfoFieldAnnotations()));
            return union.stream().map(a -> a.getClass().getSimpleName()).collect(Collectors.toList());
        }

        /**
         * Annotation group names are simple names of all interfaces that are subtypes of {@ Annotation}.
         */
        public static List<String> getAllAnnotationGroupNames() {
            return ClassUtils.knownSubInterfaceSimpleNames(Annotation.class);
        }

        public List<InfoFieldAnnotation> createInfoFieldAnnotations() {
            final List<InfoFieldAnnotation> all = makeAllInfoFieldAnnotations();
            return filterAnnotations(all);
        }

        private static List<InfoFieldAnnotation> makeAllInfoFieldAnnotations() {
            return ClassUtils.makeInstancesOfSubclasses(InfoFieldAnnotation.class, Annotation.class.getPackage());
        }

        public List<GenotypeAnnotation> createGenotypeAnnotations() {
            final List<GenotypeAnnotation> all = makeAllGenotypeAnnotations();
            return filterAnnotations(all);
        }

        private static List<GenotypeAnnotation> makeAllGenotypeAnnotations() {
            return ClassUtils.makeInstancesOfSubclasses(GenotypeAnnotation.class, Annotation.class.getPackage());
        }

        /**
         * Returns a list of annotations that either:
         *  - belong to at least one of the requested annotation groups
         *  - belong to the set of requested annotations
         *
         *  - and are NOT listed for exclusion
         *
         *  The list is sorted by simple name of the class.
         */
        private <T extends VariantAnnotation> List<T> filterAnnotations(final List<T> all) {
            final SortedSet<T> annotations = new TreeSet<>(Comparator.comparing(t -> t.getClass().getSimpleName()));

            final Set<Class<?>> knownAnnotationGroups = ClassUtils.knownSubInterfaces(Annotation.class);

            for (final T t : all){
                if (!annotationsToExclude.contains(t.getClass().getSimpleName())) {
                    //if any group matches requested groups, it's in
                    @SuppressWarnings("unchecked")
                    final Set<Class<?>> annotationGroupsForT = ReflectionUtils.getAllSuperTypes(t.getClass(), sup -> sup.isInterface() && knownAnnotationGroups.contains(sup));
                    if (annotationGroupsForT.stream().anyMatch(iface -> annotationGroupsToUse.contains(iface.getSimpleName()))) {
                        annotations.add(t);
                    } else if (annotationsToUse.contains(t.getClass().getSimpleName())) {
                        annotations.add(t);
                    }
                }
            }

            return Collections.unmodifiableList(new ArrayList<>(annotations));
        }

    }
    
}
