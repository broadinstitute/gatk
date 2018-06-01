package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

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
    private Set<String> reducibleKeys;
    private List<VAExpression> expressions = new ArrayList<>();

    private final VariantOverlapAnnotator variantOverlapAnnotator;
    private boolean expressionAlleleConcordance;
    private final boolean useRawAnnotations;

    private final static Logger logger = LogManager.getLogger(VariantAnnotatorEngine.class);

    /**
     * Creates an annotation engine from a list of selected annotations output from command line parsing
     * @param annotationList list of annotation objects (with any parameters already filled) to include
     * @param dbSNPInput input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param featureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     * @param useRaw When this is set to true, the annotation engine will call {@link ReducibleAnnotation#annotateRawData(ReferenceContext, VariantContext, ReadLikelihoods)}
     *               on annotations that extend {@link ReducibleAnnotation}, instead of {@link InfoFieldAnnotation#annotate(ReferenceContext, VariantContext, ReadLikelihoods)},
     *               which is the default for all annotations.
     */
    public VariantAnnotatorEngine(final Collection<Annotation> annotationList,
                                   final FeatureInput<VariantContext> dbSNPInput,
                                   final List<FeatureInput<VariantContext>> featureInputs,
                                   final boolean useRaw){
        Utils.nonNull(featureInputs, "comparisonFeatureInputs is null");
        infoAnnotations = new ArrayList<>();
        genotypeAnnotations = new ArrayList<>();
        for (Annotation annot : annotationList) {
            if (annot instanceof InfoFieldAnnotation) {
                infoAnnotations.add((InfoFieldAnnotation) annot);
            }
            if (annot instanceof GenotypeAnnotation) {
                genotypeAnnotations.add((GenotypeAnnotation) annot);
            }
        }
        variantOverlapAnnotator = initializeOverlapAnnotator(dbSNPInput, featureInputs);
        reducibleKeys = new HashSet<>();
        useRawAnnotations = useRaw;
        for (InfoFieldAnnotation annot : infoAnnotations) {
            if (annot instanceof ReducibleAnnotation) {
                reducibleKeys.add(((ReducibleAnnotation) annot).getRawKeyName());
            }
        }
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
        return getVCFAnnotationDescriptions(false);
    }

    /**
     * Returns the set of descriptions to be added to the VCFHeader line (for all annotations in this engine).
     * @param useRaw Whether to prefer reducible annotation raw key descriptions over their normal descriptions
     */
    public Set<VCFHeaderLine> getVCFAnnotationDescriptions(boolean useRaw) {
        final Set<VCFHeaderLine> descriptions = new LinkedHashSet<>();

        for ( final InfoFieldAnnotation annotation : infoAnnotations) {
            if (annotation instanceof ReducibleAnnotation && useRaw) {
                descriptions.addAll(((ReducibleAnnotation)annotation).getRawDescriptions());
            } else {
                descriptions.addAll(annotation.getDescriptions());
            }
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

        // Add header lines corresponding to the expression target files
        for (final VariantAnnotatorEngine.VAExpression expression : getRequestedExpressions()) {
            // special case the ID field
            if (expression.fieldName.equals("ID")) {
                descriptions.add(new VCFInfoHeaderLine(expression.fullName, 1, VCFHeaderLineType.String, "ID field transferred from external VCF resource"));
            } else {
                final VCFInfoHeaderLine targetHeaderLine = ((VCFHeader) new FeatureDataSource<>(expression.binding, 100, VariantContext.class).getHeader())
                        .getInfoHeaderLines().stream()
                        .filter(l -> l.getID().equals(expression.fieldName))
                        .findFirst().orElse(null);

                VCFInfoHeaderLine lineToAdd;
                if (targetHeaderLine != null) {
                    expression.sethInfo(targetHeaderLine);
                    if (targetHeaderLine.getCountType() == VCFHeaderLineCount.INTEGER) {
                        lineToAdd = new VCFInfoHeaderLine(expression.fullName, targetHeaderLine.getCount(), targetHeaderLine.getType(), targetHeaderLine.getDescription());
                    } else {
                        lineToAdd = new VCFInfoHeaderLine(expression.fullName, targetHeaderLine.getCountType(), targetHeaderLine.getType(), targetHeaderLine.getDescription());
                    }
                } else {
                    lineToAdd = new VCFInfoHeaderLine(expression.fullName, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Value transferred from another external VCF resource");
                    logger.warn(String.format("The requested expression attribute \"%s\" is missing from the header in its resource file %s", expression.fullName, expression.binding.getName()));
                }
                descriptions.add(lineToAdd);
                expression.sethInfo(lineToAdd);
            }
        }

        Utils.validate(!descriptions.contains(null), "getVCFAnnotationDescriptions should not contain null. This error is likely due to an incorrect implementation of getDescriptions() in one or more of the annotation classes");
        return descriptions;
    }

    /**
     * Combine (raw) data for reducible annotations (those that use raw data in gVCFs)
     * Mutates annotationMap by removing the annotations that were combined
     *
     * Additionally, will combine other annotations by parsing them as numbers and reducing them
     * down to the
     * @param allelesList   the list of merged alleles across all variants being combined
     * @param annotationMap attributes of merged variant contexts -- is modifying by removing successfully combined annotations
     * @return  a map containing the keys and raw values for the combined annotations
     */
    @SuppressWarnings({"unchecked"})
    public Map<String, Object> combineAnnotations(final List<Allele> allelesList, Map<String, List<?>> annotationMap) {
        Map<String, Object> combinedAnnotations = new HashMap<>();

        // go through all the requested reducible info annotationTypes
        for (final InfoFieldAnnotation annotationType : infoAnnotations) {
            if (annotationType instanceof ReducibleAnnotation) {
                ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;
                if (annotationMap.containsKey(currentASannotation.getRawKeyName())) {
                    final List<ReducibleAnnotationData<?>> annotationValue = (List<ReducibleAnnotationData<?>>) annotationMap.get(currentASannotation.getRawKeyName());
                    final Map<String, Object> annotationsFromCurrentType = currentASannotation.combineRawData(allelesList, annotationValue);
                    combinedAnnotations.putAll(annotationsFromCurrentType);
                    //remove the combined annotations so that the next method only processes the non-reducible ones
                    annotationMap.remove(currentASannotation.getRawKeyName());
                }
            }
        }
        return combinedAnnotations;
    }

    /**
     * Finalize reducible annotations (those that use raw data in gVCFs)
     * @param vc    the merged VC with the final set of alleles, possibly subset to the number of maxAltAlleles for genotyping
     * @param originalVC    the merged but non-subset VC that contains the full list of merged alleles
     * @return  a VariantContext with the final annotation values for reducible annotations
     */
    public VariantContext finalizeAnnotations(VariantContext vc, VariantContext originalVC) {
        final Map<String, Object> variantAnnotations = new LinkedHashMap<>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for (final InfoFieldAnnotation annotationType : infoAnnotations) {
            if (annotationType instanceof ReducibleAnnotation) {

                ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;

                final Map<String, Object> annotationsFromCurrentType = currentASannotation.finalizeRawData(vc, originalVC);
                if (annotationsFromCurrentType != null) {
                    variantAnnotations.putAll(annotationsFromCurrentType);
                }
                //clean up raw annotation data after annotations are finalized
                variantAnnotations.remove(currentASannotation.getRawKeyName());
            }
        }

        // generate a new annotated VC
        final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(variantAnnotations);

        // annotate genotypes, creating another new VC in the process
        final VariantContext annotated = builder.make();
        return annotated;
    }

    /**
     * Annotates the given variant context - adds all annotations that satisfy the predicate.
     * @param vc the variant context to annotate
     * @param features context containing the features that overlap the given variant
     * @param ref the reference context of the variant to annotate or null if there is none
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample. May be null
     * @param addAnnot function that indicates if the given annotation type should be added to the variant
     *
     */
    public VariantContext annotateContext(final VariantContext vc,
                                          final FeatureContext features,
                                          final ReferenceContext ref,
                                          final ReadLikelihoods<Allele> likelihoods,
                                          final Predicate<VariantAnnotation> addAnnot) {
        Utils.nonNull(vc, "vc cannot be null");
        Utils.nonNull(features, "features cannot be null");
        Utils.nonNull(addAnnot, "addAnnot cannot be null");

        // annotate genotypes, creating another new VC in the process
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.genotypes(annotateGenotypes(ref, vc, likelihoods, addAnnot));
        final VariantContext newGenotypeAnnotatedVC = builder.make();

        final Map<String, Object> infoAnnotMap = new LinkedHashMap<>(newGenotypeAnnotatedVC.getAttributes());
        annotateExpressions(vc, features, ref, infoAnnotMap);

        for ( final InfoFieldAnnotation annotationType : this.infoAnnotations) {
            if (addAnnot.test(annotationType)){
                final Map<String, Object> annotationsFromCurrentType;
                if (useRawAnnotations && annotationType instanceof ReducibleAnnotation) {
                    annotationsFromCurrentType = ((ReducibleAnnotation) annotationType).annotateRawData(ref, newGenotypeAnnotatedVC, likelihoods);
                } else {
                    annotationsFromCurrentType = annotationType.annotate(ref, newGenotypeAnnotatedVC, likelihoods);
                }
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

    /**
     * Method which checks if a key is a raw key of the requested reducible annotations
     * @param key annotation key to check
     * @return true if the key is the raw key for a requested annotation
     */
    public boolean isRequestedReducibleRawKey(String key) {
        return reducibleKeys.contains(key);
    }

    /**
     * A container object for storing the objects necessary for carrying over expression annotations.
     * It holds onto the source feature input as well as any relevant header lines in order to alter the vcfHeader.
     */
    public static class VAExpression {

        final private String fullName, fieldName;
        final private FeatureInput<VariantContext> binding;
        private VCFInfoHeaderLine hInfo;

        public VAExpression(String fullExpression, List<FeatureInput<VariantContext>> dataSourceList){
            final int indexOfDot = fullExpression.lastIndexOf(".");
            if ( indexOfDot == -1 ) {
                throw new UserException.BadInput("The requested expression '"+fullExpression+"' is invalid, it should be in VCFFile.value format");
            }

            fullName = fullExpression;
            fieldName = fullExpression.substring(indexOfDot+1);

            final String bindingName = fullExpression.substring(0, indexOfDot);
            Optional<FeatureInput<VariantContext>> binding = dataSourceList.stream().filter(ds -> ds.getName().equals(bindingName)).findFirst();
            if (!binding.isPresent()) {
                throw new UserException.BadInput("The requested expression '"+fullExpression+"' is invalid, could not find vcf input file");
            }
            this.binding = binding.get();
        }

        public void sethInfo(VCFInfoHeaderLine hInfo) {
            this.hInfo = hInfo;
        }
    }

    protected List<VAExpression> getRequestedExpressions() { return expressions; }

    // select specific expressions to use
    public void addExpressions(Set<String> expressionsToUse, List<FeatureInput<VariantContext>> dataSources, boolean expressionAlleleConcordance) {//, Set<VCFHeaderLines>) {
        // set up the expressions
        for ( final String expression : expressionsToUse ) {
            expressions.add(new VAExpression(expression, dataSources));
        }
        this.expressionAlleleConcordance = expressionAlleleConcordance;
    }

    /**
     * Handles logic for expressions for variant contexts. Used to add annotations from one vcf file into the fields
     * of another if the variant contexts match sufficiently between the two files.
     *
     * @param vc  VariantContext to add annotations to
     * @param features  FeatureContext object containing extra VCF features to add to vc
     * @param ref  Reference context object corresponding to the region overlapping vc
     * @param attributes  running list of attributes into which to place new annotations
     */
    private void annotateExpressions(final VariantContext vc,
                                     final FeatureContext features,
                                     final ReferenceContext ref,
                                     final Map<String, Object> attributes){
        Utils.nonNull(vc);

        // each requested expression
        for ( final VAExpression expression : expressions ) {
            List<VariantContext> variantContexts = features.getValues(expression.binding, vc.getStart());

            if (!variantContexts.isEmpty()) {
                // get the expression's variant context
                VariantContext expressionVC = variantContexts.iterator().next();

                // special-case the ID field
                if (expression.fieldName.equals("ID")) {
                    if (expressionVC.hasID()) {
                        attributes.put(expression.fullName, expressionVC.getID());
                    }
                } else if (expression.fieldName.equals("ALT")) {
                    attributes.put(expression.fullName, expressionVC.getAlternateAllele(0).getDisplayString());
                } else if (expression.fieldName.equals("FILTER")) {
                    final String filterString = expressionVC.isFiltered() ? expressionVC.getFilters().stream().collect(Collectors.joining(",")) : "PASS";
                    attributes.put(expression.fullName, filterString);
                } else if (expressionVC.hasAttribute(expression.fieldName)) {
                    // find the info field
                    final VCFInfoHeaderLine hInfo = expression.hInfo;
                    if (hInfo == null) {
                        throw new UserException("Cannot annotate expression " + expression.fullName + " at " + ref.getInterval() + " for variant allele(s) " + vc.getAlleles() + ", missing header info");
                    }

                    //
                    // Add the info field annotations
                    //
                    final boolean useRefAndAltAlleles = VCFHeaderLineCount.R == hInfo.getCountType();
                    final boolean useAltAlleles = VCFHeaderLineCount.A == hInfo.getCountType();

                    // Annotation uses ref and/or alt alleles or enforce allele concordance
                    if ((useAltAlleles || useRefAndAltAlleles) || expressionAlleleConcordance) {

                        // remove brackets and spaces from expression value
                        final String cleanedExpressionValue = expressionVC.getAttribute(expression.fieldName,"").toString().replaceAll("[\\[\\]\\s]", "");

                        // get comma separated expression values
                        final ArrayList<String> expressionValuesList = new ArrayList<>(Arrays.asList(cleanedExpressionValue.split(",")));

                        boolean canAnnotate = false;
                        // get the minimum biallelics without genotypes

                        final List<VariantContext> minBiallelicVCs = getMinRepresentationBiallelics(vc);
                        final List<VariantContext> minBiallelicExprVCs = getMinRepresentationBiallelics(expressionVC);

                        // check concordance
                        final List<String> annotationValues = new ArrayList<>();
                        for (final VariantContext biallelicVC : minBiallelicVCs) {
                            // check that ref and alt alleles are the same
                            List<Allele> exprAlleles = biallelicVC.getAlleles();
                            boolean isAlleleConcordant = false;
                            int i = 0;
                            for (final VariantContext biallelicExprVC : minBiallelicExprVCs) {
                                List<Allele> alleles = biallelicExprVC.getAlleles();
                                // concordant
                                if (alleles.equals(exprAlleles)) {
                                    // get the value for the reference if needed.
                                    if (i == 0 && useRefAndAltAlleles) {
                                        annotationValues.add(expressionValuesList.get(i++));
                                    }
                                    // use annotation expression and add to vc
                                    annotationValues.add(expressionValuesList.get(i));
                                    isAlleleConcordant = true;
                                    canAnnotate = true;
                                    break;
                                }
                                i++;
                            }

                            // can not find allele match so set to annotation value to zero
                            if (!isAlleleConcordant) {
                                annotationValues.add("0");
                            }
                        }

                        // some allele matches so add the annotation values
                        if (canAnnotate) {
                            attributes.put(expression.fullName, annotationValues);
                        }
                    } else {
                        // use all of the expression values
                        attributes.put(expression.fullName, expressionVC.getAttribute(expression.fieldName));
                    }
                }
            }
        }
    }

    /**
     * Break the variant context into bialleles (reference and alternate alleles) and trim to a minimum representation
     *
     * @param vc variant context to annotate
     * @return list of biallelics trimmed to a minimum representation
     */
    private List<VariantContext> getMinRepresentationBiallelics(final VariantContext vc) {
        final List<VariantContext> minRepresentationBiallelicVCs = new ArrayList<>();
        if (vc.getNAlleles() > 2) {
            // TODO, this doesn't actually need to be done, we can simulate it at less cost
            for (int i = 1; i < vc.getNAlleles(); i++) {
                // Determining if the biallelic would have been considered a SNP
                if (! (vc.getReference().length() == 1 && vc.getAlternateAllele(i-1).length() == 1) ) {
                    minRepresentationBiallelicVCs.add(GATKVariantContextUtils.trimAlleles(
                            new VariantContextBuilder(vc)
                                    .alleles(Arrays.asList(vc.getReference(),vc.getAlternateAllele(i-1)))
                                    .attributes(removeIrrelevantAttributes(vc.getAttributes())).make(), true, true));
                } else {
                    minRepresentationBiallelicVCs.add(new VariantContextBuilder(vc)
                            .alleles(Arrays.asList(vc.getReference(),vc.getAlternateAllele(i-1)))
                            .attributes(removeIrrelevantAttributes(vc.getAttributes())).make());
                }
            }
        } else {
            minRepresentationBiallelicVCs.add(vc);
        }

        return minRepresentationBiallelicVCs;
    }

    private Map<String, Object> removeIrrelevantAttributes(Map<String, Object> attributes) {
        // since the VC has been subset, remove the invalid attributes
        Map<String, Object> ret = new HashMap<>(attributes);
        for ( final String key : attributes.keySet() ) {
            if ( !(key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(VCFConstants.ALLELE_FREQUENCY_KEY) || key.equals(VCFConstants.ALLELE_NUMBER_KEY)) ) {
                ret.remove(key);
            }
        }

        return ret;
    }

}
