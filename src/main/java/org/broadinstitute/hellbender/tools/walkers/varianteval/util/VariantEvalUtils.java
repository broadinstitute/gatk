package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.StandardEval;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.RequiredStratification;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.StandardStratification;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.reflections.Reflections;

import java.util.*;

public class VariantEvalUtils {
    private final VariantEval variantEvalWalker;
    Logger logger;

    public VariantEvalUtils(VariantEval variantEvalWalker) {
        this.variantEvalWalker = variantEvalWalker;
        this.logger = variantEvalWalker.getLogger();
    }

    private final static Map<String, Class<? extends VariantStratifier>> stratifierClasses;
    private final static Set<String> standardStratificationNames;
    private final static Set<String> requiredStratificationNames;

    private final static Map<String, Class<? extends VariantEvaluator>> evaluatorClasses;
    private final static Set<String> standardEvaluatorNames;

    static {
        stratifierClasses = new HashMap<>();
        standardStratificationNames = new HashSet<>();
        requiredStratificationNames = new HashSet<>();

        Reflections reflectionsStrat = new Reflections("org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications");
        Set<Class<? extends VariantStratifier>> allClasses = reflectionsStrat.getSubTypesOf(VariantStratifier.class);
        for (Class<? extends VariantStratifier> clazz : allClasses) {
            stratifierClasses.put(clazz.getSimpleName(), clazz);

            if (StandardStratification.class.isAssignableFrom(clazz)) {
                standardStratificationNames.add(clazz.getSimpleName());
            }

            if (RequiredStratification.class.isAssignableFrom(clazz)) {
                requiredStratificationNames.add(clazz.getSimpleName());
            }
        }

        evaluatorClasses = new HashMap<>();
        standardEvaluatorNames= new HashSet<>();

        Reflections reflectionsEval = new Reflections("org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators");
        Set<Class<? extends VariantEvaluator>> allEvalClasses = reflectionsEval.getSubTypesOf(VariantEvaluator.class);
        for (Class<? extends VariantEvaluator> clazz : allEvalClasses) {
            evaluatorClasses.put(clazz.getSimpleName(), clazz);

            if (StandardEval.class.isAssignableFrom(clazz)) {
                standardEvaluatorNames.add(clazz.getSimpleName());
            }
        }
    }

    /**
     * List all of the available evaluation modules, then exit successfully
     */
    public void listModulesAndExit() {
        logger.info("Available stratification modules:");
        logger.info("(Standard modules are starred)");
        for (String name: stratifierClasses.keySet()) {

            logger.info("\t" + name + (requiredStratificationNames.contains(name) || standardStratificationNames.contains(name) ? "*" : ""));
        }
        logger.info("");

        logger.info("Available evaluation modules:");
        logger.info("(Standard modules are starred)");
        for (String veName : evaluatorClasses.keySet()) {
            logger.info("\t" + veName + (standardEvaluatorNames.contains(veName) ? "*" : ""));
        }
        logger.info("");

        System.exit(0);
    }

    /**
     * Initialize required, standard and user-specified stratification objects
     *
     * @param noStandardStrats  don't use the standard stratifications
     * @param modulesToUse      the list of stratification modules to use
     * @return set of stratifications to use
     */
    public List<VariantStratifier> initializeStratificationObjects(boolean noStandardStrats, List<String> modulesToUse) {
        TreeSet<VariantStratifier> strats = new TreeSet<>();
        Set<String> stratsToUse = new HashSet<>(requiredStratificationNames);

        // By default, use standard stratification modules.
        if (!noStandardStrats) {
            stratsToUse.addAll(standardStratificationNames);
        }

        // Now add the user-selected modules
        stratsToUse.addAll(modulesToUse);

        // Instantiate the stratifications
        for (String module : stratsToUse) {
            if (!stratifierClasses.containsKey(module)) {
                throw new CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if (stratifierClasses.containsKey(module)) {
                Class<? extends VariantStratifier> c = stratifierClasses.get(module);

                try {
                    VariantStratifier vs = c.newInstance();
                    vs.setVariantEvalWalker(variantEvalWalker);
                    vs.initialize();

                    strats.add(vs);
                } catch (InstantiationException e) {
                    throw new GATKException("Unable to instantiate stratification module '" + c.getSimpleName() + "'");
                } catch (IllegalAccessException e) {
                    throw new GATKException("Illegal access error when trying to instantiate stratification module '" + c.getSimpleName() + "'");
                }
            }
        }

        return new ArrayList<>(strats);
    }

    /**
     * Initialize required, standard and user-specified evaluation objects
     *
     * @param noStandardEvals don't use the standard evaluations
     * @param modulesToUse    the list of evaluation modules to use
     * @return set of evaluations to use
     */
    public Set<Class<? extends VariantEvaluator>> initializeEvaluationObjects(boolean noStandardEvals, List<String> modulesToUse) {
        Set<String> evalsToUse = new TreeSet<>(modulesToUse);

        // By default, use standard eval modules.
        if (!noStandardEvals) {
            evalsToUse.addAll(standardEvaluatorNames);
        }

        // Get the specific classes provided.
        Set<Class<? extends VariantEvaluator>> evals = new HashSet<>();
        for (String module : evalsToUse) {
            if (!evaluatorClasses.containsKey(module)) {
                throw new CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            evals.add(evaluatorClasses.get(module));
        }

        //add MetricsCollection if required modules are included
        if(evals.contains(evaluatorClasses.get("CompOverlap")) && evals.contains(evaluatorClasses.get("IndelSummary")) && evals.contains(evaluatorClasses.get("TiTvVariantEvaluator")) && evals.contains(evaluatorClasses.get("CountVariants")) && evals.contains(evaluatorClasses.get("MultiallelicSummary")) )
            evals.add(evaluatorClasses.get("MetricsCollection"));

        return evals;
    }

    /**
     * Subset a VariantContext to a single sample
     *
     * @param vc         the VariantContext object containing multiple samples
     * @param sampleName the sample to pull out of the VariantContext
     * @return a new VariantContext with just the requested sample
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, String sampleName) {
        return getSubsetOfVariantContext(vc, Collections.singleton(sampleName));
    }

    /**
     * Subset a VariantContext to a set of samples
     *
     * @param vc          the VariantContext object containing multiple samples
     * @param sampleNames the samples to pull out of the VariantContext
     * @return a new VariantContext with just the requested samples
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, Set<String> sampleNames) {
        // if we want to preserve AC0 sites as polymorphic we need to not rederive alleles
        final boolean deriveAlleles = variantEvalWalker.ignoreAC0Sites();
        return ensureAnnotations(vc, vc.subContextFromSamples(sampleNames, deriveAlleles));
    }

    public VariantContext ensureAnnotations(final VariantContext vc, final VariantContext vcsub) {
        final int originalAlleleCount = vc.getHetCount() + 2 * vc.getHomVarCount();
        final int newAlleleCount = vcsub.getHetCount() + 2 * vcsub.getHomVarCount();
        final boolean isSingleton = originalAlleleCount == newAlleleCount && newAlleleCount == 1;
        final boolean hasChrCountAnnotations = vcsub.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) &&
                vcsub.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) &&
                vcsub.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY);

        if ( ! isSingleton && hasChrCountAnnotations ) {
            // nothing to update
            return vcsub;
        } else {
            // have to do the work
            VariantContextBuilder builder = new VariantContextBuilder(vcsub);

            if ( isSingleton )
                builder.attribute(VariantEval.IS_SINGLETON_KEY, true);

            if ( ! hasChrCountAnnotations )
                VariantContextUtils.calculateChromosomeCounts(builder, true);

            return builder.make();
        }
    }

    /**
     * For a list of track names, bind the variant contexts to a trackName->sampleName->VariantContext mapping.
     * Additional variant contexts per sample are automatically generated and added to the map unless the sample name
     * matches the ALL_SAMPLE_NAME constant.
     *
     * @return the mapping of track to VC list that should be populated
     */
    public HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>>
    bindVariantContexts(ReferenceContext referenceContext,
                        FeatureContext featureContext,
                        List<FeatureInput<VariantContext>> tracks,
                        boolean byFilter,
                        boolean subsetBySample,
                        boolean trackPerSample,
                        boolean trackPerFamily,
                        boolean mergeTracks) {
        HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> bindings = new HashMap<>();

        FeatureInput<VariantContext> firstTrack = tracks.isEmpty() ? null : tracks.get(0);
        for (FeatureInput<VariantContext> track : tracks) {
            HashMap<String, Collection<VariantContext>> mapping = new HashMap<>();

            //TODO: limiting to only those w/ the same start is GATK3 behavior.
            for (VariantContext vc : featureContext.getValues(track, referenceContext.getInterval().getStart())) {

                // First, filter the VariantContext to represent only the samples for evaluation
                VariantContext vcsub = vc;

                if ((subsetBySample) && vc.hasGenotypes())
                    vcsub = getSubsetOfVariantContext(vc, variantEvalWalker.getSampleNamesForEvaluation());

                //always add a mapping for all samples together
                if ((byFilter || !vcsub.isFiltered())) {
                    addMapping(mapping, VariantEval.getAllSampleName(), vcsub);
                }

                // Now, if stratifying, split the subsetted vc per sample and add each as a new context
                if (vc.hasGenotypes() && trackPerSample) {
                    for (String sampleName : variantEvalWalker.getSampleNamesForEvaluation()) {
                        VariantContext samplevc = getSubsetOfVariantContext(vc, sampleName);

                        if (byFilter || !samplevc.isFiltered()) {
                            addMapping(mapping, sampleName, samplevc);
                        }
                    }
                }
                else if (vc.hasGenotypes() && trackPerFamily) {
                    for (final String familyName : variantEvalWalker.getFamilyNamesForEvaluation()) {
                        Set<String> familyMemberNames = new HashSet<>();
                        //if the current stratification family name is "all", then add all the families to the VC for evaluation here
                        if (familyName.equals(VariantEval.getAllFamilyName())) {
                            familyMemberNames = variantEvalWalker.getSampleNamesForEvaluation();
                        }
                        else {
                            Set<Sample> familyMembers = variantEvalWalker.getSampleDB().getFamily(familyName);
                            for (final Sample s : familyMembers) {
                                familyMemberNames.add(s.getID());
                            }
                        }
                        VariantContext samplevc = getSubsetOfVariantContext(vc, familyMemberNames);

                        if (byFilter || !samplevc.isFiltered()) {
                            addMapping(mapping, familyName, samplevc);
                        }
                    }
                }
            }

            if (mergeTracks && bindings.containsKey(firstTrack)) {
                // go through each binding of sample -> value and add all of the bindings from this entry
                HashMap<String, Collection<VariantContext>> firstMapping = bindings.get(firstTrack);
                for (Map.Entry<String, Collection<VariantContext>> elt : mapping.entrySet()) {
                    Collection<VariantContext> firstMappingSet = firstMapping.get(elt.getKey());
                    if (firstMappingSet != null) {
                        firstMappingSet.addAll(elt.getValue());
                    } else {
                        firstMapping.put(elt.getKey(), elt.getValue());
                    }
                }
            } else {
                bindings.put(track, mapping);
            }
        }

        return bindings;
    }

    private void addMapping(HashMap<String, Collection<VariantContext>> mappings, String sample, VariantContext vc) {
        if (!mappings.containsKey(sample))
            mappings.put(sample, new ArrayList<>(1));
        mappings.get(sample).add(vc);
    }
}