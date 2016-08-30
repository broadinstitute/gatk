package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

public class SomaticGenotypingEngine extends HaplotypeCallerGenotypingEngine {

    private final M2ArgumentCollection MTAC;
    private final TumorPowerCalculator strandArtifactPowerCalculator;

    private final String tumorSampleName;
    private final String matchedNormalSampleName;
    private final String DEBUG_READ_NAME;

    //Mutect2 does not run in GGA mode
    private static final List<VariantContext> NO_GIVEN_ALLELES = Collections.emptyList();

    // {@link GenotypingEngine} requires a non-null {@link AFCalculatorProvider} but this class doesn't need it.  Thus we make a dummy
    private static AFCalculatorProvider DUMMY_AF_CALCULATOR_PROVIDER = new AFCalculatorProvider() {
        public AFCalculator getInstance(final int ploidy, final int maximumAltAlleles) { return null; }
    };

    private final static Logger logger = Logger.getLogger(SomaticGenotypingEngine.class);

    public SomaticGenotypingEngine(final SampleList samples,
                                   final boolean doPhysicalPhasing,
                                   final M2ArgumentCollection MTAC,
                                   final String tumorSampleName,
                                   final String matchedNormalSampleName,
                                   final String DEBUG_READ_NAME) {
        super(MTAC, samples, DUMMY_AF_CALCULATOR_PROVIDER, doPhysicalPhasing);
        this.MTAC = MTAC;
        this.tumorSampleName = tumorSampleName;
        this.matchedNormalSampleName = matchedNormalSampleName;
        this.DEBUG_READ_NAME = DEBUG_READ_NAME;

        // coverage related initialization
        //TODO: in GATK4, use a QualityUtils method
        final double errorProbability = Math.pow(10, -MTAC.POWER_CONSTANT_QSCORE/10);
        strandArtifactPowerCalculator = new TumorPowerCalculator(errorProbability, MTAC.STRAND_ARTIFACT_LOD_THRESHOLD, 0.0f);
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param readLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     *
     * //TODO: do something about the enormous amount of code duplication between this method and
     * //TODO: HaplotypeCallerGenotypingEngine::assignGenotypeLikelihoods
     */
    public CalledHaplotypes callMutations (
            final ReadLikelihoods<Haplotype> readLikelihoods,
            final Map<String, Integer> originalNormalReadQualities,
            final Map<String, List<GATKRead>> perSampleFilteredReadList,
            final byte[] ref,
            final SimpleInterval refLoc,
            final SimpleInterval activeRegionWindow,
            final FeatureContext featureContext,
            final SAMFileHeader header) {
        Utils.nonNull(readLikelihoods, "likelihoods are null");
        Utils.validateArg(readLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(ref, "ref is null");
        Utils.validateArg(ref.length > 0, "ref is empty");
        Utils.nonNull(refLoc, "refLoc is null");
        Utils.validateArg(refLoc.size() == ref.length, "refLoc length must match ref bases");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow is null");

        final List<Haplotype> haplotypes = readLikelihoods.alleles();

        Utils.validateArg(readLikelihoods.samples().contains(tumorSampleName), "readLikelihoods does not contain the tumor sample ");
        final boolean hasNormal = matchedNormalSampleName != null;

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final List<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, ref, refLoc, NO_GIVEN_ALLELES).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = getVCsAtThisLocation(haplotypes, loc, NO_GIVEN_ALLELES);

            if( eventsAtThisLoc.isEmpty() ) { continue; }

            // Create the event mapping object which maps the original haplotype events to the events present at just this locus
            final Map<Event, List<Haplotype>> eventMapper = createEventMapper(loc, eventsAtThisLoc, haplotypes);

            // TODO: priorityList is not sorted by priority, might as well just use eventsAtThisLoc.map(VariantContext::getSource)
            final List<String> priorityList = makePriorityList(eventsAtThisLoc);

            // merge variant contexts from multiple haplotypes into one variant context
            // TODO: we should use haplotypes if possible, but that may have to wait for GATK4
            VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(eventsAtThisLoc, priorityList,
                    GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                    GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

            if( mergedVC == null ) { continue; }

            // TODO: this varaible needs a descriptive name
            final Map<VariantContext, Allele> mergeMap = new LinkedHashMap<>();

            mergeMap.put(null, mergedVC.getReference()); // the reference event (null) --> the reference allele
            for(int i = 0; i < eventsAtThisLoc.size(); i++) {
                // TODO: as noted below, this operation seems dangerous. Understand how things can go wrong.
                mergeMap.put(eventsAtThisLoc.get(i), mergedVC.getAlternateAllele(i)); // BUGBUG: This is assuming that the order of alleles is the same as the priority list given to simpleMerge function
            }

            /** TODO: the code in the for loop up to here needs refactor. The goal, as far as I can tell, is to create two things: alleleMapper and mergedVC
             * alleleMapper maps alleles to haplotypes, and we need this to create readAlleleLikelihoods.
             * To make alleleMapper we make mergeMap (of type VC -> Allele) and eventMapper (of type Event -> List(Haplotypes), where Event is essentialy Variant Context)
             * If we just want a map of Alleles to Haplotypes, we should be able to do so directly; no need for intermediate maps, which just complicates the code.
             **/

            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergeMap, eventMapper);

            // converting ReadLikelihoods<Haplotype> to ReadLikeliHoods<Allele>
            ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper,
                    new SimpleInterval(mergedVC).expandWithinContig(ALLELE_EXTENSION, header.getSequenceDictionary()));
            //LDG: do we want to do this before or after pulling out overlapping reads?
            if (MTAC.isSampleContaminationPresent()) {
                readAlleleLikelihoods.contaminationDownsampling(MTAC.getSampleContamination());
            }

            // TODO: this is a good break point for a new method
            // TODO: replace PRALM with ReadLikelihoods
            final PerReadAlleleLikelihoodMap tumorPRALM = readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(readAlleleLikelihoods.indexOfSample(tumorSampleName));
            filterPRALMForOverlappingReads(tumorPRALM, mergedVC.getReference(), loc, false);

            //TODO: uncomment after porting Mutect2.java
            //MuTect2.logReadInfo(DEBUG_READ_NAME, tumorPRALM.getLikelihoodReadMap().keySet(), "Present in Tumor PRALM after filtering for overlapping reads");
            // extend to multiple samples

            // compute tumor LOD for each alternate allele
            // TODO: somewhere we have to ensure that the all the alleles in the variant context is in alleleFractions passed to getHetGenotypeLogLikelihoods. getHetGenotypeLogLikelihoods will not check that for you
            final PerAlleleCollection<Double> altAlleleFractions = estimateAlleleFraction(mergedVC, tumorPRALM, false);
            final PerAlleleCollection<Double> tumorHetGenotypeLLs = getHetGenotypeLogLikelihoods(mergedVC, tumorPRALM, originalNormalReadQualities, altAlleleFractions);

            final PerAlleleCollection<Double> tumorLods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
            for (final Allele altAllele : mergedVC.getAlternateAlleles()){
                tumorLods.set(altAllele, tumorHetGenotypeLLs.get(altAllele) - tumorHetGenotypeLLs.getRef());
            }

            // TODO: another good breakpoint e.g. compute normal LOD/set thresholds
            // TODO: anything related to normal should be encapsulated in Optional

            // A variant candidate whose normal LOD is below this threshold will be filtered as 'germline_risk'
            // This is a more stringent threshold than normalLodThresholdForVCF
            double normalLodFilterThreshold = -Double.MAX_VALUE;
            PerReadAlleleLikelihoodMap normalPRALM = null;
            final PerAlleleCollection<Double> normalLods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);

            // if normal bam is available, compute normal LOD
            // TODO: this if statement should be a standalone method for computing normal LOD
            // TODO: then we can do something like normalLodThreshold = hasNormal ? thisMethod() : Optional.empty()
            if (hasNormal) {
                normalPRALM = readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(readAlleleLikelihoods.indexOfSample(matchedNormalSampleName));
                filterPRALMForOverlappingReads(normalPRALM, mergedVC.getReference(), loc, true);

                //TODO: uncomment after porting Mutect2.java
                //MuTect2.logReadInfo(DEBUG_READ_NAME, normalPRALM.getLikelihoodReadMap().keySet(), "Present after in Nomral PRALM filtering for overlapping reads");

                final Collection<VariantContext> cosmicVC = featureContext.getValues(MTAC.cosmicFeatureInput, loc);
                final Collection<VariantContext> dbsnpVC = featureContext.getValues(MTAC.dbsnp.dbsnp, loc);
                final boolean germlineAtRisk = !dbsnpVC.isEmpty() && cosmicVC.isEmpty();

                normalLodFilterThreshold = germlineAtRisk ? MTAC.NORMAL_DBSNP_LOD_THRESHOLD : MTAC.NORMAL_LOD_THRESHOLD;

                // compute normal LOD = LL(X|REF)/LL(X|ALT) where REF is the diploid HET with AF = 0.5
                // note normal LOD is REF over ALT, the reciprocal of the tumor LOD
                final PerAlleleCollection<Double> diploidHetAlleleFractions = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
                for (final Allele allele : mergedVC.getAlternateAlleles()){
                    diploidHetAlleleFractions.setAlt(allele, 0.5);
                }

                final PerAlleleCollection<Double> normalGenotypeLLs = getHetGenotypeLogLikelihoods(mergedVC, normalPRALM, originalNormalReadQualities, diploidHetAlleleFractions);

                for (final Allele altAllele : mergedVC.getAlternateAlleles()){
                    normalLods.setAlt(altAllele, normalGenotypeLLs.getRef() - normalGenotypeLLs.getAlt(altAllele));
                }
            }

            int numPassingAlts = 0;
            final Set<Allele> allelesThatPassThreshold = new HashSet<>();
            Allele alleleWithHighestTumorLOD = null;

            for (final Allele altAllele : mergedVC.getAlternateAlleles()) {
                final boolean passesTumorLodThreshold = tumorLods.getAlt(altAllele) >= MTAC.INITIAL_TUMOR_LOD_THRESHOLD;
                final boolean passesNormalLodThreshold = hasNormal ? normalLods.getAlt(altAllele) >= MTAC.INITIAL_NORMAL_LOD_THRESHOLD : true;
                if (passesTumorLodThreshold && passesNormalLodThreshold) {
                    numPassingAlts++;
                    allelesThatPassThreshold.add(altAllele);
                    if (alleleWithHighestTumorLOD == null || tumorLods.getAlt(altAllele) > tumorLods.getAlt(alleleWithHighestTumorLOD)){
                        alleleWithHighestTumorLOD = altAllele;
                    }
                }
            }

            if (numPassingAlts == 0) {
                continue;
            }

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC);
            final int haplotypeCount = alleleMapper.get(alleleWithHighestTumorLOD).size();
            callVcb.attribute(GATKVCFConstants.HAPLOTYPE_COUNT_KEY, haplotypeCount);
            callVcb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, tumorLods.getAlt(alleleWithHighestTumorLOD));

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY, normalLods.getAlt(alleleWithHighestTumorLOD));
                if (normalLods.getAlt(alleleWithHighestTumorLOD) < normalLodFilterThreshold) {
                    callVcb.filter(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
                }
            }

            // TODO: this should be a separate method
            // TODO: move code to MuTect2::calculateFilters()
            if (MTAC.ENABLE_STRAND_ARTIFACT_FILTER && numPassingAlts == 1) {
                final PerReadAlleleLikelihoodMap forwardPRALM = new PerReadAlleleLikelihoodMap();
                final PerReadAlleleLikelihoodMap reversePRALM = new PerReadAlleleLikelihoodMap();
                splitPRALMintoForwardAndReverseReads(tumorPRALM, forwardPRALM, reversePRALM);

                //TODO: uncomment after porting Mutect2.java
                //MuTect2.logReadInfo(DEBUG_READ_NAME, tumorPRALM.getLikelihoodReadMap().keySet(), "Present in tumor PRALM after PRALM is split");
                //MuTect2.logReadInfo(DEBUG_READ_NAME, forwardPRALM.getLikelihoodReadMap().keySet(), "Present in forward PRALM after PRALM is split");
                //MuTect2.logReadInfo(DEBUG_READ_NAME, reversePRALM.getLikelihoodReadMap().keySet(), "Present in reverse PRALM after PRALM is split");

                // TODO: build a new type for probability, likelihood, and log_likelihood. e.g. f_fwd :: probability[], tumorGLs_fwd :: likelihood[]
                // TODO: don't want to call getHetGenotypeLogLikelihoods on more than one alternate alelle. May need to overload it to take a scalar f_fwd.
                final PerAlleleCollection<Double> alleleFractionsForward = estimateAlleleFraction(mergedVC, forwardPRALM, true);
                final PerAlleleCollection<Double> tumorGenotypeLLForward = getHetGenotypeLogLikelihoods(mergedVC, forwardPRALM, originalNormalReadQualities, alleleFractionsForward);

                final PerAlleleCollection<Double> alleleFractionsReverse = estimateAlleleFraction(mergedVC, reversePRALM, true);
                final PerAlleleCollection<Double> tumorGenotypeLLReverse = getHetGenotypeLogLikelihoods(mergedVC, reversePRALM, originalNormalReadQualities, alleleFractionsReverse);

                final double tumorLod_fwd = tumorGenotypeLLForward.getAlt(alleleWithHighestTumorLOD) - tumorGenotypeLLForward.getRef();
                final double tumorLod_rev = tumorGenotypeLLReverse.getAlt(alleleWithHighestTumorLOD) - tumorGenotypeLLReverse.getRef();

                // Note that we use the observed combined (+ and -) allele fraction for power calculation in either direction
                final double tumorSBpower_fwd = strandArtifactPowerCalculator.cachedPowerCalculation(forwardPRALM.size(), altAlleleFractions.getAlt(alleleWithHighestTumorLOD));
                final double tumorSBpower_rev = strandArtifactPowerCalculator.cachedPowerCalculation(reversePRALM.size(), altAlleleFractions.getAlt(alleleWithHighestTumorLOD));

                callVcb.attribute(MuTect2.TLOD_FWD_KEY, tumorLod_fwd);
                callVcb.attribute(MuTect2.TLOD_REV_KEY, tumorLod_rev);
                callVcb.attribute(MuTect2.TUMOR_SB_POWER_FWD_KEY, tumorSBpower_fwd);
                callVcb.attribute(MuTect2.TUMOR_SB_POWER_REV_KEY, tumorSBpower_rev);

                if ((tumorSBpower_fwd > MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && tumorLod_fwd < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                        (tumorSBpower_rev > MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && tumorLod_rev < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD))
                    callVcb.filter(MuTect2.STRAND_ARTIFACT_FILTER_NAME);
            }

            // TODO: this probably belongs in M2::calculateFilters()
            if (numPassingAlts > 1) {
                callVcb.filter(MuTect2.TRIALLELIC_SITE_FILTER_NAME);
            }

            // build genotypes TODO: this part needs review and refactor
            final List<Allele> tumorAlleles = Arrays.asList(mergedVC.getReference(), alleleWithHighestTumorLOD);
            // TODO: estimateAlleleFraction should not repeat counting allele depths
            final PerAlleleCollection<Integer> tumorAlleleDepths = getRefAltCount(mergedVC, tumorPRALM, false);
            final int tumorRefAlleleDepth = tumorAlleleDepths.getRef();
            final int tumorAltAlleleDepth = tumorAlleleDepths.getAlt(alleleWithHighestTumorLOD);
            final Genotype tumorGenotype = new GenotypeBuilder(tumorSampleName, tumorAlleles)
                    .AD(new int[] { tumorRefAlleleDepth, tumorAltAlleleDepth })
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, altAlleleFractions.getAlt(alleleWithHighestTumorLOD))
                    .make();

            final List<Genotype> genotypes = new ArrayList<>();
            genotypes.add(tumorGenotype);

            // We assume that the genotype in the normal is 0/0
            // TODO: is normal always homozygous reference?
            final List<Allele> homRefAllelesforNormalGenotype = Collections.nCopies(2, mergedVC.getReference());

            // if we are calling with a normal, build the genotype for the sample to appear in vcf
            if (hasNormal) {
                final PerAlleleCollection<Integer> normalAlleleDepths = getRefAltCount(mergedVC, normalPRALM, false);
                final int normalRefAlleleDepth = normalAlleleDepths.getRef();
                final int normalAltAlleleDepth = normalAlleleDepths.getAlt(alleleWithHighestTumorLOD);
                final double normalAlleleFraction = (double) normalAltAlleleDepth / ( normalRefAlleleDepth + normalAltAlleleDepth);

                final Genotype normalGenotype = new GenotypeBuilder(matchedNormalSampleName, homRefAllelesforNormalGenotype)
                        .AD(new int[] { normalRefAlleleDepth, normalAltAlleleDepth })
                        .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, normalAlleleFraction)
                        .make();
                genotypes.add(normalGenotype);
            }

            final VariantContext call = new VariantContextBuilder(callVcb).alleles(tumorAlleles).genotypes(genotypes).make();

            // how should we be making use of _perSampleFilteredReadList_?
            readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
                    false, alleleMapper, readAlleleLikelihoods, call);

            final SimpleInterval locus = new SimpleInterval(mergedVC.getContig(), mergedVC.getStart(), mergedVC.getEnd());
            final SimpleInterval refLocInterval= new SimpleInterval(refLoc);
            final ReferenceDataSource refData = new ReferenceMemorySource(new ReferenceBases(ref, refLocInterval), header.getSequenceDictionary());
            final ReferenceContext referenceContext = new ReferenceContext(refData, locus, refLocInterval);

            VariantContext annotatedCall = annotationEngine.annotateContext(call, featureContext, referenceContext, readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(), a -> true);

            if( call.getAlleles().size() != mergedVC.getAlleles().size() )
                annotatedCall = GATKVariantContextUtils.reverseTrimAlleles(annotatedCall);

            // maintain the set of all called haplotypes
            call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            returnCalls.add( annotatedCall );
        }

        // TODO: understand effect of enabling this for somatic calling...
        final List<VariantContext> outputCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        return new CalledHaplotypes(outputCalls, calledHaplotypes);
    }

    /** Calculate the likelihoods of hom ref and each het genotype of the form ref/alt
     *
     * @param mergedVC                              input VC
     * @param tumorPRALM                            read likelihoods
     * @param originalNormalMQs                     original MQs, before boosting normals to avoid qual capping
     * @param alleleFractions                       allele fraction(s) for alternate allele(s)
     *
     * @return                                      genotype likelihoods for homRef and het for each alternate allele
     */
    private PerAlleleCollection<Double> getHetGenotypeLogLikelihoods(final VariantContext mergedVC,
                                                                     final PerReadAlleleLikelihoodMap tumorPRALM,
                                                                     final Map<String, Integer> originalNormalMQs,
                                                                     final PerAlleleCollection<Double> alleleFractions) {
        Utils.validateArg(mergedVC.getAlternateAlleles().containsAll(alleleFractions.getAltAlleles()), "alleleFractions has alleles that are not in the variant context");

        final PerAlleleCollection<MutableDouble> genotypeLogLikelihoods = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        mergedVC.getAlleles().forEach(a -> genotypeLogLikelihoods.set(a, new MutableDouble(0)));

        final Allele refAllele = mergedVC.getReference();
        for(Map.Entry<GATKRead,Map<Allele, Double>> readAlleleLikelihoodMap : tumorPRALM.getLikelihoodReadMap().entrySet()) {
            final Map<Allele, Double> alleleLikelihoodMap = readAlleleLikelihoodMap.getValue();
            if (originalNormalMQs.get(readAlleleLikelihoodMap.getKey().getName()) == 0) {
                continue;
            }

            final double readRefLogLikelihood = alleleLikelihoodMap.get(refAllele);
            genotypeLogLikelihoods.getRef().add(readRefLogLikelihood);

            for (final Allele altAllele : alleleFractions.getAltAlleles()) {
                final double readAltLogLikelihood = alleleLikelihoodMap.get(altAllele);
                final double adjustedReadAltLL = Math.log10(
                        Math.pow(10, readRefLogLikelihood) * (1 - alleleFractions.getAlt(altAllele)) +
                                Math.pow(10, readAltLogLikelihood) * alleleFractions.getAlt(altAllele)
                );
                genotypeLogLikelihoods.get(altAllele).add(adjustedReadAltLL);
            }
        }

        final PerAlleleCollection<Double> result = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        mergedVC.getAlleles().stream().forEach(a -> result.set(a,genotypeLogLikelihoods.get(a).toDouble()));

        return result;
    }

    /**
     * Find the allele fractions for each alternate allele
     *
     * @param vc                        input VC, for alleles
     * @param pralm                     read likelihoods
     * @return                          estimated AF for each alt
     */
    // FIXME: calculate using the uncertainty rather than this cheap approach
    private PerAlleleCollection<Double> estimateAlleleFraction(final VariantContext vc,
                                                               final PerReadAlleleLikelihoodMap pralm,
                                                               final boolean oneStrandOnly) {
        final PerAlleleCollection<Integer> alleleCounts = getRefAltCount(vc, pralm, oneStrandOnly);
        final PerAlleleCollection<Double> alleleFractions = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);

        final int refCount = alleleCounts.getRef();
        for ( final Allele altAllele : vc.getAlternateAlleles() ) {
            final int altCount = alleleCounts.getAlt(altAllele);
            double alleleFraction = (double) altCount / (refCount + altCount);
            // weird case, but I've seen it happen in one strand cases
            if (refCount == 0 && altCount == refCount ) {
                alleleFraction = 0;
            }
            alleleFractions.setAlt(altAllele, alleleFraction);
            // logger.info("Counted " + refCount + " ref and " + altCount + " alt " );
        }

        return alleleFractions;
    }

    /**
     *  Go through the PRALM and tally the most likely allele in each read. Only count informative reads.
     *
     * @param vc                      input VC, for alleles
     * @param pralm                         read likelihoods
     * @return                              an array giving the read counts for the ref and each alt allele
     */
    private PerAlleleCollection<Integer> getRefAltCount(final VariantContext vc,
                                                        final PerReadAlleleLikelihoodMap pralm,
                                                        final boolean oneStrandOnly) {
        // Check that the alleles in Variant Context are in PRALM
        // Skip the check for strand-conscious PRALM; + reads may not have alleles in - reads, for example.
        final Set<Allele> vcAlleles = new HashSet<>(vc.getAlleles());
        if ( ! oneStrandOnly && ! pralm.getAllelesSet().containsAll( vcAlleles ) ) {
            StringBuilder message = new StringBuilder();
            message.append("At Locus chr" + vc.getContig() + ":" + vc.getStart() + ", we detected that variant context had alleles that not in PRALM. ");
            message.append("VC alleles = " + vcAlleles + ", PRALM alleles = " + pralm.getAllelesSet());
            logger.warn(message);
        }


        final PerAlleleCollection<MutableInt> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        vcAlleles.stream().forEach(a -> alleleCounts.set(a, new MutableInt(0)));

        for (final Map.Entry<GATKRead, Map<Allele, Double>> readAlleleLikelihoodMap : pralm.getLikelihoodReadMap().entrySet()) {
            final GATKRead read = readAlleleLikelihoodMap.getKey();
            final Map<Allele, Double> alleleLikelihoodMap = readAlleleLikelihoodMap.getValue();
            final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(alleleLikelihoodMap, vcAlleles);

            if (read.getMappingQuality() > 0 && mostLikelyAllele.isInformative()) {
                alleleCounts.get(mostLikelyAllele.getMostLikelyAllele()).increment();
            }

        }

        final PerAlleleCollection<Integer> result = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        vc.getAlleles().stream().forEach(a -> result.set(a, alleleCounts.get(a).toInteger()));

        return(result);
    }

    private void logM2Debug(String s) {
        if (MTAC.M2_DEBUG) {
            logger.info(s);
        }
    }

    private void filterPRALMForOverlappingReads(final PerReadAlleleLikelihoodMap pralm, final Allele ref, final int location, final boolean retainMismatches) {
        final Map<GATKRead, Map<Allele, Double>> m = pralm.getLikelihoodReadMap();

        // iterate through the reads, if the name has been seen before we have overlapping (potentially) fragments, so handle them
        final Map<String, GATKRead> nameToRead = new HashMap<>();
        final Set<GATKRead> readsToKeep = new HashSet<>();

        for(final GATKRead rec : m.keySet()) {
            // if we haven't seen it... just record the name and add it to the list of reads to keep
            final GATKRead existing = nameToRead.get(rec.getName());
            if (existing == null) {
                nameToRead.put(rec.getName(), rec);
                readsToKeep.add(rec);
            } else {
                logM2Debug("Found a paired read for " + rec.getName());

                // NOTE: Can we use FragmentUtils to do all of this processing (to find overlapping pairs?)
                // seems like maybe, but it has some requirements about the order of the reads supplied which may be painful to meet
                // TODO: CHECK IF THE READS BOTH OVERLAP THE POSITION!!!!
                if ( ReadUtils.isInsideRead(existing, location) && ReadUtils.isInsideRead(rec, location) ) {

                    final MostLikelyAllele existingMLA = PerReadAlleleLikelihoodMap.getMostLikelyAllele(pralm.getLikelihoodReadMap().get(existing));
                    final Allele existingAllele = existingMLA.getMostLikelyAllele();

                    final MostLikelyAllele recMLA = PerReadAlleleLikelihoodMap.getMostLikelyAllele(pralm.getLikelihoodReadMap().get(rec));
                    final Allele recAllele = recMLA.getMostLikelyAllele();

                    // if the reads disagree at this position...
                    if (!existingAllele.equals(recAllele)) {
                        //... and we're not retaining mismatches, throw them both out
                        if (!retainMismatches) {
                            logM2Debug("Discarding read-pair due to disagreement" + rec.getName() + " and allele " + existingAllele);
                            readsToKeep.remove(existing);

                            //... and we are retaining mismatches, keep the mismatching one
                        } else {
                            if (existingAllele.equals(ref)) {
                                logM2Debug("Discarding read to keep mismatching " + rec.getName() + " and allele " + existingAllele);
                                readsToKeep.remove(existing);
                                readsToKeep.add(rec);
                            }
                        }
                        // Otherwise, keep the element with the higher quality score
                    } else {
                        logM2Debug("Discarding lower quality read of overlapping pair " + rec.getName() + " and allele " + existingAllele);
                        if (existingMLA.getLog10LikelihoodOfMostLikely() < recMLA.getLog10LikelihoodOfMostLikely()) {
                            readsToKeep.remove(existing);
                            readsToKeep.add(rec);
                        }
                    }
                } else {
                    // although these are overlapping fragments, they don't overlap at the position in question
                    // so keep the read
                    readsToKeep.add(rec);
                }
            }

        }

        // perhaps moved into PRALM
        final Iterator<Map.Entry<GATKRead, Map<Allele, Double>>> it = m.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<GATKRead, Map<Allele, Double>> record = it.next();
            if(!readsToKeep.contains(record.getKey())) {
                it.remove();
                logM2Debug("Dropping read " + record.getKey() + " due to overlapping read fragment rules");
            }
        }
    }

    private void splitPRALMintoForwardAndReverseReads(final PerReadAlleleLikelihoodMap originalPRALM, final PerReadAlleleLikelihoodMap forwardPRALM, final PerReadAlleleLikelihoodMap reversePRALM) {
        final Map<GATKRead, Map<Allele, Double>> origReadAlleleLikelihoodMap = originalPRALM.getLikelihoodReadMap();
        for (final GATKRead read : origReadAlleleLikelihoodMap.keySet()) {

            //TODO: does GATK4 have strandless reads?
            //if (read.isStrandless())
            //    continue;

            for (final Map.Entry<Allele, Double> alleleLikelihoodMap : origReadAlleleLikelihoodMap.get(read).entrySet()) {
                final Allele allele = alleleLikelihoodMap.getKey();
                final Double likelihood = alleleLikelihoodMap.getValue();
                if (read.isReverseStrand())
                    reversePRALM.add(read, allele, likelihood);
                else
                    forwardPRALM.add(read, allele, likelihood);
            }
        }
    }
}
