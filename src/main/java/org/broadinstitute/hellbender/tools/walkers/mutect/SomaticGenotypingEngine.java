package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticGenotypingEngine {

    private final M2ArgumentCollection MTAC;
    private final Set<String> normalSamples;
    final boolean hasNormal;
    protected VariantAnnotatorEngine annotationEngine;

    // If MTAC.minAF is non-zero we softly cut off allele fractions below minAF with a Beta prior of the form Beta(1+epsilon, 1); that is
    // the prior on allele fraction f is proportional to f^epsilon.  If epsilon is small this prior vanishes as f -> 0
    // and very rapidly becomes flat.  We choose epsilon such that minAF^epsilon = 0.5.
    private final double refPseudocount = 1;
    private final double altPseudocount;

    public SomaticGenotypingEngine(final M2ArgumentCollection MTAC, final Set<String> normalSamples, final VariantAnnotatorEngine annotationEngine) {
        this.MTAC = MTAC;
        altPseudocount = MTAC.minAF == 0.0 ? 1 : 1 - Math.log(2)/Math.log(MTAC.minAF);

        this.normalSamples = normalSamples;
        hasNormal = !normalSamples.isEmpty();
        this.annotationEngine = annotationEngine;
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param logReadLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param activeRegionWindow                     Active window
     * @param withBamOut                            whether to annotate reads in readLikelihoods for future writing to bamout
     * @param emitRefConf                           generate reference confidence (GVCF) data?
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     */
    public CalledHaplotypes callMutations(
            final AlleleLikelihoods<GATKRead, Haplotype> logReadLikelihoods,
            final AssemblyResultSet assemblyResultSet,
            final ReferenceContext referenceContext,
            final SimpleInterval activeRegionWindow,
            final FeatureContext featureContext,
            final List<VariantContext> givenAlleles,
            final SAMFileHeader header,
            final boolean withBamOut,
            final boolean emitRefConf) {
        Utils.nonNull(logReadLikelihoods);
        Utils.validateArg(logReadLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(activeRegionWindow);

        final List<Haplotype> haplotypes = logReadLikelihoods.alleles();

        final List<Integer> startPosKeySet = EventMap.buildEventMapsForHaplotypes(haplotypes, assemblyResultSet.getFullReferenceWithPadding(),
                assemblyResultSet.getPaddedReferenceLoc(), MTAC.assemblerArgs.debugAssembly, MTAC.maxMnpDistance).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        if(withBamOut){
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(logReadLikelihoods, activeRegionWindow);
        }

        if (MTAC.likelihoodArgs.phredScaledGlobalReadMismappingRate > 0) {
            logReadLikelihoods.normalizeLikelihoods(NaturalLogUtils.qualToLogErrorProb(MTAC.likelihoodArgs.phredScaledGlobalReadMismappingRate), true);
        }
        final AlleleLikelihoods<Fragment, Haplotype> logFragmentLikelihoods = logReadLikelihoods.groupEvidence(MTAC.independentMates ? read -> read : GATKRead::getName, Fragment::createAndAvoidFailure);

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(loc, haplotypes, false);
            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( mergedVC == null ) {
                continue;
            }

            // converting haplotype likelihoods to allele likelihoods
            final Map<Allele, List<Haplotype>> alleleMapper = AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, haplotypes, true);
            final AlleleLikelihoods<Fragment, Allele> logLikelihoods = logFragmentLikelihoods.marginalize(alleleMapper);
            final SimpleInterval variantCallingRelevantFragmentOverlap = new SimpleInterval(mergedVC).expandWithinContig(MTAC.informativeReadOverlapMargin, header.getSequenceDictionary());
            logLikelihoods.retainEvidence(variantCallingRelevantFragmentOverlap::overlaps);

            if (emitRefConf) {
                mergedVC = ReferenceConfidenceUtils.addNonRefSymbolicAllele(mergedVC);
                logLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }
            final List<LikelihoodMatrix<Fragment, Allele>> tumorMatrices = IntStream.range(0, logLikelihoods.numberOfSamples())
                    .filter(n -> !normalSamples.contains(logLikelihoods.getSample(n)))
                    .mapToObj(logLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final AlleleList<Allele> alleleList = tumorMatrices.get(0);
            final LikelihoodMatrix<Fragment, Allele> logTumorMatrix = combinedLikelihoodMatrix(tumorMatrices, alleleList);
            final PerAlleleCollection<Double> tumorLogOdds = somaticLogOdds(logTumorMatrix);

            final List<LikelihoodMatrix<Fragment, Allele>> normalMatrices = IntStream.range(0, logLikelihoods.numberOfSamples())
                    .filter(n -> normalSamples.contains(logLikelihoods.getSample(n)))
                    .mapToObj(logLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final LikelihoodMatrix<Fragment, Allele> logNormalMatrix = combinedLikelihoodMatrix(normalMatrices, alleleList);
            final PerAlleleCollection<Double> normalLogOdds = diploidAltLogOdds(logNormalMatrix);
            final PerAlleleCollection<Double> normalArtifactLogOdds = somaticLogOdds(logNormalMatrix);


            final Set<Allele> forcedAlleles = AssemblyBasedCallerUtils.getAllelesConsistentWithGivenAlleles(givenAlleles, mergedVC);
            final List<Allele> tumorAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> forcedAlleles.contains(allele) || tumorLogOdds.getAlt(allele) > MTAC.getEmissionLogOdds())
                    .collect(Collectors.toList());

            final long somaticAltCount = tumorAltAlleles.stream()
                    .filter(allele -> forcedAlleles.contains(allele) || !hasNormal || MTAC.genotypeGermlineSites || normalLogOdds.getAlt(allele) > MathUtils.log10ToLog(MTAC.normalLog10Odds))
                    .count();

            // if every alt allele is germline, skip this variant.  However, if some alt alleles are germline and others
            // are not we emit them all so that the filtering engine can see them
            if (somaticAltCount == 0) {
                continue;
            }

            final List<Allele> allAllelesToEmit = ListUtils.union(Arrays.asList(mergedVC.getReference()), tumorAltAlleles);


            final Map<String, Object> negativeLogPopulationAFAnnotation = getNegativeLogPopulationAFAnnotation(featureContext.getValues(MTAC.germlineResource, loc), tumorAltAlleles, MTAC.getDefaultAlleleFrequency());

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC)
                    .alleles(allAllelesToEmit)
                    .attributes(negativeLogPopulationAFAnnotation)
                    .attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, tumorAltAlleles.stream().mapToDouble(a -> MathUtils.logToLog10(tumorLogOdds.getAlt(a))).toArray());

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_ARTIFACT_LOG_10_ODDS_KEY,
                        Arrays.stream(normalArtifactLogOdds.asDoubleArray(tumorAltAlleles)).map(x->-MathUtils.logToLog10(x)).toArray());
                callVcb.attribute(GATKVCFConstants.NORMAL_LOG_10_ODDS_KEY,
                        Arrays.stream(normalLogOdds.asDoubleArray(tumorAltAlleles)).map(MathUtils::logToLog10).toArray());
            }

            if (!featureContext.getValues(MTAC.pon, mergedVC.getStart()).isEmpty()) {
                callVcb.attribute(GATKVCFConstants.IN_PON_KEY, true);
            }

            addGenotypes(logLikelihoods, allAllelesToEmit, callVcb);
            final VariantContext call = callVcb.make();
            final VariantContext trimmedCall = GATKVariantContextUtils.trimAlleles(call, true, true);
            final List<Allele> trimmedAlleles = trimmedCall.getAlleles();
            final List<Allele> untrimmedAlleles = call.getAlleles();
            final Map<Allele, List<Allele>> trimmedToUntrimmedAlleleMap = IntStream.range(0, trimmedCall.getNAlleles()).boxed()
                    .collect(Collectors.toMap(n -> trimmedAlleles.get(n), n -> Arrays.asList(untrimmedAlleles.get(n))));
            final AlleleLikelihoods<Fragment, Allele> trimmedLikelihoods = logLikelihoods.marginalize(trimmedToUntrimmedAlleleMap);

            // AlleleLikelihoods for annotation only
            final AlleleLikelihoods<GATKRead, Allele> logReadAlleleLikelihoods = logReadLikelihoods.marginalize(alleleMapper);
            logReadAlleleLikelihoods.retainEvidence(variantCallingRelevantFragmentOverlap::overlaps);

            if (emitRefConf) {
                logReadAlleleLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }

            final AlleleLikelihoods<GATKRead, Allele> trimmedLikelihoodsForAnnotation = logReadAlleleLikelihoods.marginalize(trimmedToUntrimmedAlleleMap);


            final VariantContext annotatedCall =  annotationEngine.annotateContext(trimmedCall, featureContext, referenceContext, trimmedLikelihoodsForAnnotation, a -> true);
            if(withBamOut) {
                AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(trimmedCall, trimmedLikelihoods, Fragment::getReads);
            }

            call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            returnCalls.add( annotatedCall );
        }

        final List<VariantContext> outputCalls = AssemblyBasedCallerUtils.phaseCalls(returnCalls, calledHaplotypes);
        final int eventCount = outputCalls.size();
        final List<VariantContext> outputCallsWithEventCountAnnotation = outputCalls.stream()
                .map(vc -> new VariantContextBuilder(vc).attribute(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount).make())
                .collect(Collectors.toList());
        return new CalledHaplotypes(outputCallsWithEventCountAnnotation, calledHaplotypes);
    }

    private double[] makePriorPseudocounts(final int numAlleles) {
        return new IndexRange(0, numAlleles).mapToDouble(n -> n == 0 ? refPseudocount : altPseudocount);
    }

    // compute the likelihoods that each allele is contained at some allele fraction in the sample
    protected <EVIDENCE extends Locatable> PerAlleleCollection<Double> somaticLogOdds(final LikelihoodMatrix<EVIDENCE, Allele> logMatrix) {
        final int alleleListEnd = logMatrix.alleles().size()-1;
        if (logMatrix.alleles().contains(Allele.NON_REF_ALLELE) && !(logMatrix.alleles().get(alleleListEnd).equals(Allele.NON_REF_ALLELE))) {
            throw new IllegalStateException("<NON_REF> must be last in the allele list.");
        }

        final double logEvidenceWithAllAlleles = logMatrix.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(getAsRealMatrix(logMatrix), makePriorPseudocounts(logMatrix.numberOfAlleles()));

        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int refIndex = getRefIndex(logMatrix);
        IntStream.range(0, logMatrix.numberOfAlleles()).filter(a -> a != refIndex).forEach( a -> {
            final Allele allele = logMatrix.getAllele(a);
            final LikelihoodMatrix<EVIDENCE, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(logMatrix, allele);
            final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.evidenceCount() == 0 ? 0 :
                    SomaticLikelihoodsEngine.logEvidence(getAsRealMatrix(logMatrixWithoutThisAllele), makePriorPseudocounts(logMatrixWithoutThisAllele.numberOfAlleles()));
            lods.setAlt(allele, logEvidenceWithAllAlleles - logEvidenceWithoutThisAllele);
        });
        return lods;
    }

    private <EVIDENCE extends Locatable> void addGenotypes(final AlleleLikelihoods<EVIDENCE, Allele> logLikelihoods,
                              final List<Allele> allelesToEmit,
                              final VariantContextBuilder callVcb) {
        final List<Genotype> genotypes = IntStream.range(0, logLikelihoods.numberOfSamples()).mapToObj(n -> {
            final String sample = logLikelihoods.getSample(n);
            final LikelihoodMatrix<EVIDENCE, Allele> logMatrix = new SubsettedLikelihoodMatrix<>(logLikelihoods.sampleMatrix(n), allelesToEmit);
            final double[] alleleCounts = getEffectiveCounts(logMatrix);
            final double[] flatPriorPseudocounts = new IndexRange(0, logMatrix.numberOfAlleles()).mapToDouble(a -> 1);
            final double[] alleleFractionsPosterior = logMatrix.evidenceCount() == 0 ? flatPriorPseudocounts :
                    SomaticLikelihoodsEngine.alleleFractionsPosterior(getAsRealMatrix(logMatrix), flatPriorPseudocounts);
            final double[] tumorAlleleFractionsMean = MathUtils.normalizeSumToOne(alleleFractionsPosterior);

            // TODO: We shouldn't always assume that the genotype in the normal is hom ref
            final Allele ref = logMatrix.getAllele(getRefIndex(logMatrix));
            return new GenotypeBuilder(sample, normalSamples.contains(sample) ? Collections.nCopies(2, ref) : logMatrix.alleles())
                    .AD(Arrays.stream(alleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, Arrays.copyOfRange(tumorAlleleFractionsMean, 1, tumorAlleleFractionsMean.length))
                    .make();
        }).collect(Collectors.toList());

        callVcb.genotypes(genotypes);
    }

    private static <EVIDENCE> double[] getEffectiveCounts(final LikelihoodMatrix<EVIDENCE, Allele> logLikelihoodMatrix) {
        if (logLikelihoodMatrix.evidenceCount() == 0) {
            return new double[logLikelihoodMatrix.numberOfAlleles()]; // zero counts for each allele
        }
        final RealMatrix logLikelihoods = getAsRealMatrix(logLikelihoodMatrix);
        return MathUtils.sumArrayFunction(0, logLikelihoods.getColumnDimension(),
                read -> NaturalLogUtils.normalizeFromLogToLinearSpace(logLikelihoods.getColumn(read)));
    }

    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * @param matrix a matrix of log likelihoods
     */
    private <EVIDENCE extends Locatable> PerAlleleCollection<Double> diploidAltLogOdds(final LikelihoodMatrix<EVIDENCE, Allele> matrix) {
        final int refIndex = getRefIndex(matrix);
        final int numReads = matrix.evidenceCount();
        final double homRefLogLikelihood = new IndexRange(0, numReads).sum(r -> matrix.get(refIndex,r));

        final PerAlleleCollection<Double> result = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        // hom ref likelihood for the ref allele, het likelihood for alt alleles
        IntStream.range(0, matrix.numberOfAlleles()).filter(a -> a != refIndex)
                .forEach( a -> {
                    final double hetLogLikelihood = new IndexRange(0, numReads)
                            .sum(r -> NaturalLogUtils.logSumExp(matrix.get(refIndex, r), matrix.get(a,r)) + NaturalLogUtils.LOG_ONE_HALF);
                    result.setAlt(matrix.getAllele(a), homRefLogLikelihood - hetLogLikelihood);
                });
        return result;
    }

    private <EVIDENCE> int getRefIndex(LikelihoodMatrix<EVIDENCE, Allele> matrix) {
        final OptionalInt optionalRefIndex = IntStream.range(0, matrix.numberOfAlleles()).filter(a -> matrix.getAllele(a).isReference()).findFirst();
        Utils.validateArg(optionalRefIndex.isPresent(), "No ref allele found in likelihoods");
        return optionalRefIndex.getAsInt();
    }

    //convert a likelihood matrix of alleles x reads into a RealMatrix
    public static <EVIDENCE> RealMatrix getAsRealMatrix(final LikelihoodMatrix<EVIDENCE, Allele> matrix) {
        final RealMatrix result = new Array2DRowRealMatrix(matrix.numberOfAlleles(), matrix.evidenceCount());
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return matrix.get(row, column);
            }
        });
        return result;
    }

    private static <EVIDENCE extends Locatable> LikelihoodMatrix<EVIDENCE, Allele> combinedLikelihoodMatrix(final List<LikelihoodMatrix<EVIDENCE, Allele>> matrices, final AlleleList<Allele> alleleList) {
        final List<EVIDENCE> reads = matrices.stream().flatMap(m -> m.evidence().stream()).collect(Collectors.toList());
        final AlleleLikelihoods<EVIDENCE, Allele> combinedLikelihoods = new AlleleLikelihoods<>(SampleList.singletonSampleList("COMBINED"), alleleList, ImmutableMap.of("COMBINED", reads));

        int combinedReadIndex = 0;
        final LikelihoodMatrix<EVIDENCE, Allele> result = combinedLikelihoods.sampleMatrix(0);
        final int alleleCount = result.numberOfAlleles();
        for (final LikelihoodMatrix<EVIDENCE, Allele> matrix : matrices) {
            final int readCount = matrix.evidenceCount();
            for (int r = 0; r < readCount; r++) {
                for (int a = 0; a < alleleCount; a++) {
                    result.set(a, combinedReadIndex, matrix.get(a, r));
                }
                combinedReadIndex++;
            }
        }
        return result;
    }

    private <E> Optional<E> getForNormal(final Supplier<E> supplier) {
        return hasNormal ? Optional.of(supplier.get()) : Optional.empty();
    }

    private static Map<String, Object> getNegativeLogPopulationAFAnnotation(List<VariantContext> germlineResourceVariants,
                                                                            final List<Allele> altAlleles,
                                                                            final double afOfAllelesNotInGermlineResource) {
        final Optional<VariantContext> germlineVC = germlineResourceVariants.isEmpty() ? Optional.empty()
                : Optional.of(germlineResourceVariants.get(0));  // assume only one VC per site
        final double[] populationAlleleFrequencies = getGermlineAltAlleleFrequencies(altAlleles, germlineVC, afOfAllelesNotInGermlineResource);

        return ImmutableMap.of(GATKVCFConstants.POPULATION_AF_KEY, MathUtils.applyToArray(populationAlleleFrequencies, x -> - Math.log10(x)));
    }

    @VisibleForTesting
    static double[] getGermlineAltAlleleFrequencies(final List<Allele> altAlleles, final Optional<VariantContext> germlineVC, final double afOfAllelesNotInGermlineResource) {
        if (germlineVC.isPresent())  {
            final List<Double> germlineAltAFs = Mutect2Engine.getAttributeAsDoubleList(germlineVC.get(), VCFConstants.ALLELE_FREQUENCY_KEY, afOfAllelesNotInGermlineResource);
            return altAlleles.stream()
                    .mapToDouble(allele -> {
                        final VariantContext vc = germlineVC.get();
                        final OptionalInt germlineAltIndex = IntStream.range(0, vc.getNAlleles() - 1)
                                .filter(n -> vc.getAlternateAllele(n).basesMatch(allele))
                                .findAny();
                        return germlineAltIndex.isPresent() ? germlineAltAFs.get(germlineAltIndex.getAsInt())
                                : afOfAllelesNotInGermlineResource;
                    }).toArray();
        }

        return Doubles.toArray(Collections.nCopies(altAlleles.size(), afOfAllelesNotInGermlineResource));
    }
}
