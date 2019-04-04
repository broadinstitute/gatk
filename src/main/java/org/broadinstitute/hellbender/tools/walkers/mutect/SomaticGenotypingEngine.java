package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
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

    public SomaticGenotypingEngine(final M2ArgumentCollection MTAC, final Set<String> normalSamples, final VariantAnnotatorEngine annotationEngine) {
        this.MTAC = MTAC;
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
            final ReadLikelihoods<Haplotype> logReadLikelihoods,
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

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(loc, haplotypes, false);
            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( mergedVC == null ) {
                continue;
            }

            // converting ReadLikelihoods<Haplotype> to ReadLikelihoods<Allele>
            final Map<Allele, List<Haplotype>> alleleMapper = AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, haplotypes, null);
            final ReadLikelihoods<Allele> logLikelihoods = logReadLikelihoods.marginalize(alleleMapper,
                    new SimpleInterval(mergedVC).expandWithinContig(HaplotypeCallerGenotypingEngine.ALLELE_EXTENSION, header.getSequenceDictionary()));

            if (emitRefConf) {
                mergedVC = ReferenceConfidenceUtils.addNonRefSymbolicAllele(mergedVC);
                logLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }
            final List<LikelihoodMatrix<Allele>> tumorMatrices = IntStream.range(0, logLikelihoods.numberOfSamples())
                    .filter(n -> !normalSamples.contains(logLikelihoods.getSample(n)))
                    .mapToObj(logLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final AlleleList<Allele> alleleList = tumorMatrices.get(0);
            final LikelihoodMatrix<Allele> logTumorMatrix = combinedLikelihoodMatrix(tumorMatrices, alleleList);
            final PerAlleleCollection<Double> tumorLogOdds = somaticLogOdds(logTumorMatrix);

            final List<LikelihoodMatrix<Allele>> normalMatrices = IntStream.range(0, logLikelihoods.numberOfSamples())
                    .filter(n -> normalSamples.contains(logLikelihoods.getSample(n)))
                    .mapToObj(logLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final LikelihoodMatrix<Allele> logNormalMatrix = combinedLikelihoodMatrix(normalMatrices, alleleList);
            final PerAlleleCollection<Double> normalLogOdds = diploidAltLogOdds(logNormalMatrix);
            final PerAlleleCollection<Double> normalArtifactLogOdds = somaticLogOdds(logNormalMatrix);


            final Set<Allele> forcedAlleles = getAllelesConsistentWithGivenAlleles(givenAlleles, loc, mergedVC);
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
            final ReadLikelihoods<Allele> trimmedLikelihoods = logLikelihoods.marginalize(trimmedToUntrimmedAlleleMap);

            final VariantContext annotatedCall =  annotationEngine.annotateContext(trimmedCall, featureContext, referenceContext, trimmedLikelihoods, a -> true);
            if(withBamOut) {
                AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(trimmedCall, trimmedLikelihoods);
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

    private Set<Allele> getAllelesConsistentWithGivenAlleles(List<VariantContext> givenAlleles, int loc, VariantContext mergedVC) {
        final List<Pair<Allele, Allele>> givenAltAndRefAllelesInOriginalContext =  AssemblyBasedCallerUtils.getVariantContextsFromGivenAlleles(loc, givenAlleles, false).stream()
                .flatMap(vc -> vc.getAlternateAlleles().stream().map(allele -> ImmutablePair.of(allele, vc.getReference()))).collect(Collectors.toList());

        return mergedVC.getAlternateAlleles().stream()
                .map(allele -> ImmutablePair.of(allele, mergedVC.getReference()))
                .filter(altAndRef -> givenAltAndRefAllelesInOriginalContext.stream().anyMatch(givenAltAndRef -> allelesAreConsistent(givenAltAndRef, altAndRef)))
                .map(altAndRefPair -> altAndRefPair.getLeft())
                .collect(Collectors.toSet());
    }

    // check whether two alleles coming from different variant contexts and with possibly different reference alleles
    // could in fact be the same.  The condition is that one is a prefix of the other
    private boolean allelesAreConsistent(final Pair<Allele,Allele> altAndRef1, final Pair<Allele,Allele> altAndRef2) {
        final Allele alt1 = altAndRef1.getLeft();
        final Allele alt2 = altAndRef2.getLeft();
        if (alt1.isSymbolic() || alt2.isSymbolic()) {
            return false;
        } else {
            final int sizeDiff1 = alt1.length() - altAndRef1.getRight().length();
            final int sizeDiff2 = alt2.length() - altAndRef2.getRight().length();
            return (sizeDiff1 == sizeDiff2) && (alt1.length() < alt2.length() ?
                    alt1.basesMatch(Arrays.copyOf(alt2.getBases(), alt1.length())) :
                    alt2.basesMatch(Arrays.copyOf(alt1.getBases(), alt2.length())));
        }
    }

    // compute the likelihoods that each allele is contained at some allele fraction in the sample
    protected PerAlleleCollection<Double> somaticLogOdds(final LikelihoodMatrix<Allele> logMatrix) {
        int alleleListEnd = logMatrix.alleles().size()-1;
        final int nonRefIndex = logMatrix.alleles().contains(Allele.NON_REF_ALLELE)
            && logMatrix.alleles().get(alleleListEnd).equals(Allele.NON_REF_ALLELE) ? alleleListEnd : -1;
        if (logMatrix.alleles().contains(Allele.NON_REF_ALLELE) && !(logMatrix.alleles().get(alleleListEnd).equals(Allele.NON_REF_ALLELE))) {
            throw new IllegalStateException("<NON_REF> must be last in the allele list.");
        }
        final double logEvidenceWithAllAlleles = logMatrix.numberOfReads() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(getAsRealMatrix(logMatrix), MTAC.minAF, nonRefIndex);

        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int refIndex = getRefIndex(logMatrix);
        IntStream.range(0, logMatrix.numberOfAlleles()).filter(a -> a != refIndex).forEach( a -> {
            final Allele allele = logMatrix.getAllele(a);
            final LikelihoodMatrix<Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(logMatrix, allele);
            final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.numberOfReads() == 0 ? 0 :
                    SomaticLikelihoodsEngine.logEvidence(getAsRealMatrix(logMatrixWithoutThisAllele), MTAC.minAF, logMatrixWithoutThisAllele.numberOfAlleles() > 1 ? nonRefIndex-1 : -1);  //nonRefIndex-1 because we're evaluating without one allele; if th
            lods.setAlt(allele, logEvidenceWithAllAlleles - logEvidenceWithoutThisAllele);
        });
        return lods;
    }

    private void addGenotypes(final ReadLikelihoods<Allele> logLikelihoods,
                              final List<Allele> allelesToEmit,
                              final VariantContextBuilder callVcb) {
        final List<Genotype> genotypes = IntStream.range(0, logLikelihoods.numberOfSamples()).mapToObj(n -> {
            final String sample = logLikelihoods.getSample(n);
            final LikelihoodMatrix<Allele> logMatrix = new SubsettedLikelihoodMatrix<>(logLikelihoods.sampleMatrix(n), allelesToEmit);
            final double[] alleleCounts = getEffectiveCounts(logMatrix);
            final double[] flatPriorPseudocounts = new IndexRange(0, logMatrix.numberOfAlleles()).mapToDouble(a -> 1);
            final double[] alleleFractionsPosterior = logMatrix.numberOfReads() == 0 ? flatPriorPseudocounts :
                    SomaticLikelihoodsEngine.alleleFractionsPosterior(getAsRealMatrix(logMatrix), flatPriorPseudocounts);
            final double[] tumorAlleleFractionsMean = MathUtils.normalizeFromRealSpace(alleleFractionsPosterior);

            // TODO: We shouldn't always assume that the genotype in the normal is hom ref
            final Allele ref = logMatrix.getAllele(getRefIndex(logMatrix));
            return new GenotypeBuilder(sample, normalSamples.contains(sample) ? Collections.nCopies(2, ref) : logMatrix.alleles())
                    .AD(Arrays.stream(alleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, Arrays.copyOfRange(tumorAlleleFractionsMean, 1, tumorAlleleFractionsMean.length))
                    .make();
        }).collect(Collectors.toList());

        callVcb.genotypes(genotypes);
    }

    private static double[] getEffectiveCounts(final LikelihoodMatrix<Allele> logLikelihoodMatrix) {
        if (logLikelihoodMatrix.numberOfReads() == 0) {
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
    private PerAlleleCollection<Double> diploidAltLogOdds(final LikelihoodMatrix<Allele> matrix) {
        final int refIndex = getRefIndex(matrix);
        final int numReads = matrix.numberOfReads();
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

    private int getRefIndex(LikelihoodMatrix<Allele> matrix) {
        final OptionalInt optionalRefIndex = IntStream.range(0, matrix.numberOfAlleles()).filter(a -> matrix.getAllele(a).isReference()).findFirst();
        Utils.validateArg(optionalRefIndex.isPresent(), "No ref allele found in likelihoods");
        return optionalRefIndex.getAsInt();
    }

    //convert a likelihood matrix of alleles x reads into a RealMatrix
    public static RealMatrix getAsRealMatrix(final LikelihoodMatrix<Allele> matrix) {
        final RealMatrix result = new Array2DRowRealMatrix(matrix.numberOfAlleles(), matrix.numberOfReads());
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return matrix.get(row, column);
            }
        });
        return result;
    }

    private static LikelihoodMatrix<Allele> combinedLikelihoodMatrix(final List<LikelihoodMatrix<Allele>> matrices, final AlleleList<Allele> alleleList) {
        final List<GATKRead> reads = matrices.stream().flatMap(m -> m.reads().stream()).collect(Collectors.toList());
        final ReadLikelihoods<Allele> combinedLikelihoods = new ReadLikelihoods<>(SampleList.singletonSampleList("COMBINED"), alleleList, ImmutableMap.of("COMBINED", reads));

        int combinedReadIndex = 0;
        final LikelihoodMatrix<Allele> result = combinedLikelihoods.sampleMatrix(0);
        final int alleleCount = result.numberOfAlleles();
        for (final LikelihoodMatrix<Allele> matrix : matrices) {
            final int readCount = matrix.numberOfReads();
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
