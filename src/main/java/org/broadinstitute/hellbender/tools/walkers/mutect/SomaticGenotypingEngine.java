package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
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
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
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
    public static final String DISCARDED_MATE_READ_TAG = "DM";
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
     * @param log10ReadLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param activeRegionWindow                     Active window
     * @param withBamOut                            whether to annotate reads in readLikelihoods for future writing to bamout
     * @param emitRefConf                           generate reference confidence (GVCF) data?
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     */
    public CalledHaplotypes callMutations(
            final ReadLikelihoods<Haplotype> log10ReadLikelihoods,
            final AssemblyResultSet assemblyResultSet,
            final ReferenceContext referenceContext,
            final SimpleInterval activeRegionWindow,
            final FeatureContext featureContext,
            final List<VariantContext> givenAlleles,
            final SAMFileHeader header,
            final boolean withBamOut,
            final boolean emitRefConf) {
        Utils.nonNull(log10ReadLikelihoods);
        Utils.validateArg(log10ReadLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(activeRegionWindow);

        final List<Haplotype> haplotypes = log10ReadLikelihoods.alleles();

        final List<Integer> startPosKeySet = EventMap.buildEventMapsForHaplotypes(haplotypes, assemblyResultSet.getFullReferenceWithPadding(),
                assemblyResultSet.getPaddedReferenceLoc(), MTAC.debug, MTAC.maxMnpDistance).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        if(withBamOut){
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(log10ReadLikelihoods, activeRegionWindow);
        }

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(loc, haplotypes, false);
            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( mergedVC == null ) {
                continue;
            }

            // converting ReadLikelihoods<Haplotype> to ReadLikelihoods<Allele>
            final Map<Allele, List<Haplotype>> alleleMapper = AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, haplotypes, null);
            final ReadLikelihoods<Allele> log10Likelihoods = log10ReadLikelihoods.marginalize(alleleMapper,
                    new SimpleInterval(mergedVC).expandWithinContig(HaplotypeCallerGenotypingEngine.ALLELE_EXTENSION, header.getSequenceDictionary()));
            filterOverlappingReads(log10Likelihoods, mergedVC.getReference(), loc, false);

            if (emitRefConf) {
                mergedVC = ReferenceConfidenceUtils.addNonRefSymbolicAllele(mergedVC);
                log10Likelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }
            final List<LikelihoodMatrix<Allele>> tumorMatrices = IntStream.range(0, log10Likelihoods.numberOfSamples())
                    .filter(n -> !normalSamples.contains(log10Likelihoods.getSample(n)))
                    .mapToObj(log10Likelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final AlleleList<Allele> alleleList = tumorMatrices.get(0);
            final LikelihoodMatrix<Allele> log10TumorMatrix = combinedLikelihoodMatrix(tumorMatrices, alleleList);
            final PerAlleleCollection<Double> tumorLog10Odds = somaticLog10Odds(log10TumorMatrix);

            final List<LikelihoodMatrix<Allele>> normalMatrices = IntStream.range(0, log10Likelihoods.numberOfSamples())
                    .filter(n -> normalSamples.contains(log10Likelihoods.getSample(n)))
                    .mapToObj(log10Likelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final LikelihoodMatrix<Allele> log10NormalMatrix = combinedLikelihoodMatrix(normalMatrices, alleleList);
            final PerAlleleCollection<Double> normalLog10Odds = diploidAltLog10Odds(log10NormalMatrix);
            final PerAlleleCollection<Double> normalArtifactLog10Odds = somaticLog10Odds(log10NormalMatrix);


            final Set<Allele> forcedAlleles = getAllelesConsistentWithGivenAlleles(givenAlleles, loc, mergedVC);
            final List<Allele> tumorAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> forcedAlleles.contains(allele) || tumorLog10Odds.getAlt(allele) > MTAC.getEmissionLod())
                    .collect(Collectors.toList());

            final long somaticAltCount = tumorAltAlleles.stream()
                    .filter(allele -> forcedAlleles.contains(allele) || !hasNormal || MTAC.genotypeGermlineSites || normalLog10Odds.getAlt(allele) > MTAC.normalLod)
                    .count();

            // if every alt allele is germline, skip this variant.  However, if some alt alleles are germline and others
            // are not we emit them all so that the filtering engine can see them
            if (somaticAltCount == 0) {
                continue;
            }

            final List<Allele> allAllelesToEmit = ListUtils.union(Arrays.asList(mergedVC.getReference()), tumorAltAlleles);


            final Map<String, Object> negativeLog10PopulationAFAnnotation = GermlineProbabilityCalculator.getNegativeLog10PopulationAFAnnotation(featureContext.getValues(MTAC.germlineResource, loc), tumorAltAlleles, MTAC.getDefaultAlleleFrequency());

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC)
                    .alleles(allAllelesToEmit)
                    .attributes(negativeLog10PopulationAFAnnotation)
                    .attribute(GATKVCFConstants.TUMOR_LOD_KEY, tumorAltAlleles.stream().mapToDouble(tumorLog10Odds::getAlt).toArray());

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE,
                        Arrays.stream(normalArtifactLog10Odds.asDoubleArray(tumorAltAlleles)).map(x->-x).toArray());
                callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY,
                        Arrays.stream(normalLog10Odds.asDoubleArray(tumorAltAlleles)).map(x->x).toArray());
            }

            if (!featureContext.getValues(MTAC.pon, mergedVC.getStart()).isEmpty()) {
                callVcb.attribute(GATKVCFConstants.IN_PON_VCF_ATTRIBUTE, true);
            }

            addGenotypes(log10Likelihoods, allAllelesToEmit, callVcb);
            final VariantContext call = callVcb.make();
            final VariantContext trimmedCall = GATKVariantContextUtils.trimAlleles(call, true, true);
            final List<Allele> trimmedAlleles = trimmedCall.getAlleles();
            final List<Allele> untrimmedAlleles = call.getAlleles();
            final Map<Allele, List<Allele>> trimmedToUntrimmedAlleleMap = IntStream.range(0, trimmedCall.getNAlleles()).boxed()
                    .collect(Collectors.toMap(n -> trimmedAlleles.get(n), n -> Arrays.asList(untrimmedAlleles.get(n))));
            final ReadLikelihoods<Allele> trimmedLikelihoods = log10Likelihoods.marginalize(trimmedToUntrimmedAlleleMap);

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
    protected PerAlleleCollection<Double> somaticLog10Odds(final LikelihoodMatrix<Allele> log10Matrix) {
        int alleleListEnd = log10Matrix.alleles().size()-1;
        final int nonRefIndex = log10Matrix.alleles().contains(Allele.NON_REF_ALLELE)
            && log10Matrix.alleles().get(alleleListEnd).equals(Allele.NON_REF_ALLELE) ? alleleListEnd : -1;
        if (log10Matrix.alleles().contains(Allele.NON_REF_ALLELE) && !(log10Matrix.alleles().get(alleleListEnd).equals(Allele.NON_REF_ALLELE))) {
            throw new IllegalStateException("<NON_REF> must be last in the allele list.");
        }
        final double log10EvidenceWithAllAlleles = log10Matrix.numberOfReads() == 0 ? 0 :
                SomaticLikelihoodsEngine.log10Evidence(getAsRealMatrix(log10Matrix), MTAC.minAF, nonRefIndex);

        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int refIndex = getRefIndex(log10Matrix);
        IntStream.range(0, log10Matrix.numberOfAlleles()).filter(a -> a != refIndex).forEach( a -> {
            final Allele allele = log10Matrix.getAllele(a);
            final LikelihoodMatrix<Allele> log10MatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(log10Matrix, allele);
            final double log10EvidenceWithoutThisAllele = log10MatrixWithoutThisAllele.numberOfReads() == 0 ? 0 :
                    SomaticLikelihoodsEngine.log10Evidence(getAsRealMatrix(log10MatrixWithoutThisAllele), MTAC.minAF, log10MatrixWithoutThisAllele.numberOfAlleles() > 1 ? nonRefIndex-1 : -1);  //nonRefIndex-1 because we're evaluating without one allele; if th
            lods.setAlt(allele, log10EvidenceWithAllAlleles - log10EvidenceWithoutThisAllele);
        });
        return lods;
    }

    private void addGenotypes(final ReadLikelihoods<Allele> log10Likelihoods,
                              final List<Allele> allelesToEmit,
                              final VariantContextBuilder callVcb) {
        final List<Genotype> genotypes = IntStream.range(0, log10Likelihoods.numberOfSamples()).mapToObj(n -> {
            final String sample = log10Likelihoods.getSample(n);
            final LikelihoodMatrix<Allele> log10Matrix = new SubsettedLikelihoodMatrix<>(log10Likelihoods.sampleMatrix(n), allelesToEmit);
            final double[] alleleCounts = getEffectiveCounts(log10Matrix);
            final double[] flatPriorPseudocounts = new IndexRange(0, log10Matrix.numberOfAlleles()).mapToDouble(a -> 1);
            final double[] alleleFractionsPosterior = log10Matrix.numberOfReads() == 0 ? flatPriorPseudocounts :
                    SomaticLikelihoodsEngine.alleleFractionsPosterior(getAsRealMatrix(log10Matrix), flatPriorPseudocounts);
            final double[] tumorAlleleFractionsMean = MathUtils.normalizeFromRealSpace(alleleFractionsPosterior);

            // TODO: We shouldn't always assume that the genotype in the normal is hom ref
            final Allele ref = log10Matrix.getAllele(getRefIndex(log10Matrix));
            return new GenotypeBuilder(sample, normalSamples.contains(sample) ? Collections.nCopies(2, ref) : log10Matrix.alleles())
                    .AD(Arrays.stream(alleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, Arrays.copyOfRange(tumorAlleleFractionsMean, 1, tumorAlleleFractionsMean.length))
                    .make();
        }).collect(Collectors.toList());

        callVcb.genotypes(genotypes);
    }

    private static double[] getEffectiveCounts(final LikelihoodMatrix<Allele> log10LikelihoodMatrix) {
        if (log10LikelihoodMatrix.numberOfReads() == 0) {
            return new double[log10LikelihoodMatrix.numberOfAlleles()]; // zero counts for each allele
        }
        final RealMatrix log10Likelihoods = getAsRealMatrix(log10LikelihoodMatrix);
        return MathUtils.sumArrayFunction(0, log10Likelihoods.getColumnDimension(),
                read -> MathUtils.normalizeFromLog10ToLinearSpace(log10Likelihoods.getColumn(read)));
    }

    /**
     * Calculate the log10 likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log10 likelihood to get the log-odds.
     *
     * @param matrix a matrix of log10 likelihoods
     */
    private PerAlleleCollection<Double> diploidAltLog10Odds(final LikelihoodMatrix<Allele> matrix) {
        final int refIndex = getRefIndex(matrix);
        final int numReads = matrix.numberOfReads();
        final double homRefLog10Likelihood = new IndexRange(0, numReads).sum(r -> matrix.get(refIndex,r));

        final PerAlleleCollection<Double> result = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        // hom ref likelihood for the ref allele, het likelihood for alt alleles
        IntStream.range(0, matrix.numberOfAlleles()).filter(a -> a != refIndex)
                .forEach( a -> {
                    final double hetLog10Likelihood = new IndexRange(0, numReads)
                            .sum(r -> MathUtils.log10SumLog10(matrix.get(refIndex, r), matrix.get(a,r)) + MathUtils.LOG10_ONE_HALF);
                    result.setAlt(matrix.getAllele(a), homRefLog10Likelihood - hetLog10Likelihood);
                });
        return result;
    }

    private int getRefIndex(LikelihoodMatrix<Allele> matrix) {
        final OptionalInt optionalRefIndex = IntStream.range(0, matrix.numberOfAlleles()).filter(a -> matrix.getAllele(a).isReference()).findFirst();
        Utils.validateArg(optionalRefIndex.isPresent(), "No ref allele found in likelihoods");
        return optionalRefIndex.getAsInt();
    }

    private void filterOverlappingReads(final ReadLikelihoods<Allele> likelihoods, final Allele ref, final int location, final boolean retainMismatches) {
        for (final String sample : likelihoods.samples()) {
            // Get the best alleles of each read and group them by the read name.
            // This puts paired reads from the same fragment together
            final Map<String, List<ReadLikelihoods<Allele>.BestAllele>> fragments = likelihoods.bestAllelesBreakingTies(sample).stream()
                    .collect(Collectors.groupingBy(ba -> ba.read.getName()));

            // We only potentially filter read pairs that overlap at this position
            final List<Pair<ReadLikelihoods<Allele>.BestAllele, ReadLikelihoods<Allele>.BestAllele>> overlappingReadPairs =
                    fragments.values().stream()
                            .filter(l -> l.size() == 2)
                            .map(l -> new ImmutablePair<>(l.get(0), l.get(1)))
                            .filter(p -> ReadUtils.isInsideRead(p.getLeft().read, location) && ReadUtils.isInsideRead(p.getRight().read, location))
                            .collect(Collectors.toList());

            final Set<GATKRead> readsToDiscard = new HashSet<>();

            for (final Pair<ReadLikelihoods<Allele>.BestAllele, ReadLikelihoods<Allele>.BestAllele> pair : overlappingReadPairs) {
                final ReadLikelihoods<Allele>.BestAllele read = pair.getLeft();
                final ReadLikelihoods<Allele>.BestAllele mate = pair.getRight();

                if (read.allele.equals(mate.allele)) {
                    // keep the higher-quality read
                    readsToDiscard.add(read.likelihood < mate.likelihood ? read.read : mate.read);

                    // mark the read to indicate that its mate was dropped - so that we can account for it in {@link StrandArtifact}
                    // and {@link StrandBiasBySample}
                    if (MTAC.annotateBasedOnReads){
                        final GATKRead readToKeep = read.likelihood >= mate.likelihood ? read.read : mate.read;
                        readToKeep.setAttribute(DISCARDED_MATE_READ_TAG, 1);
                    }
                } else if (retainMismatches) {
                    // keep the alt read
                    readsToDiscard.add(read.allele.equals(ref) ? read.read : mate.read);
                } else {
                    // throw out both
                    readsToDiscard.add(read.read);
                    readsToDiscard.add(mate.read);
                }
            }

            likelihoods.removeSampleReads(likelihoods.indexOfSample(sample), readsToDiscard, likelihoods.numberOfAlleles());
        }
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
}
