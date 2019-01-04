package org.broadinstitute.hellbender.tools.walkers.mutect;

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
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticGenotypingEngine extends AssemblyBasedCallerGenotypingEngine {

    private final M2ArgumentCollection MTAC;
    public final String tumorSample;
    private final String normalSample;
    final boolean hasNormal;
    public static final String DISCARDED_MATE_READ_TAG = "DM";

    // {@link GenotypingEngine} requires a non-null {@link AFCalculatorProvider} but this class doesn't need it.  Thus we make a dummy
    private static final AFCalculatorProvider DUMMY_AF_CALCULATOR_PROVIDER = new AFCalculatorProvider() {
        @Override
        public AFCalculator getInstance(final int ploidy, final int maximumAltAlleles) { return null; }
    };

    @Override
    protected String callSourceString() {
        return "M2_call";
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) { return false; }

    public SomaticGenotypingEngine(final SampleList samples,
                                   final M2ArgumentCollection MTAC,
                                   final String tumorSample,
                                   final String normalSample) {
        super(MTAC, samples, DUMMY_AF_CALCULATOR_PROVIDER, !MTAC.doNotRunPhysicalPhasing);
        this.MTAC = MTAC;
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;
        hasNormal = normalSample != null;
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param log10ReadLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param activeRegionWindow                     Active window
     * @param withBamOut                            whether to annotate reads in readLikelihoods for future writing to bamout
     *
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
        Utils.validateArg(log10ReadLikelihoods.samples().contains(tumorSample), "readLikelihoods does not contain the tumor sample ");

        final List<Haplotype> haplotypes = log10ReadLikelihoods.alleles();

        // Note: passing an empty list of activeAllelesToGenotype is the correct behavior even when givenAlleles is
        // non-empty.  At this point any given alleles have already been injected into the haplotypes, and passing
        // givenAlleles to {@code decomposeHaplotypesIntoVariantContexts} actually overrides any non-given (discovery) alleles, which
        // is not what we want.
        final List<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, assemblyResultSet.getFullReferenceWithPadding(), assemblyResultSet.getPaddedReferenceLoc(), Collections.emptyList(), MTAC.maxMnpDistance).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        if(withBamOut){
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(log10ReadLikelihoods, activeRegionWindow);
        }

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = getVariantContextsFromActiveHaplotypes(loc, haplotypes, false);
            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( mergedVC == null ) {
                continue;
            }
            
            // converting ReadLikelihoods<Haplotype> to ReadLikelihoods<Allele>
            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergedVC, loc, haplotypes);
            final ReadLikelihoods<Allele> log10Likelihoods = log10ReadLikelihoods.marginalize(alleleMapper,
                    new SimpleInterval(mergedVC).expandWithinContig(ALLELE_EXTENSION, header.getSequenceDictionary()));
            filterOverlappingReads(log10Likelihoods, mergedVC.getReference(), loc, false);

            if (emitRefConf) {
                mergedVC = ReferenceConfidenceUtils.addNonRefSymbolicAllele(mergedVC);
                log10Likelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }

            final LikelihoodMatrix<Allele> log10TumorMatrix = log10Likelihoods.sampleMatrix(log10Likelihoods.indexOfSample(tumorSample));
            final Optional<LikelihoodMatrix<Allele>> log10NormalMatrix =
                    getForNormal(() -> log10Likelihoods.sampleMatrix(log10Likelihoods.indexOfSample(normalSample)));

            final PerAlleleCollection<Double> tumorLog10Odds = somaticLog10Odds(log10TumorMatrix);
            final Optional<PerAlleleCollection<Double>> normalLog10Odds = getForNormal(() -> diploidAltLog10Odds(log10NormalMatrix.get()));
            final Optional<PerAlleleCollection<Double>> normalArtifactLog10Odds = getForNormal(() -> somaticLog10Odds(log10NormalMatrix.get()));

            final Set<Allele> forcedAlleles = getAllelesConsistentWithGivenAlleles(givenAlleles, loc, mergedVC);

            final List<Allele> tumorAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> forcedAlleles.contains(allele) || tumorLog10Odds.getAlt(allele) > MTAC.getEmissionLod())
                    .collect(Collectors.toList());

            final long somaticAltCount = tumorAltAlleles.stream()
                    .filter(allele -> forcedAlleles.contains(allele) || !hasNormal || MTAC.genotypeGermlineSites || normalLog10Odds.get().getAlt(allele) > MTAC.normalLod)
                    .count();

            // if every alt allele is germline, skip this variant.  However, if some alt alleles are germline and others
            // are not we emit them all so that the filtering engine can see them
            if (somaticAltCount == 0) {
                continue;
            }


            final List<Allele> allAllelesToEmit = ListUtils.union(Arrays.asList(mergedVC.getReference()), tumorAltAlleles);


            final LikelihoodMatrix<Allele> subsettedLog10TumorMatrix = new SubsettedLikelihoodMatrix<>(log10TumorMatrix, allAllelesToEmit);
            final Optional<LikelihoodMatrix<Allele>> subsettedLog10NormalMatrix =
                    getForNormal(() -> new SubsettedLikelihoodMatrix<>(log10NormalMatrix.get(), allAllelesToEmit));

            final Map<String, Object> negativeLog10PopulationAFAnnotation = GermlineProbabilityCalculator.getNegativeLog10PopulationAFAnnotation(featureContext.getValues(MTAC.germlineResource, loc), tumorAltAlleles, MTAC.getDefaultAlleleFrequency());

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC)
                    .alleles(allAllelesToEmit)
                    .attributes(negativeLog10PopulationAFAnnotation)
                    .attribute(GATKVCFConstants.TUMOR_LOD_KEY, tumorAltAlleles.stream().mapToDouble(tumorLog10Odds::getAlt).toArray());

            normalArtifactLog10Odds.ifPresent(values -> callVcb.attribute(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE, Arrays.stream(values.asDoubleArray(tumorAltAlleles)).map(x->-x).toArray()));

            normalLog10Odds.ifPresent(values -> callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY, values.asDoubleArray(tumorAltAlleles)));

            if (!featureContext.getValues(MTAC.pon, mergedVC.getStart()).isEmpty()) {
                callVcb.attribute(GATKVCFConstants.IN_PON_VCF_ATTRIBUTE, true);
            }

            addGenotypes(subsettedLog10TumorMatrix, subsettedLog10NormalMatrix, callVcb);
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

        final List<VariantContext> outputCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        final int eventCount = outputCalls.size();
        final List<VariantContext> outputCallsWithEventCountAnnotation = outputCalls.stream()
                .map(vc -> new VariantContextBuilder(vc).attribute(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount).make())
                .collect(Collectors.toList());
        return new CalledHaplotypes(outputCallsWithEventCountAnnotation, calledHaplotypes);
    }

    private Set<Allele> getAllelesConsistentWithGivenAlleles(List<VariantContext> givenAlleles, int loc, VariantContext mergedVC) {
        final List<Pair<Allele, Allele>> givenAltAndRefAllelesInOriginalContext =  getVariantContextsFromGivenAlleles(loc, givenAlleles, false).stream()
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

    private void addGenotypes(final LikelihoodMatrix<Allele> tumorLog10Matrix,
                              final Optional<LikelihoodMatrix<Allele>> normalLog10Matrix,
                              final VariantContextBuilder callVcb) {
        final GenotypeBuilder gb = new GenotypeBuilder(tumorSample, tumorLog10Matrix.alleles());
        final double[] flatPriorPseudocounts = new IndexRange(0, tumorLog10Matrix.numberOfAlleles()).mapToDouble(n -> 1);
        final double[] alleleFractionsPosterior = tumorLog10Matrix.numberOfReads() == 0 ? flatPriorPseudocounts :
                SomaticLikelihoodsEngine.alleleFractionsPosterior(getAsRealMatrix(tumorLog10Matrix), flatPriorPseudocounts);
        // Use mean of the allele fraction posterior distribution
        final double[] tumorAlleleFractionsMean = MathUtils.normalizeFromRealSpace(alleleFractionsPosterior);
        gb.attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, Arrays.copyOfRange(tumorAlleleFractionsMean, 1, tumorAlleleFractionsMean.length));
        final Genotype tumorGenotype = gb.make();
        final List<Genotype> genotypes = new ArrayList<>(Arrays.asList(tumorGenotype));

        // TODO: We shouldn't always assume that the genotype in the normal is hom ref
        final Allele ref = tumorLog10Matrix.getAllele(getRefIndex(tumorLog10Matrix));
        final List<Allele> homRefAllelesforNormalGenotype = Collections.nCopies(2, ref);

        // if we are calling with a normal, build the genotype for the sample to appear in vcf
        if (hasNormal) {
            final double[] normalAlleleCounts = getEffectiveCounts(normalLog10Matrix.get());
            final Genotype normalGenotype = new GenotypeBuilder(normalSample, homRefAllelesforNormalGenotype)
                    .AD(Arrays.stream(normalAlleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, getAltAlleleFractions(normalAlleleCounts))
                    .make();
            genotypes.add(normalGenotype);
        }

        callVcb.genotypes(genotypes);
    }

    public static double[] getEffectiveCounts(final LikelihoodMatrix<Allele> log10LikelihoodMatrix) {
        if (log10LikelihoodMatrix.numberOfReads() == 0) {
            return new double[log10LikelihoodMatrix.numberOfAlleles()]; // zero counts for each allele
        }
        final RealMatrix log10Likelihoods = getAsRealMatrix(log10LikelihoodMatrix);
        return MathUtils.sumArrayFunction(0, log10Likelihoods.getColumnDimension(),
                read -> MathUtils.normalizeFromLog10ToLinearSpace(log10Likelihoods.getColumn(read)));
    }

    private static double[] getAltAlleleFractions(final double[] alleleCounts) {
        final double[] allAlleleFractions = MathUtils.normalizeFromRealSpace(alleleCounts);
        return Arrays.copyOfRange(allAlleleFractions, 1, allAlleleFractions.length); //omit the first entry of the array corresponding to the reference
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

    private <E> Optional<E> getForNormal(final Supplier<E> supplier) {
        return hasNormal ? Optional.of(supplier.get()) : Optional.empty();
    }
}
