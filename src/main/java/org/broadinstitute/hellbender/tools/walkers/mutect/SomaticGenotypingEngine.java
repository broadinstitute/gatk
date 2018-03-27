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
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import static org.broadinstitute.hellbender.utils.OptimizationUtils.argmax;

public class SomaticGenotypingEngine extends AssemblyBasedCallerGenotypingEngine {

    private final M2ArgumentCollection MTAC;

    public final String tumorSampleName;
    private final String matchedNormalSampleName;
    final boolean hasNormal;

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
                                   final String tumorSampleName,
                                   final String matchedNormalSampleName) {
        super(MTAC, samples, DUMMY_AF_CALCULATOR_PROVIDER, !MTAC.doNotRunPhysicalPhasing);
        this.MTAC = MTAC;
        this.tumorSampleName = tumorSampleName;
        this.matchedNormalSampleName = matchedNormalSampleName;
        hasNormal = matchedNormalSampleName != null;
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param log10ReadLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param activeRegionWindow                     Active window
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
            final SAMFileHeader header) {
        Utils.nonNull(log10ReadLikelihoods, "likelihoods are null");
        Utils.validateArg(log10ReadLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow is null");
        Utils.validateArg(log10ReadLikelihoods.samples().contains(tumorSampleName), "readLikelihoods does not contain the tumor sample ");

        final List<Haplotype> haplotypes = log10ReadLikelihoods.alleles();

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final List<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, assemblyResultSet.getFullReferenceWithPadding(), assemblyResultSet.getPaddedReferenceLoc(), givenAlleles).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        for( final int loc : startPosKeySet ) {
            // Note: passing an empty list of activeAllelesToGenotype is the correct behavior even when givenAlleles is
            // non-empty.  At this point any given alleles have already been injected into the haplotypes, and passing
            // givenAlleles to getVCsAtThisLocation actually overrides any non-given (discovery) alleles, which
            // is not what we want.
            final List<VariantContext> eventsAtThisLoc = getVCsAtThisLocation(haplotypes, loc, Collections.emptyList());
            final VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( mergedVC == null ) {
                continue;
            }
            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(eventsAtThisLoc, mergedVC, loc, haplotypes);

            // converting ReadLikelihoods<Haplotype> to ReadLikeliHoods<Allele>
            final ReadLikelihoods<Allele> log10Likelihoods = log10ReadLikelihoods.marginalize(alleleMapper,
                    new SimpleInterval(mergedVC).expandWithinContig(ALLELE_EXTENSION, header.getSequenceDictionary()));
            filterOverlappingReads(log10Likelihoods, mergedVC.getReference(), loc, false);

            final LikelihoodMatrix<Allele> log10TumorMatrix = log10Likelihoods.sampleMatrix(log10Likelihoods.indexOfSample(tumorSampleName));
            final Optional<LikelihoodMatrix<Allele>> log10NormalMatrix =
                    getForNormal(() -> log10Likelihoods.sampleMatrix(log10Likelihoods.indexOfSample(matchedNormalSampleName)));

            final PerAlleleCollection<Double> tumorLog10Odds = somaticLog10Odds(log10TumorMatrix);
            final Optional<PerAlleleCollection<Double>> normalLog10Odds = getForNormal(() -> diploidAltLog10Odds(log10NormalMatrix.get()));
            final Optional<PerAlleleCollection<Double>> normalArtifactLog10Odds = getForNormal(() -> somaticLog10Odds(log10NormalMatrix.get()));

            final List<Allele> givenAllelesInOriginalContext =  getVCsAtThisLocation(Collections.emptyList(), loc, givenAlleles).stream()
                    .flatMap(vc -> vc.getAlternateAlleles().stream()).collect(Collectors.toList());

            final Set<Allele> allelesConsistentWithGivenAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> givenAllelesInOriginalContext.stream().anyMatch(givenAllele -> allelesAreConsistent(givenAllele, allele)))
                    .collect(Collectors.toSet());

            final List<Allele> somaticAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> allelesConsistentWithGivenAlleles.contains(allele) ||
                            ((tumorLog10Odds.getAlt(allele) > MTAC.emissionLodThreshold) &&
                            (!hasNormal || MTAC.genotypeGermlineSites || normalLog10Odds.get().getAlt(allele) > MTAC.normalLodThreshold)))
                    .collect(Collectors.toList());
            final List<Allele> allSomaticAlleles = ListUtils.union(Arrays.asList(mergedVC.getReference()), somaticAltAlleles);
            if (somaticAltAlleles.isEmpty()) {
                continue;
            }

            final LikelihoodMatrix<Allele> subsettedLog10TumorMatrix = new SubsettedLikelihoodMatrix<>(log10TumorMatrix, allSomaticAlleles);
            final Optional<LikelihoodMatrix<Allele>> subsettedLog10NormalMatrix =
                    getForNormal(() -> new SubsettedLikelihoodMatrix<>(log10NormalMatrix.get(), allSomaticAlleles));

            final Map<String, Object> populationAlleleFreqeuncyAnnotation = GermlineProbabilityCalculator.getPopulationAlleleFrequencyAnnotation(featureContext.getValues(MTAC.germlineResource, loc), somaticAltAlleles, MTAC.afOfAllelesNotInGermlineResource);

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC)
                    .alleles(allSomaticAlleles)
                    .attributes(populationAlleleFreqeuncyAnnotation)
                    .attribute(GATKVCFConstants.TUMOR_LOD_KEY, somaticAltAlleles.stream().mapToDouble(tumorLog10Odds::getAlt).toArray());

            normalLog10Odds.ifPresent(values -> callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY, values.asDoubleArray(somaticAltAlleles)));
            normalArtifactLog10Odds.ifPresent(values -> callVcb.attribute(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE, values.asDoubleArray(somaticAltAlleles)));

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

    // check whether two alleles coming from different variant contexts and with possibly different reference alleles
    // could in fact be the same.  The condition is that one is a prefix of the other
    private boolean allelesAreConsistent(final Allele allele1, final Allele allele2) {
        if (allele1.isSymbolic() || allele2.isSymbolic()) {
            return false;
        } else {
            return allele1.length() < allele2.length() ?
                    allele1.basesMatch(Arrays.copyOf(allele2.getBases(), allele1.length())) :
                    allele2.basesMatch(Arrays.copyOf(allele1.getBases(), allele2.length()));
        }
    }

    // compute the likelihoods that each allele is contained at some allele fraction in the sample
    private PerAlleleCollection<Double> somaticLog10Odds(final LikelihoodMatrix<Allele> log10Matrix) {
        final double log10EvidenceWithAllAlleles = log10Matrix.numberOfReads() == 0 ? 0 :
                SomaticLikelihoodsEngine.log10Evidence(getAsRealMatrix(log10Matrix));

        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int refIndex = getRefIndex(log10Matrix);
        IntStream.range(0, log10Matrix.numberOfAlleles()).filter(a -> a != refIndex).forEach( a -> {
            final Allele allele = log10Matrix.getAllele(a);
            final LikelihoodMatrix<Allele> log10MatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(log10Matrix, allele);
            final double log10EvidenceWithoutThisAllele = log10MatrixWithoutThisAllele.numberOfReads() == 0 ? 0 :
                    SomaticLikelihoodsEngine.log10Evidence(getAsRealMatrix(log10MatrixWithoutThisAllele));
            lods.setAlt(allele, log10EvidenceWithAllAlleles - log10EvidenceWithoutThisAllele);
        });
        return lods;
    }

    private void addGenotypes(final LikelihoodMatrix<Allele> tumorLog10Matrix,
                                        final Optional<LikelihoodMatrix<Allele>> normalLog10Matrix,
                                        final VariantContextBuilder callVcb) {
        final double[] tumorAlleleCounts = getEffectiveCounts(tumorLog10Matrix);
        final Genotype tumorGenotype = new GenotypeBuilder(tumorSampleName, tumorLog10Matrix.alleles())
                .AD(Arrays.stream(tumorAlleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, getAltAlleleFractions(tumorAlleleCounts))
                .make();
        final List<Genotype> genotypes = new ArrayList<>(Arrays.asList(tumorGenotype));

        // TODO: We shouldn't always assume that the genotype in the normal is hom ref
        final Allele ref = tumorLog10Matrix.getAllele(getRefIndex(tumorLog10Matrix));
        final List<Allele> homRefAllelesforNormalGenotype = Collections.nCopies(2, ref);

        // if we are calling with a normal, build the genotype for the sample to appear in vcf
        if (hasNormal) {
            final double[] normalAlleleCounts = getEffectiveCounts(normalLog10Matrix.get());
            final Genotype normalGenotype = new GenotypeBuilder(matchedNormalSampleName, homRefAllelesforNormalGenotype)
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
        return Arrays.copyOfRange(allAlleleFractions, 1, allAlleleFractions.length);
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
            final Map<String, List<ReadLikelihoods<Allele>.BestAllele>> fragments = likelihoods.bestAlleles(sample).stream()
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
