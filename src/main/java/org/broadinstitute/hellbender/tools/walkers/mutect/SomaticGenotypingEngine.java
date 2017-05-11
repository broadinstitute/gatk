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

import static org.broadinstitute.hellbender.utils.OptimizationUtils.argmax;

public class SomaticGenotypingEngine extends AssemblyBasedCallerGenotypingEngine {

    public static final String IN_COSMIC_VCF_ATTRIBUTE = "IN_COSMIC";
    public static final String IN_DBSNP_VCF_ATTRIBUTE = "IN_DBSNP";
    public static final String IN_PON_VCF_ATTRIBUTE = "IN_PON";
    public static final String NORMAL_ARTIFACT_LOD_ATTRIBUTE = "N_ART_LOD";

    private final M2ArgumentCollection MTAC;

    public final String tumorSampleName;
    private final String matchedNormalSampleName;
    final boolean hasNormal;

    //Mutect2 does not run in GGA mode
    private static final List<VariantContext> NO_GIVEN_ALLELES = Collections.emptyList();

    // {@link GenotypingEngine} requires a non-null {@link AFCalculatorProvider} but this class doesn't need it.  Thus we make a dummy
    private static final AFCalculatorProvider DUMMY_AF_CALCULATOR_PROVIDER = new AFCalculatorProvider() {
        @Override
        public AFCalculator getInstance(final int ploidy, final int maximumAltAlleles) { return null; }
    };

    private static final Logger logger = Logger.getLogger(SomaticGenotypingEngine.class);

    @Override
    protected String callSourceString() {
        return "M2_call";
    }

    public SomaticGenotypingEngine(final SampleList samples,
                                   final M2ArgumentCollection MTAC,
                                   final String tumorSampleName,
                                   final String matchedNormalSampleName) {
        super(MTAC, samples, DUMMY_AF_CALCULATOR_PROVIDER, !MTAC.doNotRunPhysicalPhasing);
        this.MTAC = MTAC;
        this.tumorSampleName = tumorSampleName;
        this.matchedNormalSampleName = matchedNormalSampleName;
        hasNormal = matchedNormalSampleName != null;

        final double errorProbability = QualityUtils.qualToErrorProb(MTAC.POWER_CONSTANT_QSCORE);
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
            final SAMFileHeader header) {
        Utils.nonNull(log10ReadLikelihoods, "likelihoods are null");
        Utils.validateArg(log10ReadLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow is null");
        Utils.validateArg(log10ReadLikelihoods.samples().contains(tumorSampleName), "readLikelihoods does not contain the tumor sample ");

        final List<Haplotype> haplotypes = log10ReadLikelihoods.alleles();

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final List<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, assemblyResultSet.getFullReferenceWithPadding(), assemblyResultSet.getPaddedReferenceLoc(), NO_GIVEN_ALLELES).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = getVCsAtThisLocation(haplotypes, loc, NO_GIVEN_ALLELES);
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
            final Optional<PerAlleleCollection<Double>> normalArtifactLods = getForNormal(() -> somaticLog10Odds(log10NormalMatrix.get()));

            final List<Allele> somaticAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> tumorLog10Odds.getAlt(allele) > MTAC.TUMOR_LOD_THRESHOLD)
                    .filter(allele -> !hasNormal || normalLog10Odds.get().getAlt(allele) > MTAC.NORMAL_LOD_THRESHOLD)
                    .collect(Collectors.toList());
            final List<Allele> allSomaticAlleles = ListUtils.union(Arrays.asList(mergedVC.getReference()), somaticAltAlleles);
            if (somaticAltAlleles.isEmpty()) {
                continue;
            }

            final LikelihoodMatrix<Allele> subsettedLog10TumorMatrix = new SubsettedLikelihoodMatrix<>(log10TumorMatrix, allSomaticAlleles);
            final Optional<LikelihoodMatrix<Allele>> subsettedLog10NormalMatrix =
                    getForNormal(() -> new SubsettedLikelihoodMatrix<>(log10NormalMatrix.get(), allSomaticAlleles));

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC);
            callVcb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, somaticAltAlleles.stream().mapToDouble(tumorLog10Odds::getAlt).toArray());

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY, somaticAltAlleles.stream().mapToDouble(a -> normalLog10Odds.get().getAlt(a)).toArray());
                callVcb.attribute(NORMAL_ARTIFACT_LOD_ATTRIBUTE, somaticAltAlleles.stream().mapToDouble(a -> normalArtifactLods.get().getAlt(a)).toArray());
            }

            if (!featureContext.getValues(MTAC.cosmicFeatureInput, loc).isEmpty()) {
                callVcb.attribute(IN_COSMIC_VCF_ATTRIBUTE, true);
            }

            if (!featureContext.getValues(MTAC.dbsnp.dbsnp, loc).isEmpty()) {
                callVcb.attribute(IN_DBSNP_VCF_ATTRIBUTE, true);
            }

            if (!featureContext.getValues(MTAC.normalPanelFeatureInput, mergedVC.getStart()).isEmpty()) {
                callVcb.attribute(IN_PON_VCF_ATTRIBUTE, true);
            }

            final VariantContext call = addGenotypes(subsettedLog10TumorMatrix, subsettedLog10NormalMatrix, callVcb);
            final ReadLikelihoods<Allele> log10LikelihoodsForAnotations = prepareReadAlleleLikelihoodsForAnnotation(log10ReadLikelihoods, Collections.emptyMap(),
                    false, alleleMapper, log10Likelihoods, call);
            final VariantContext annotatedCall =  annotationEngine.annotateContext(call, featureContext, referenceContext, log10LikelihoodsForAnotations, a -> true);

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

    private VariantContext addGenotypes(final LikelihoodMatrix<Allele> tumorLog10Matrix,
                                        final Optional<LikelihoodMatrix<Allele>> normalLog10Matrix,
                                        final VariantContextBuilder callVcb) {
        final double[] tumorAlleleCounts = SomaticLikelihoodsEngine.getEffectiveCounts(getAsRealMatrix(tumorLog10Matrix));
        final Genotype tumorGenotype = new GenotypeBuilder(tumorSampleName, tumorLog10Matrix.alleles())
                .AD(Arrays.stream(tumorAlleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, MathUtils.normalizeFromRealSpace(tumorAlleleCounts))
                .make();
        final List<Genotype> genotypes = new ArrayList<>(Arrays.asList(tumorGenotype));

        // TODO: We shouldn't always assume that the genotype in the normal is hom ref
        final Allele ref = tumorLog10Matrix.getAllele(getRefIndex(tumorLog10Matrix));
        final List<Allele> homRefAllelesforNormalGenotype = Collections.nCopies(2, ref);

        // if we are calling with a normal, build the genotype for the sample to appear in vcf
        if (hasNormal) {
            final double[] normalAlleleCounts = SomaticLikelihoodsEngine.getEffectiveCounts(getAsRealMatrix(normalLog10Matrix.get()));
            final Genotype normalGenotype = new GenotypeBuilder(matchedNormalSampleName, homRefAllelesforNormalGenotype)
                    .AD(Arrays.stream(normalAlleleCounts).mapToInt(x -> (int) FastMath.round(x)).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, MathUtils.normalizeFromRealSpace(normalAlleleCounts))
                    .make();
            genotypes.add(normalGenotype);
        }

        return new VariantContextBuilder(callVcb).alleles(tumorLog10Matrix.alleles()).genotypes(genotypes).make();
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
