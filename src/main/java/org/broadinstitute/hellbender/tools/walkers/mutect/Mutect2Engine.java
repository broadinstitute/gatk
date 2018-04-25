package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingGivenAllelesUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by davidben on 9/15/16.
 */
public final class Mutect2Engine implements AssemblyRegionEvaluator {

    private static final String MUTECT_VERSION = "2.1";

    public static final String TUMOR_SAMPLE_KEY_IN_VCF_HEADER = "tumor_sample";
    public static final String NORMAL_SAMPLE_KEY_IN_VCF_HEADER = "normal_sample";

    private static final Logger logger = LogManager.getLogger(Mutect2Engine.class);
    private final static List<VariantContext> NO_CALLS = Collections.emptyList();
    public static final int INDEL_START_QUAL = 30;
    public static final int INDEL_CONTINUATION_QUAL = 10;
    public static final double MAX_ALT_FRACTION_IN_NORMAL = 0.3;
    public static final int MAX_NORMAL_QUAL_SUM = 100;

    private M2ArgumentCollection MTAC;
    private SAMFileHeader header;

    private static final int MIN_READ_LENGTH = 30;
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;
    public static final int MINIMUM_BASE_QUALITY = 6;   // for active region determination

    private final SampleList samplesList;
    private final String tumorSample;
    private final String normalSample;

    private CachingIndexedFastaSequenceFile referenceReader;
    private ReadThreadingAssembler assemblyEngine;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private SomaticGenotypingEngine genotypingEngine;
    private Optional<HaplotypeBAMWriter> haplotypeBAMWriter;
    private VariantAnnotatorEngine annotationEngine;
    private final SmithWatermanAligner aligner;
    private AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer();

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     *
     * @param MTAC command-line arguments for the HaplotypeCaller
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param header header for the reads
     * @param reference path to the reference
     */
    public Mutect2Engine(final M2ArgumentCollection MTAC, final boolean createBamOutIndex, final boolean createBamOutMD5, final SAMFileHeader header, final String reference ) {
        this.MTAC = Utils.nonNull(MTAC);
        this.header = Utils.nonNull(header);
        referenceReader = AssemblyBasedCallerUtils.createReferenceReader(Utils.nonNull(reference));
        aligner = SmithWatermanAligner.getAligner(MTAC.smithWatermanImplementation);
        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(header)));
        tumorSample = decodeSampleNameIfNecessary(MTAC.tumorSample);
        normalSample = MTAC.normalSample == null ? null : decodeSampleNameIfNecessary(MTAC.normalSample);
        checkSampleInBamHeader(tumorSample);
        checkSampleInBamHeader(normalSample);

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(MTAC.defaultGATKVariantAnnotationArgumentCollection, null, Collections.emptyList(), false);
        assemblyEngine = AssemblyBasedCallerUtils.createReadThreadingAssembler(MTAC);
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs);
        genotypingEngine = new SomaticGenotypingEngine(samplesList, MTAC, tumorSample, normalSample);
        genotypingEngine.setAnnotationEngine(annotationEngine);
        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(MTAC, createBamOutIndex, createBamOutMD5, header);
        trimmer.initialize(MTAC.assemblyRegionTrimmerArgs, header.getSequenceDictionary(), MTAC.debug,
                MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, false);
    }

    //default M2 read filters.  Cheap ones come first in order to fail fast.
    public static List<ReadFilter> makeStandardMutect2ReadFilters() {
        return Arrays.asList(new MappingQualityReadFilter(READ_QUALITY_FILTER_THRESHOLD),
                ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE,
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO,
                ReadFilterLibrary.MAPPED,
                ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT,
                ReadFilterLibrary.NOT_DUPLICATE,
                ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK,
                ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT,
                new ReadLengthReadFilter(MIN_READ_LENGTH, Integer.MAX_VALUE),
                ReadFilterLibrary.GOOD_CIGAR,
                new WellformedReadFilter());
    }

    public void writeHeader(final VariantContextWriter vcfWriter, final Set<VCFHeaderLine> defaultToolHeaderLines) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(new VCFHeaderLine("Mutect Version", MUTECT_VERSION));
        headerInfo.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "Warning: unfiltered Mutect 2 calls.  Please run " + FilterMutectCalls.class.getSimpleName() + " to remove false positives."));
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerInfo.addAll(defaultToolHeaderLines);
        GATKVCFConstants.STANDARD_MUTECT_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);

        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ALLELE_FRACTION_KEY));
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));

        headerInfo.add(new VCFHeaderLine(TUMOR_SAMPLE_KEY_IN_VCF_HEADER, tumorSample));
        if (hasNormal()) {
            headerInfo.add(new VCFHeaderLine(NORMAL_SAMPLE_KEY_IN_VCF_HEADER, normalSample));
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo, samplesList.asListOfSamples());
        vcfHeader.setSequenceDictionary(header.getSequenceDictionary());
        vcfWriter.writeHeader(vcfHeader);
    }

    public List<VariantContext> callRegion(final AssemblyRegion originalAssemblyRegion, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        if ( !originalAssemblyRegion.isActive() || originalAssemblyRegion.size() == 0 ) {
            return NO_CALLS;
        }

        final List<VariantContext> givenAlleles = MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ?
                featureContext.getValues(MTAC.alleles).stream().filter(vc -> MTAC.genotypeFilteredAlleles || vc.isNotFiltered()).collect(Collectors.toList()) :
                Collections.emptyList();

        final AssemblyRegion assemblyActiveRegion = AssemblyBasedCallerUtils.assemblyRegionWithWellMappedReads(originalAssemblyRegion, READ_QUALITY_FILTER_THRESHOLD, header);
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyActiveRegion, givenAlleles, MTAC, header, samplesList, logger, referenceReader, assemblyEngine, aligner);
        final SortedSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(originalAssemblyRegion,allVariationEvents);
        if (!trimmingResult.isVariationPresent()) {
            return NO_CALLS;
        }

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        // we might find out after assembly that the "active" region actually has no variants
        if( ! assemblyResult.isVariationPresent() ) {
            return NO_CALLS;
        }

        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();

        final Map<String,List<GATKRead>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        final ReadLikelihoods<Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner);
        readLikelihoods.changeReads(readRealignments);

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods, assemblyResult, referenceContext, regionForGenotyping.getSpan(), featureContext, givenAlleles, header);
        writeBamOutput(assemblyResult, readLikelihoods, calledHaplotypes);
        return calledHaplotypes.getCalls();
    }

    private void writeBamOutput(AssemblyResultSet assemblyResult, ReadLikelihoods<Haplotype> readLikelihoods, HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        if ( haplotypeBAMWriter.isPresent() ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            haplotypeBAMWriter.get().writeReadsAlignedToHaplotypes(
                    assemblyResult.getHaplotypeList(),
                    assemblyResult.getPaddedReferenceLoc(),
                    assemblyResult.getHaplotypeList(),
                    calledHaplotypeSet,
                    readLikelihoods);
        }
    }

    //TODO: should be a variable, not a function
    private boolean hasNormal() {
        return (normalSample != null);
    }

    protected Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return AssemblyBasedCallerUtils.splitReadsBySample(samplesList, header, reads);
    }

    public void shutdown() {
        likelihoodCalculationEngine.close();
        aligner.close();
        haplotypeBAMWriter.ifPresent(writer -> writer.close());
    }

    @Override
    public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext featureContext) {
        if ( MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vcFromAllelesRod = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(featureContext, ref.getInterval(), false, MTAC.genotypeFilteredAlleles, logger, MTAC.alleles);
            if( vcFromAllelesRod != null ) {
                return new ActivityProfileState(ref.getInterval(), 1.0);
            }
        }

        final byte refBase = ref.getBase();
        final SimpleInterval refInterval = ref.getInterval();

        if( context == null || context.getBasePileup().isEmpty() ) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        final ReadPileup pileup = context.getBasePileup();
        final ReadPileup tumorPileup = pileup.getPileupForSample(tumorSample, header);
        final Pair<Integer, Double> tumorAltCountAndQualSum = altCountAndQualSum(tumorPileup, refBase);
        final int tumorAltCount = tumorAltCountAndQualSum.getFirst();
        final int tumorRefCount = tumorPileup.size() - tumorAltCount;

        final double tumorLog10Odds = -QualityUtils.qualToErrorProbLog10(tumorAltCountAndQualSum.getSecond()) +
                MathUtils.log10Factorial(tumorAltCount) + MathUtils.log10Factorial(tumorRefCount) - MathUtils.log10Factorial(tumorPileup.size() + 1);

        if (tumorLog10Odds < MTAC.initialTumorLodThreshold) {
            return new ActivityProfileState(refInterval, 0.0);
        } else if (hasNormal() && !MTAC.genotypeGermlineSites) {
            final ReadPileup normalPileup = pileup.getPileupForSample(normalSample, header);
            final Pair<Integer, Double> normalAltCountAndQualSum = altCountAndQualSum(normalPileup, refBase);
            final int normalAltCount = normalAltCountAndQualSum.getFirst();
            final double normalQualSum = normalAltCountAndQualSum.getSecond();
            if (normalAltCount > normalPileup.size() * MAX_ALT_FRACTION_IN_NORMAL && normalQualSum > MAX_NORMAL_QUAL_SUM) {
                return new ActivityProfileState(refInterval, 0.0);
            }
        } else if (!MTAC.genotypeGermlineSites) {
            final List<VariantContext> germline = featureContext.getValues(MTAC.germlineResource, refInterval);
            if (!germline.isEmpty()){
                final List<Double> germlineAlleleFrequencies = germline.get(0).getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
                if (! germlineAlleleFrequencies.isEmpty() && germlineAlleleFrequencies.get(0) > MTAC.maxPopulationAlleleFrequency) {
                    return new ActivityProfileState(refInterval, 0.0);
                }
            }
        }

        if (!MTAC.genotypePonSites && !featureContext.getValues(MTAC.pon, new SimpleInterval(context.getContig(), (int) context.getPosition(), (int) context.getPosition())).isEmpty()) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        return new ActivityProfileState( refInterval, 1.0, ActivityProfileState.Type.NONE, null);
    }

    private static int getCurrentOrFollowingIndelLength(final PileupElement pe) {
        return pe.isDeletion() ? pe.getCurrentCigarElement().getLength() : pe.getLengthOfImmediatelyFollowingIndel();
    }

    private static double indelQual(final int indelLength) {
        return INDEL_START_QUAL + (indelLength - 1) * INDEL_CONTINUATION_QUAL;
    }

    private static Pair<Integer, Double> altCountAndQualSum(final ReadPileup pileup, final byte refBase) {
        int altCount = 0;
        double qualSum = 0;

        for (final PileupElement pe : pileup) {
            final int indelLength = getCurrentOrFollowingIndelLength(pe);
            if (indelLength > 0) {
                altCount++;
                qualSum += indelQual(indelLength);
            } else if (isNextToUsefulSoftClip(pe)) {
                altCount++;
                qualSum += indelQual(1);
            } else if (pe.getBase() != refBase && pe.getQual() > MINIMUM_BASE_QUALITY) {
                altCount++;
                qualSum += pe.getQual();
            }
        }

        return new Pair<>(altCount, qualSum);
    }

    // check that we're next to a soft clip that is not due to a read that got out of sync and ended in a bunch of BQ2's
    // we only need to check the next base's quality
    private static boolean isNextToUsefulSoftClip(final PileupElement pe) {
        final int offset = pe.getOffset();
        return pe.getQual() > MINIMUM_BASE_QUALITY &&
                ((pe.isBeforeSoftClip() && pe.getRead().getBaseQuality(offset + 1) > MINIMUM_BASE_QUALITY)
                        || (pe.isAfterSoftClip() && pe.getRead().getBaseQuality(offset - 1) > MINIMUM_BASE_QUALITY));
    }

    private void checkSampleInBamHeader(final String sample) {
        if (sample != null && !samplesList.asListOfSamples().contains(sample)) {
            throw new UserException.BadInput("Sample " + sample + " is not in BAM header: " + samplesList.asListOfSamples());
        }
    }

    private String decodeSampleNameIfNecessary(final String name) {
        return samplesList.asListOfSamples().contains(name) ? name : IOUtils.urlDecode(name);
    }
}
