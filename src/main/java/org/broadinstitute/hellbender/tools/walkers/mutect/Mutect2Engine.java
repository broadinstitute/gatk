package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.StreamSupport;

/**
 * Created by davidben on 9/15/16.
 */
public final class Mutect2Engine implements AssemblyRegionEvaluator {

    private static final String MUTECT_VERSION = "2.1-beta";

    public static final String TUMOR_SAMPLE_KEY_IN_VCF_HEADER = "tumor_sample";
    public static final String NORMAL_SAMPLE_KEY_IN_VCF_HEADER = "normal_sample";

    private static final Logger logger = LogManager.getLogger(Mutect2Engine.class);
    private final static List<VariantContext> NO_CALLS = Collections.emptyList();

    private M2ArgumentCollection MTAC;
    private SAMFileHeader header;

    private static final int MIN_READ_LENGTH = 30;
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;

    private SampleList samplesList;

    private CachingIndexedFastaSequenceFile referenceReader;
    private ReadThreadingAssembler assemblyEngine;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private SomaticGenotypingEngine genotypingEngine;
    private Optional<HaplotypeBAMWriter> haplotypeBAMWriter;
    private VariantAnnotatorEngine annotationEngine;

    private AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer();

    private static ReadFilter GOOD_READ_LENGTH_FILTER = new ReadFilter() {
        private static final long serialVersionUID = 985763L;
        @Override
        public boolean test(final GATKRead read) {
            return read.getLength() >= MIN_READ_LENGTH;
        }
    };

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     *
     * @param MTAC command-line arguments for the HaplotypeCaller
     * @param header header for the reads
     * @param reference path to the reference
     */
    public Mutect2Engine(final M2ArgumentCollection MTAC, final SAMFileHeader header, final String reference ) {
        this.MTAC = Utils.nonNull(MTAC);
        this.header = Utils.nonNull(header);
        Utils.nonNull(reference);
        referenceReader = AssemblyBasedCallerUtils.createReferenceReader(reference);
        initialize();
    }

    private void initialize() {

        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(header)));
        if (!samplesList.asListOfSamples().contains(MTAC.tumorSampleName)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given tumor" +
                    " sample name " + MTAC.tumorSampleName);
        } else if (MTAC.normalSampleName != null && !samplesList.asListOfSamples().contains(MTAC.normalSampleName)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given normal" +
                    " sample name " + MTAC.normalSampleName);
        }

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(MTAC.annotationGroupsToUse,
                MTAC.annotationsToUse,
                MTAC.annotationsToExclude,
                MTAC.dbsnp.dbsnp,
                MTAC.comps);

        assemblyEngine = AssemblyBasedCallerUtils.createReadThreadingAssembler(MTAC);
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs);
        genotypingEngine = new SomaticGenotypingEngine(samplesList, MTAC, MTAC.tumorSampleName, MTAC.normalSampleName);
        genotypingEngine.setAnnotationEngine(annotationEngine);
        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(MTAC, header);

        trimmer.initialize(MTAC.assemblyRegionTrimmerArgs, header.getSequenceDictionary(), MTAC.debug,
                MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, false);

        if( MTAC.CONTAMINATION_FRACTION_FILE != null ) {
            MTAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(MTAC.CONTAMINATION_FRACTION_FILE, MTAC.CONTAMINATION_FRACTION, samplesList.asSetOfSamples(), logger));
        }
    }

    /**
     * @return the default set of read filters for use with Mutect2
     */
    public static List<ReadFilter> makeStandardMutect2ReadFilters() {
        // The order in which we apply filters is important. Cheap filters come first so we fail fast
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(new MappingQualityReadFilter(READ_QUALITY_FILTER_THRESHOLD));
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.PRIMARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(GOOD_READ_LENGTH_FILTER);
        filters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());

        return filters;
    }

    public void writeHeader(final VariantContextWriter vcfWriter, final SAMSequenceDictionary sequenceDictionary,
                            final Set<VCFHeaderLine>  defaultToolHeaderLines) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "Warning: unfiltered Mutect 2 calls.  Please run " + FilterMutectCalls.class.getSimpleName() + " to remove false positives."));

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());
        headerInfo.addAll(defaultToolHeaderLines);

        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        headerInfo.addAll(getM2HeaderLines());
        headerInfo.addAll(getSampleHeaderLines());

        final VCFHeader vcfHeader = new VCFHeader(headerInfo, samplesList.asListOfSamples());
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        vcfWriter.writeHeader(vcfHeader);
    }

    private Set<VCFHeaderLine> getM2HeaderLines(){
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        headerInfo.add(new VCFHeaderLine("Mutect Version", MUTECT_VERSION));

        GATKVCFConstants.STANDARD_MUTECT_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);

        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ALLELE_FRACTION_KEY));

        if ( ! MTAC.doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }
        return headerInfo;
    }

    private Set<VCFHeaderLine> getSampleHeaderLines(){
        final Set<VCFHeaderLine> sampleLines = new HashSet<>();
        if (hasNormal()) {
            sampleLines.add(new VCFHeaderLine(NORMAL_SAMPLE_KEY_IN_VCF_HEADER, MTAC.normalSampleName));
        }
        sampleLines.add(new VCFHeaderLine(TUMOR_SAMPLE_KEY_IN_VCF_HEADER, MTAC.tumorSampleName));
        return sampleLines;
    }

    public List<VariantContext> callRegion(final AssemblyRegion originalAssemblyRegion, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        if ( MTAC.justDetermineActiveRegions || !originalAssemblyRegion.isActive() || originalAssemblyRegion.size() == 0 ) {
            return NO_CALLS;
        }

        final AssemblyRegion assemblyActiveRegion = AssemblyBasedCallerUtils.assemblyRegionWithWellMappedReads(originalAssemblyRegion, READ_QUALITY_FILTER_THRESHOLD, header);
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyActiveRegion, Collections.emptyList(), MTAC, header, samplesList, logger, referenceReader, assemblyEngine);
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
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods,
                assemblyResult,
                referenceContext,
                regionForGenotyping.getSpan(),
                featureContext,
                header);

        writeBamOutput(assemblyResult, readLikelihoods, calledHaplotypes);

        if( MTAC.debug) { logger.info("----------------------------------------------------------------------------------"); }
        return calledHaplotypes.getCalls();
    }

    private void writeBamOutput(AssemblyResultSet assemblyResult, ReadLikelihoods<Haplotype> readLikelihoods, HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        if ( haplotypeBAMWriter.isPresent() ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (MTAC.disableOptimizations) {
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            }
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
        return (MTAC.normalSampleName != null);
    }

    protected Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return AssemblyBasedCallerUtils.splitReadsBySample(samplesList, header, reads);
    }

    /**
     * Shutdown this M2 engine, closing resources as appropriate
     */
    public void shutdown() {
        likelihoodCalculationEngine.close();

        if ( haplotypeBAMWriter.isPresent() ) {
            haplotypeBAMWriter.get().close();
        }
    }

    @Override
    public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext featureContext) {
        final byte refBase = ref.getBase();
        final SimpleInterval refInterval = ref.getInterval();
        if( context == null || context.getBasePileup().isEmpty() ) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        // because new pileups must be allocated when getting the tumor and normal pileups, we first
        // opportunistically check whether the combined pileup has no evidence of variation, which we will define as
        // having at most one variant read
        final int totalNonRef = countNonRef(refBase, context);
        if (totalNonRef < MTAC.minVariantsInPileup) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        final Map<String, AlignmentContext> splitContexts = context.splitContextBySampleName(header);
        final AlignmentContext tumorContext = splitContexts.get(MTAC.tumorSampleName);
        final AlignmentContext normalContext = splitContexts.get(MTAC.normalSampleName);

        // if there are no tumor reads... there is no activity!
        if (tumorContext == null) {
            return new ActivityProfileState(refInterval, 0);
        }

        final int tumorNonRef = countNonRef(refBase, tumorContext);
        if (tumorNonRef < MTAC.minVariantsInPileup) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        // since errors are rare, the number of errors (if reads are independent) is approximately a Poisson random variable,
        // with mean equal to its variance
        final double expectedTumorNonRefDueToError = StreamSupport.stream(tumorContext.getBasePileup().spliterator(), false)
                .mapToDouble(pe -> QualityUtils.qualToErrorProb(pe.getQual()))
                .sum();
        final double tumorNonRefStdev = Math.sqrt(expectedTumorNonRefDueToError);

        if (tumorNonRef < expectedTumorNonRefDueToError + MTAC.tumorStandardDeviationsThreshold * tumorNonRefStdev) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        if (hasNormal() && normalContext != null && countNonRef(refBase, normalContext) > normalContext.getBasePileup().size() * MTAC.minNormalVariantFraction) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        return new ActivityProfileState( refInterval, 1.0, ActivityProfileState.Type.NONE, null);
    }

    private boolean isNonRef(final byte refBase, final PileupElement p) {
        return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
    }

    private int countNonRef(byte refBase, AlignmentContext context) {
        return context.getBasePileup().getNumberOfElements(p -> isNonRef(refBase, p));
    }
}
