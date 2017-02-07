package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
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
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static java.lang.Math.pow;

/**
 * Created by davidben on 9/15/16.
 */
public final class Mutect2Engine implements AssemblyRegionEvaluator {
    //TODO: move these lists to GATKVCFConstants
    public static final List<String> STANDARD_M_2_INFO_FIELDS = Arrays.asList(GATKVCFConstants.NORMAL_LOD_KEY, GATKVCFConstants.TUMOR_LOD_KEY,
            GATKVCFConstants.PANEL_OF_NORMALS_COUNT_KEY, GATKVCFConstants.HAPLOTYPE_COUNT_KEY, GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY);
    public static final List<String> STRAND_ARTIFACT_INFO_FIELDS = Arrays.asList(GATKVCFConstants.TLOD_FWD_KEY, GATKVCFConstants.TLOD_REV_KEY,
            GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY, GATKVCFConstants.TUMOR_SB_POWER_REV_KEY);

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

    private final Predicate<GATKRead> useReadForGenotyping;

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

        final Predicate<GATKRead> goodReadLengthForGenotyping = read -> read.getLength() >= MIN_READ_LENGTH;
        final Predicate<GATKRead> goodMappingQuality = read -> read.getMappingQuality() >= MTAC.MIN_MAPPING_QUALITY_SCORE;
        final Predicate<GATKRead> isInReadGroupsToKeep = read ->  MTAC.keepRG == null || read.getReadGroup().equals(MTAC.keepRG);
        useReadForGenotyping = goodReadLengthForGenotyping.and(goodMappingQuality).
                and(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE).
                and(isInReadGroupsToKeep);

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

        //TODO: add required annotations for all the filtering steps I propose
        //TODO: how about add this annotation regardless but *filter* optionally?
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            MTAC.annotationsToUse.add("ClusteredReadPosition");
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
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.PRIMARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());

        return filters;
    }


    public void writeHeader(final VariantContextWriter vcfWriter, final SAMSequenceDictionary sequenceDictionary) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

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

        STANDARD_M_2_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);
        STRAND_ARTIFACT_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);

        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.IN_COSMIC_VCF_ATTRIBUTE, 0, VCFHeaderLineType.Flag, "site found in COSMIC database"));
        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.IN_DBSNP_VCF_ATTRIBUTE, 0, VCFHeaderLineType.Flag, "site found in dbSNP database"));
        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.IN_PON_VCF_ATTRIBUTE, 0, VCFHeaderLineType.Flag, "site found in panel of normals"));
        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.NORMAL_ARTIFACT_LOD_ATTRIBUTE, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "log odds of artifact in normal with same allele fraction as tumor"));


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

        final AssemblyRegion assemblyActiveRegion = AssemblyBasedCallerUtils.assemblyRegionWithWellMappedReads(originalAssemblyRegion, MTAC.MIN_MAPPING_QUALITY_SCORE, header);
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

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalAssemblyRegion?
        final Collection<GATKRead> readsToRemove = regionForGenotyping.getReads().stream().filter(useReadForGenotyping.negate()).collect(Collectors.toList());
        regionForGenotyping.removeAll(readsToRemove);

        // we might have no reads left after filtering
        if (regionForGenotyping.getReads().isEmpty()){
            return NO_CALLS;
        }

        // TODO: this quantity and all downstream uses of it seems like it can be obtained from
        // TODO: ReadLikelihoods<Allele>::sampleReads
        final Map<String, List<GATKRead>> perSampleFilteredReadList = splitReadsBySample(readsToRemove);

        final Map<String,List<GATKRead>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        final ReadLikelihoods<Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods,
                perSampleFilteredReadList,
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

    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @return genotype likelihoods of [AA,AB]
     */
    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadPileup pileup, final byte refBase, final byte minBaseQual, final double f) {
        double homRefLog10Likelihood = 0;
        double hetLog10Likelihood = 0;
        for( final PileupElement p : pileup ) {
            if( p.isDeletion() || p.getQual() > minBaseQual ) {
                // TODO: why not use base qualities here?
                //double pobs = QualityUtils.qualToErrorProbLog10(qual);
                final double pobs = 1 - pow(10, (30 / -10.0));
                if( isNonRef(refBase, p)) {
                    hetLog10Likelihood += Math.log10(f*pobs + (1-f)*pobs/3);
                    homRefLog10Likelihood += Math.log10((1-pobs)/3);
                } else {
                    hetLog10Likelihood += Math.log10(f*(1-pobs)/3 + (1-f)*pobs);
                    homRefLog10Likelihood += Math.log10(pobs);
                }
            }
        }

        return new double[]{homRefLog10Likelihood, hetLog10Likelihood};
    }

    protected int getCountOfNonRefEvents(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        return (int) StreamSupport.stream(pileup.spliterator(), false)
                .filter(p -> isNonRef(refBase, p))
                .filter(p -> p.isDeletion() || p.getQual() >= minBaseQual)
                .count();
    }

    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        final double f = calculateF(pileup, refBase, minBaseQual);
        return calcGenotypeLikelihoodsOfRefVsAny(pileup, refBase, minBaseQual, f);
    }

    private double calculateF(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        int totalCount = 0, altCount = 0;
        for( final PileupElement p : pileup ) {
            if( p.isDeletion() || p.getQual() > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    altCount++;
                }
                totalCount++;
            }
        }
        return (double) altCount / totalCount;
    }

    private boolean isNonRef(final byte refBase, final PileupElement p) {
        return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
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
        if( context == null || context.getBasePileup().isEmpty() ) {
            return new ActivityProfileState(ref.getInterval(), 0.0);
        }

        final Map<String, AlignmentContext> splitContexts = context.splitContextBySampleName(header);
        final AlignmentContext tumorContext = splitContexts.get(MTAC.tumorSampleName);
        final AlignmentContext normalContext = splitContexts.get(MTAC.normalSampleName);

        // if there are no tumor reads... there is no activity!
        if (tumorContext == null) {
            return new ActivityProfileState(ref.getInterval(), 0);
        }

        final ReadPileup tumorPileup = tumorContext.getBasePileup().makeFilteredPileup(el -> el.getMappingQual() >= MTAC.MIN_MAPPING_QUALITY_SCORE);
        final double[] tumorGLs = calcGenotypeLikelihoodsOfRefVsAny(tumorPileup, ref.getBase(), MTAC.minBaseQualityScore);
        final double tumorLod = tumorGLs[1] - tumorGLs[0];

        // NOTE: do I want to convert to a probability (or just keep this as a LOD score)

        // also at this point, we should throw out noisy sites (hence the nonRefInNormalCheck) but this is non-optimal
        double prob = 0;
        if (tumorLod > MTAC.INITIAL_TUMOR_LOD_THRESHOLD) {

            // TODO: should we even do this performance optimization?
            // in any case, we have to handle the case where there is no normal (and thus no normal context) which is
            // different than having a normal but having no reads (where we should not enter the active region)
            // TODO: as mentioned above, when we provide the normal bam but there are no normal reads in the region, we call it active. This does not seem right to me.
            if (hasNormal() && normalContext != null) {
                final int nonRefInNormal = getCountOfNonRefEvents(normalContext.getBasePileup(), ref.getBase(), MTAC.minBaseQualityScore);

                final double[] normalGLs = calcGenotypeLikelihoodsOfRefVsAny(normalContext.getBasePileup(), ref.getBase(), MTAC.minBaseQualityScore, 0.5f);
                final double normalLod = normalGLs[0] - normalGLs[1];

                // TODO: parameterize these
                if (normalLod > 1.0 && nonRefInNormal < 4) {
                    prob = 1;
                    logger.debug("At " + ref.getInterval().toString() + " tlod: " + tumorLod + " nlod: " + normalLod + " with normal non-ref of " + nonRefInNormal);
                }
            } else {
                prob = 1;
                logger.debug("At " + ref.getInterval().toString() + " tlod: " + tumorLod + " and no-normal calling");
            }
        }

        return new ActivityProfileState( ref.getInterval(), prob, ActivityProfileState.Type.NONE, null);
    }
}
