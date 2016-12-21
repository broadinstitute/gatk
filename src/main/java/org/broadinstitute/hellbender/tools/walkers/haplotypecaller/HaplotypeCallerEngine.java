package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;

import java.io.File;
import java.util.*;

/**
 * The core engine for the HaplotypeCaller that does all of the actual work of the tool.
 *
 * Usage:
 * -Pass the HaplotypeCaller args into the constructor, which will initialize the HC engine completely.
 * -Get the appropriate VCF or GVCF writer (depending on our arguments) from {@link #makeVCFWriter}
 * -Write the appropriate VCF header via {@link #writeHeader}
 * -Repeatedly call {@link #isActive} to identify active vs. inactive regions
 * -Repeatedly call {@link #callRegion} to call variants in each region, and add them to your writer
 * -When done, call {@link #shutdown}. Close the writer you got from {@link #makeVCFWriter} yourself.
 */
public final class HaplotypeCallerEngine implements AssemblyRegionEvaluator {

    private static final Logger logger = LogManager.getLogger(HaplotypeCallerEngine.class);

    private final HaplotypeCallerArgumentCollection hcArgs;

    private final SAMFileHeader readsHeader;

    private ReferenceConfidenceModel referenceConfidenceModel = null;

    private AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer();

    // the genotyping engine for the isActive() determination
    private MinimalGenotypingEngine activeRegionEvaluationGenotyperEngine = null;

    private ReadThreadingAssembler assemblyEngine = null;

    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    private HaplotypeCallerGenotypingEngine genotypingEngine = null;

    private VariantAnnotatorEngine annotationEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    private final ReferenceSequenceFile referenceReader;

    // writes Haplotypes to a bam file when the -bamout option is specified
    private Optional<HaplotypeBAMWriter> haplotypeBAMWriter;

    private Set<String> sampleSet;
    private SampleList samplesList;

    private byte minTailQuality;
    public static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    /**
     * Minimum (exclusive) average number of high quality bases per soft-clip to consider that a set of soft-clips is a
     * high quality set.
     */
    private static final double AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD = 6.0;

    /**
     * Maximum-mininum confidence on a variant to exist to consider the position as a potential variant harbouring locus
     * when looking for active regions.
     */
    private static final double MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY = 4.0;

    /**
     * Minimum ploidy assumed when looking for loci that may harbour variation to identify active regions.
     * <p>
     * By default we take the putative ploidy provided by the user, but this turned out to be too insensitive
     * for low ploidy, notoriously with haploid samples. Therefore we impose this minimum.
     * </p>
     */
    private static final int MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY = 2;


    /**
     * Reads with length lower than this number, after clipping off overhands outside the active region,
     * won't be considered for genotyping.
     */
    private static final int READ_LENGTH_FILTER_THRESHOLD = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;

    private static final List<VariantContext> NO_CALLS = Collections.emptyList();

    private static final Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private static final Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     *
     * @param hcArgs command-line arguments for the HaplotypeCaller
     * @param readsHeader header for the reads
     * @param referenceReader reader to provide reference data
     */
    public HaplotypeCallerEngine( final HaplotypeCallerArgumentCollection hcArgs, final SAMFileHeader readsHeader, ReferenceSequenceFile referenceReader ) {
        this.hcArgs = Utils.nonNull(hcArgs);
        this.readsHeader = Utils.nonNull(readsHeader);
        this.referenceReader = Utils.nonNull(referenceReader);

        initialize();
    }

    private void initialize() {
        // Note: order of operations matters here!

        initializeSamples();

        // Must be called after initializeSamples()
        validateAndInitializeArgs();
        minTailQuality = (byte)(hcArgs.minBaseQualityScore - 1);

        initializeActiveRegionEvaluationGenotyperEngine();

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(hcArgs.annotationGroupsToUse, hcArgs.annotationsToUse, hcArgs.annotationsToExclude, hcArgs.dbsnp.dbsnp, hcArgs.comps);

        genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, FixedAFCalculatorProvider.createThreadSafeProvider(hcArgs), ! hcArgs.doNotRunPhysicalPhasing);
        genotypingEngine.setAnnotationEngine(annotationEngine);

        referenceConfidenceModel = new ReferenceConfidenceModel(samplesList, readsHeader, hcArgs.indelSizeToEliminateInRefModel);

        //Allele-specific annotations are not yet supported in the VCF mode
        if (isAlleleSpecificMode(annotationEngine) && isVCFMode()){
           throw new UserException("Allele-specific annotations are not yet supported in the VCF mode");
        }

        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(hcArgs, readsHeader);
        assemblyEngine = AssemblyBasedCallerUtils.createReadThreadingAssembler(hcArgs);
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs);

        trimmer.initialize(hcArgs.assemblyRegionTrimmerArgs, readsHeader.getSequenceDictionary(), hcArgs.debug,
                hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, emitReferenceConfidence());
    }

    private boolean isVCFMode() {
        return hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.NONE;
    }

    private boolean isAlleleSpecificMode(final VariantAnnotatorEngine annotationEngine) {
        //HACK. Note: we can't use subclass information from ReducibleAnnotation (which would be the obvious choice)
        // because RMSMappingQuality is both a reducible annotation and a standard annotation.

        return annotationEngine.getInfoAnnotations().stream()
                .anyMatch(infoFieldAnnotation -> infoFieldAnnotation.getClass().getSimpleName().startsWith("AS_")) ||
                annotationEngine.getGenotypeAnnotations().stream()
                .anyMatch(genotypeAnnotation -> genotypeAnnotation.getClass().getSimpleName().startsWith("AS_"));
    }

    private void validateAndInitializeArgs() {
        if ( hcArgs.genotypeArgs.samplePloidy != HomoSapiensConstants.DEFAULT_PLOIDY && ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Currently, physical phasing is only available for diploid samples.");
        }

        if ( hcArgs.dontGenotype && emitReferenceConfidence() ) {
            throw new UserException("You cannot request gVCF output and 'do not genotype' at the same time");
        }

        if ( emitReferenceConfidence() ) {
            if ( hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) {
                throw new CommandLineException.BadArgumentValue("ERC/gt_mode", "you cannot request reference confidence output and GENOTYPE_GIVEN_ALLELES at the same time");
            }

            hcArgs.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = -0.0;
            hcArgs.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = -0.0;

            // also, we don't need to output several of the annotations
            hcArgs.annotationsToExclude.add(ChromosomeCounts.class.getSimpleName());
            hcArgs.annotationsToExclude.add(FisherStrand.class.getSimpleName());
            hcArgs.annotationsToExclude.add(StrandOddsRatio.class.getSimpleName());
            hcArgs.annotationsToExclude.add(QualByDepth.class.getSimpleName());

            // but we definitely want certain other ones
            hcArgs.annotationsToUse.add(StrandBiasBySample.class.getSimpleName());
            logger.info("Standard Emitting and Calling confidence set to 0.0 for reference-model confidence output");
            if ( ! hcArgs.annotateAllSitesWithPLs ) {
                logger.info("All sites annotated with PLs forced to true for reference-model confidence output");
            }
            hcArgs.annotateAllSitesWithPLs = true;
        } else if ( ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Disabling physical phasing, which is supported only for reference-model confidence output");
        }

        if ( hcArgs.CONTAMINATION_FRACTION_FILE != null ) {
            hcArgs.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(hcArgs.CONTAMINATION_FRACTION_FILE, hcArgs.CONTAMINATION_FRACTION, sampleSet, logger));
        }

        if ( hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES && hcArgs.assemblerArgs.consensusMode ) {
            throw new UserException("HaplotypeCaller cannot be run in both GENOTYPE_GIVEN_ALLELES mode and in consensus mode at the same time. Please choose one or the other.");
        }

        Utils.validateArg(hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD >= QualityUtils.MIN_USABLE_Q_SCORE, "BASE_QUALITY_SCORE_THRESHOLD must be greater than or equal to " + QualityUtils.MIN_USABLE_Q_SCORE + " (QualityUtils.MIN_USABLE_Q_SCORE)");

        if ( emitReferenceConfidence() && samplesList.numberOfSamples() != 1 ) {
            throw new CommandLineException.BadArgumentValue("--emitRefConfidence", "Can only be used in single sample mode currently. Use the sample_name argument to run on a single sample out of a multi-sample BAM file.");
        }
    }

    private void initializeSamples() {
        sampleSet = ReadUtils.getSamplesFromHeader(readsHeader);
        samplesList = new IndexedSampleList(sampleSet);

        if ( hcArgs.sampleNameToUse != null ) {
            if ( ! sampleSet.contains(hcArgs.sampleNameToUse) ) {
                throw new CommandLineException.BadArgumentValue("--sample_name", "Specified name does not exist in input bam files");
            }
            if (sampleSet.size() == 1) {
                //No reason to incur performance penalty associated with filtering if they specified the name of the only sample
                hcArgs.sampleNameToUse = null;
            } else {
                samplesList = new IndexedSampleList(hcArgs.sampleNameToUse);
                sampleSet.clear();
                sampleSet.add(hcArgs.sampleNameToUse);
            }
        }
    }

    private void initializeActiveRegionEvaluationGenotyperEngine() {
        // create a UAC but with the exactCallsLog = null, so we only output the log for the HC caller itself, if requested
        final UnifiedArgumentCollection simpleUAC = new UnifiedArgumentCollection();
        simpleUAC.copyStandardCallerArgsFrom(hcArgs);

        simpleUAC.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
        simpleUAC.genotypingOutputMode = GenotypingOutputMode.DISCOVERY;
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, hcArgs.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, hcArgs.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.CONTAMINATION_FRACTION = 0.0;
        simpleUAC.CONTAMINATION_FRACTION_FILE = null;
        simpleUAC.exactCallsLog = null;
        // Seems that at least with some test data we can lose genuine haploid variation if we use
        // UGs engine with ploidy == 1
        simpleUAC.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, hcArgs.genotypeArgs.samplePloidy);

        activeRegionEvaluationGenotyperEngine = new MinimalGenotypingEngine(simpleUAC, samplesList,
                FixedAFCalculatorProvider.createThreadSafeProvider(simpleUAC));
        activeRegionEvaluationGenotyperEngine.setLogger(logger);
    }

    /**
     * @return the default set of read filters for use with the HaplotypeCaller
     */
    public static List<ReadFilter> makeStandardHCReadFilters() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(new MappingQualityReadFilter(READ_QUALITY_FILTER_THRESHOLD));
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.PRIMARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());

        return filters;
    }

    /**
     * Create a VCF or GVCF writer as appropriate, given our arguments
     *
     * @param outputVCF location to which the vcf should be written
     * @param readsDictionary sequence dictionary for the reads
     * @return a VCF or GVCF writer as appropriate, ready to use
     */
    public VariantContextWriter makeVCFWriter( final String outputVCF, final SAMSequenceDictionary readsDictionary ) {
        Utils.nonNull(outputVCF);
        Utils.nonNull(readsDictionary);

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(new File(outputVCF), readsDictionary, false);

        if ( hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
            try {
                writer = new GVCFWriter(writer, hcArgs.GVCFGQBands, hcArgs.genotypeArgs.samplePloidy);
            } catch ( IllegalArgumentException e ) {
                throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }
        }

        return writer;
    }

    /**
     * Writes an appropriate VCF header, given our arguments, to the provided writer
     *
     * @param vcfWriter writer to which the header should be written
     */
    public void writeHeader( final VariantContextWriter vcfWriter, final SAMSequenceDictionary sequenceDictionary ) {
        Utils.nonNull(vcfWriter);

        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        headerInfo.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());
        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());
        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.DOWNSAMPLED_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        if ( ! hcArgs.doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        if ( emitReferenceConfidence() ) {
            headerInfo.addAll(referenceConfidenceModel.getVCFHeaderLines());
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo, sampleSet);
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        vcfWriter.writeHeader(vcfHeader);
    }


    /**
     * Given a pileup, returns an ActivityProfileState containing the probability (0.0 to 1.0) that it's an "active" site.
     *
     * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
     * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
     * of probabilities. This is consistent with GATK 3.
     *
     * @param context reads pileup to examine
     * @param ref reference base overlapping the pileup locus
     * @param features features overlapping the pileup locus
     * @return probability between 0.0 and 1.0 that the site is active (in practice with this implementation: either 0.0 or 1.0)
     */
    @Override
    public ActivityProfileState isActive( final AlignmentContext context, final ReferenceContext ref, final FeatureContext features ) {

        if ( hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vcFromAllelesRod = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(features, ref.getInterval(), false, logger, hcArgs.alleles);
            if( vcFromAllelesRod != null ) {
                return new ActivityProfileState(ref.getInterval(), 1.0);
            }
        }

        if ( hcArgs.USE_ALLELES_TRIGGER ) {
            return new ActivityProfileState( ref.getInterval(), features.getValues(hcArgs.alleles, ref.getInterval()).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null || context.getBasePileup().isEmpty() ) {
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getInterval(), 0.0);
        }

        final int ploidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied

        final Map<String, AlignmentContext> splitContexts;
        if (samplesList.numberOfSamples() == 1) {
            //If we know a priori that there's just one sample, take a shortcut and dont examine each read in the pileup
            splitContexts = context.splitContextBySampleName(samplesList.getSample(0), readsHeader);
        } else {
            splitContexts = context.splitContextBySampleName(readsHeader);
        }

        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        for( final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet() ) {
            // The ploidy here is not dictated by the sample but by the simple genotyping-engine used to determine whether regions are active or not.
            final int activeRegionDetectionHackishSamplePloidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
            final double[] genotypeLikelihoods = referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(activeRegionDetectionHackishSamplePloidy,sample.getValue().getBasePileup(), ref.getBase(), hcArgs.minBaseQualityScore, averageHQSoftClips).getGenotypeLikelihoods();
            genotypes.add( new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make() );
        }

        final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE , FAKE_ALT_ALLELE);
        final double isActiveProb;

        if (genotypes.size() == 1) {
            // Faster implementation avoiding the costly and over complicated Exact AFCalculator machinery:
            // This is the case when doing GVCF output.
            isActiveProb = activeRegionEvaluationGenotyperEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(genotypes.get(0).getLikelihoods().getAsVector());
        } else {
            final VariantCallContext vcOut = activeRegionEvaluationGenotyperEngine.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getEnd(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.SNP, readsHeader);
            isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb(vcOut.getPhredScaledQual());
        }
        return new ActivityProfileState(ref.getInterval(), isActiveProb, averageHQSoftClips.mean() > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean() );
    }

    /**
     * Generate variant calls for an assembly region
     *
     * @param region region to assemble and perform variant calling on
     * @param features Features overlapping the assembly region
     * @return List of variants discovered in the region (may be empty)
     */
    public List<VariantContext> callRegion( final AssemblyRegion region, final FeatureContext features ) {
        if ( hcArgs.justDetermineActiveRegions ) {
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;
        }

        if ( hcArgs.sampleNameToUse != null ) {
            removeReadsFromAllSamplesExcept(hcArgs.sampleNameToUse, region);
        }

        if( ! region.isActive() ) {
            // Not active so nothing to do!
            return referenceModelForNoVariation(region, true);
        }

        final List<VariantContext> givenAlleles = new ArrayList<>();
        if ( hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            features.getValues(hcArgs.alleles).stream().filter(VariantContext::isNotFiltered).forEach(givenAlleles::add);

            // No alleles found in this region so nothing to do!
            if ( givenAlleles.isEmpty() ) {
                return referenceModelForNoVariation(region, true);
            }
        } else if( region.size() == 0 ) {
            // No reads here so nothing to do!
            return referenceModelForNoVariation(region, true);
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult =  AssemblyBasedCallerUtils.assembleReads(region, givenAlleles, hcArgs, readsHeader, samplesList, logger, referenceReader, assemblyEngine);

        final SortedSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        // TODO - line bellow might be unnecessary : it might be that assemblyResult will always have those alleles anyway
        // TODO - so check and remove if that is the case:
        allVariationEvents.addAll(givenAlleles);

        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(region, allVariationEvents);

        if ( ! trimmingResult.isVariationPresent() && ! hcArgs.disableOptimizations ) {
            return referenceModelForNoVariation(region, false);
        }

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);
        final Map<String, List<GATKRead>> perSampleFilteredReadList = splitReadsBySample(filteredReads);

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( ! assemblyResult.isVariationPresent() && ! hcArgs.disableOptimizations ) {
            return referenceModelForNoVariation(region, false);
        }

        // For sure this is not true if gVCF is on.
        if ( hcArgs.dontGenotype ) {
            return NO_CALLS; // user requested we not proceed
        }

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if ( regionForGenotyping.size() == 0 && ! hcArgs.disableOptimizations ) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(region, false);
        }

        // evaluate each sample's reads against all haplotypes
        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKRead>> reads = splitReadsBySample(regionForGenotyping.getReads());

        // Calculate the likelihoods: CPU intensive part.
        final ReadLikelihoods<Haplotype> readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads);

        // Realign reads to their best haplotype.
        final Map<GATKRead, GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                readLikelihoods,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getSpan(),
                features,
                (hcArgs.assemblerArgs.consensusMode ? Collections.<VariantContext>emptyList() : givenAlleles),
                emitReferenceConfidence(),
                readsHeader);

        if ( haplotypeBAMWriter.isPresent() ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if ( hcArgs.disableOptimizations ) {
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            }
            haplotypeBAMWriter.get().writeReadsAlignedToHaplotypes(haplotypes, assemblyResult.getPaddedReferenceLoc(), haplotypes,
                                                             calledHaplotypeSet, readLikelihoods);
        }

        if( hcArgs.debug) {
            logger.info("----------------------------------------------------------------------------------");
        }

        if ( emitReferenceConfidence() ) {
            if ( !containsCalls(calledHaplotypes) ) {
                // no called all of the potential haplotypes
                return referenceModelForNoVariation(region, false);
            }
            else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section:
                if (trimmingResult.hasLeftFlankingRegion()) {
                    result.addAll(referenceModelForNoVariation(trimmingResult.nonVariantLeftFlankRegion(), false));
                }
                // output variant containing region.
                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        readLikelihoods, genotypingEngine.getPloidyModel(), calledHaplotypes.getCalls()));
                // output right-flanking non-variant section:
                if (trimmingResult.hasRightFlankingRegion()) {
                    result.addAll(referenceModelForNoVariation(trimmingResult.nonVariantRightFlankRegion(), false));
                }
                return result;
            }
        }
        else {
            return calledHaplotypes.getCalls();
        }
    }

    private boolean containsCalls(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        return calledHaplotypes.getCalls().stream()
                .flatMap(call -> call.getGenotypes().stream())
                .anyMatch(Genotype::isCalled);
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    private List<VariantContext> referenceModelForNoVariation(final AssemblyRegion region, final boolean needsToBeFinalized) {
        if ( emitReferenceConfidence() ) {
            //TODO - why the activeRegion cannot manage its own one-time finalization and filtering?
            //TODO - perhaps we can remove the last parameter of this method and the three lines bellow?
            if ( needsToBeFinalized ) {
                finalizeRegion(region);
            }
            filterNonPassingReads(region);

            final SimpleInterval paddedLoc = region.getExtendedSpan();
            final Haplotype refHaplotype = AssemblyBasedCallerUtils.createReferenceHaplotype(region, paddedLoc, referenceReader);
            final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
            return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                    paddedLoc, region, createDummyStratifiedReadMap(refHaplotype, samplesList, region),
                    genotypingEngine.getPloidyModel(), Collections.emptyList());
        }
        else {
            return NO_CALLS;
        }
    }

    /**
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     * @param refHaplotype a non-null reference haplotype
     * @param samples a list of all samples
     * @param region the assembly region containing reads
     * @return a map from sample -> PerReadAlleleLikelihoodMap that maps each read to ref
     */
    public ReadLikelihoods<Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                   final SampleList samples,
                                                                   final AssemblyRegion region) {
        return new ReadLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                                     splitReadsBySample(samples, region.getReads()));
    }

    /**
     * Shutdown this HC engine, closing resources as appropriate
     */
    public void shutdown() {
        likelihoodCalculationEngine.close();

        if ( haplotypeBAMWriter.isPresent() ) {
            haplotypeBAMWriter.get().close();
        }
    }

    private void finalizeRegion(final AssemblyRegion region) {
        AssemblyBasedCallerUtils.finalizeRegion(region, hcArgs.errorCorrectReads, hcArgs.dontUseSoftClippedBases, minTailQuality, readsHeader, samplesList);
    }

    private Set<GATKRead> filterNonPassingReads( final AssemblyRegion activeRegion ) {
        // TODO: can we unify this additional filtering with makeStandardHCReadFilter()?

        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : activeRegion.getReads() ) {
            if( rec.getLength() < READ_LENGTH_FILTER_THRESHOLD || rec.getMappingQuality() < READ_QUALITY_FILTER_THRESHOLD || ! ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(rec) || (hcArgs.keepRG != null && !rec.getReadGroup().equals(hcArgs.keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    private Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return splitReadsBySample(samplesList, reads);
    }

    private Map<String, List<GATKRead>> splitReadsBySample( final SampleList samplesList, final Collection<GATKRead> reads ) {
        return AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, reads);
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    private void removeReadsFromAllSamplesExcept(final String targetSample, final AssemblyRegion activeRegion) {
        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : activeRegion.getReads() ) {
            if( ! ReadUtils.getSampleName(rec, readsHeader).equals(targetSample) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
    }
}
