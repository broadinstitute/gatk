package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
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
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.fragments.FragmentUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;

import java.io.File;
import java.io.FileNotFoundException;
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

    // path to a reference file
    private final String reference;

    private ReferenceConfidenceModel referenceConfidenceModel = null;

    private double log10GlobalReadMismappingRate;

    private AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer();

    // the genotyping engine for the isActive() determination
    private MinimalGenotypingEngine activeRegionEvaluationGenotyperEngine = null;

    private ReadThreadingAssembler assemblyEngine = null;

    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    private HaplotypeCallerGenotypingEngine genotypingEngine = null;

    private VariantAnnotatorEngine annotationEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    protected CachingIndexedFastaSequenceFile referenceReader;

    // writes Haplotypes to a bam file when the -bamout option is specified
    private HaplotypeBAMWriter haplotypeBAMWriter;

    private Set<String> sampleSet;
    private SampleList samplesList;

    // reference base padding size
    private static final int REFERENCE_PADDING = 500;

    private byte minTailQuality;
    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

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
    private final static int READ_LENGTH_FILTER_THRESHOLD = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;

    private static final List<VariantContext> NO_CALLS = Collections.emptyList();

    private final static Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private final static Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     *
     * @param hcArgs command-line arguments for the HaplotypeCaller
     * @param readsHeader header for the reads
     * @param reference path to the reference
     */
    public HaplotypeCallerEngine( final HaplotypeCallerArgumentCollection hcArgs, final SAMFileHeader readsHeader, final String reference ) {
        Utils.nonNull(hcArgs);
        Utils.nonNull(readsHeader);
        Utils.nonNull(reference);

        this.hcArgs = hcArgs;
        this.readsHeader = readsHeader;
        this.reference = reference;

        initialize();
    }

    private void initialize() {
        // Note: order of operations matters here!

        initializeSamples();

        // Must be called after initializeSamples()
        validateAndInitializeArgs();
        minTailQuality = (byte)(hcArgs.MIN_BASE_QUALTY_SCORE - 1);

        initializeActiveRegionEvaluationGenotyperEngine();

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(hcArgs.annotationGroupsToUse, hcArgs.annotationsToUse, hcArgs.annotationsToExclude, hcArgs.dbsnp.dbsnp, hcArgs.comps);

        genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, FixedAFCalculatorProvider.createThreadSafeProvider(hcArgs), ! hcArgs.doNotRunPhysicalPhasing);
        genotypingEngine.setAnnotationEngine(annotationEngine);

        referenceConfidenceModel = new ReferenceConfidenceModel(samplesList, readsHeader, hcArgs.indelSizeToEliminateInRefModel);

        //Allele-specific annotations are not yet supported in the VCF mode
        if (isAlleleSpecficMode(annotationEngine) && isVCFMode()){
           throw new UserException("Allele-specific annotations are not yet supported in the VCF mode");
        }
        initializeReferenceReader();

        initializeHaplotypeBamWriter();

        initializeAssemblyEngine();

        // create our likelihood calculation engine -- must be done after setting log10GlobalReadMismappingRate in validateAndInitializeArgs()
        likelihoodCalculationEngine = createLikelihoodCalculationEngine();

        trimmer.initialize(hcArgs.assemblyRegionTrimmerArgs, readsHeader.getSequenceDictionary(), hcArgs.DEBUG,
                hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, emitReferenceConfidence());
    }

    private boolean isVCFMode() {
        return hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.NONE;
    }

    private boolean isAlleleSpecficMode(final VariantAnnotatorEngine annotationEngine) {
        //HACK. Note: we can't use subclass information from ReducibleAnnotation (which would be the obvious choice)
        // because RMSMappingQuality is both a reducible annotation and a standard annotation.
        final boolean infoAnnAlleleSpefic = annotationEngine.getInfoAnnotations().stream().anyMatch(infoFieldAnnotation -> infoFieldAnnotation.getClass().getSimpleName().startsWith("AS_"));
        if (infoAnnAlleleSpefic){
            return true;
        }
        final boolean genoAnnAlleleSpefic = annotationEngine.getGenotypeAnnotations().stream().anyMatch(genotypeAnnotation -> genotypeAnnotation.getClass().getSimpleName().startsWith("AS_"));
        return genoAnnAlleleSpefic;
    }

    private void validateAndInitializeArgs() {
        if ( hcArgs.genotypeArgs.samplePloidy != HomoSapiensConstants.DEFAULT_PLOIDY && ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Currently, physical phasing is not available when ploidy is different than " + HomoSapiensConstants.DEFAULT_PLOIDY + "; therefore it won't be performed");
        }

        if ( hcArgs.dontGenotype && emitReferenceConfidence() ) {
            throw new UserException("You cannot request gVCF output and 'do not genotype' at the same time");
        }

        if ( emitReferenceConfidence() ) {
            if ( hcArgs.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) {
                throw new UserException.BadArgumentValue("ERC/gt_mode", "you cannot request reference confidence output and GENOTYPE_GIVEN_ALLELES at the same time");
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

        if ( hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate < 0 ) {
            hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate = -1;
        }

        // configure the global mismapping rate
        if ( hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate < 0 ) {
            log10GlobalReadMismappingRate = - Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate);
            logger.info("Using global mismapping rate of " + hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate + " => " + log10GlobalReadMismappingRate + " in log10 likelihood units");
        }

        if ( emitReferenceConfidence() && samplesList.numberOfSamples() != 1 ) {
            throw new UserException.BadArgumentValue("--emitRefConfidence", "Can only be used in single sample mode currently. Use the sample_name argument to run on a single sample out of a multi-sample BAM file.");
        }
    }

    private void initializeSamples() {
        sampleSet = ReadUtils.getSamplesFromHeader(readsHeader);
        samplesList = new IndexedSampleList(sampleSet);

        if ( hcArgs.sampleNameToUse != null ) {
            if ( ! sampleSet.contains(hcArgs.sampleNameToUse) ) {
                throw new UserException.BadArgumentValue("--sample_name", "Specified name does not exist in input bam files");
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

    private void initializeReferenceReader() {
        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(new File(reference));
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(new File(reference), e);
        }
    }

    private void initializeHaplotypeBamWriter() {
        if ( hcArgs.bamOutputPath != null ) {
            haplotypeBAMWriter = HaplotypeBAMWriter.create(hcArgs.bamWriterType, new File(hcArgs.bamOutputPath), readsHeader);
        }
    }

    private void initializeAssemblyEngine() {
        // create and setup the assembler
        assemblyEngine = new ReadThreadingAssembler(hcArgs.assemblerArgs.maxNumHaplotypesInPopulation, hcArgs.assemblerArgs.kmerSizes, hcArgs.assemblerArgs.dontIncreaseKmerSizesForCycles, hcArgs.assemblerArgs.allowNonUniqueKmersInRef, hcArgs.assemblerArgs.numPruningSamples);

        assemblyEngine.setErrorCorrectKmers(hcArgs.assemblerArgs.errorCorrectKmers);
        assemblyEngine.setPruneFactor(hcArgs.assemblerArgs.MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(hcArgs.DEBUG);
        assemblyEngine.setDebugGraphTransformations(hcArgs.assemblerArgs.debugGraphTransformations);
        assemblyEngine.setRecoverDanglingBranches(!hcArgs.assemblerArgs.doNotRecoverDanglingBranches);
        assemblyEngine.setMinDanglingBranchLength(hcArgs.assemblerArgs.minDanglingBranchLength);
        assemblyEngine.setMinBaseQualityToUseInAssembly(hcArgs.MIN_BASE_QUALTY_SCORE);

        if ( hcArgs.assemblerArgs.graphOutput != null ) {
            assemblyEngine.setGraphWriter(new File(hcArgs.assemblerArgs.graphOutput));
        }
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        switch ( hcArgs.likelihoodEngineImplementation) {
            case PairHMM:
                return new PairHMMLikelihoodCalculationEngine((byte) hcArgs.likelihoodArgs.gcpHMM, hcArgs.likelihoodArgs.pairHMM, log10GlobalReadMismappingRate, hcArgs.pcrErrorModel);
            case Random:
                return new RandomLikelihoodCalculationEngine();
            default:
                //Note: we do not include in the error message list as it is of no grand public interest.
                throw new UserException("Unsupported likelihood calculation engine: " + likelihoodCalculationEngine);
        }
    }

    /**
     * Create the default read filter for use with the HaplotypeCaller
     *
     * @param hcArgs HaplotypeCaller arguments
     * @param header reads header
     * @return the default read filter for use with the HaplotypeCaller
     */
    public static CountingReadFilter makeStandardHCReadFilter( final HaplotypeCallerArgumentCollection hcArgs, final SAMFileHeader header ) {
        Utils.nonNull(hcArgs);
        Utils.nonNull(header);

        return new CountingReadFilter("MAPPING_QUALITY", new MappingQualityReadFilter(hcArgs.MIN_MAPPING_QUALITY_SCORE))
                .and(new CountingReadFilter("MAPPING_QUALITY_AVAILABLE", ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE))
                .and(new CountingReadFilter("MAPPED", ReadFilterLibrary.MAPPED))
                .and(new CountingReadFilter("PRIMARY_ALIGNMENT", ReadFilterLibrary.PRIMARY_ALIGNMENT))
                .and(new CountingReadFilter("NOT_DUPLICATE", ReadFilterLibrary.NOT_DUPLICATE))
                .and(new CountingReadFilter("PASSES_VENDOR_QUALITY_CHECK", ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK))
                .and(new CountingReadFilter("GOOD_CIGAR", ReadFilterLibrary.GOOD_CIGAR))
                .and(new CountingReadFilter("WELLFORMED", new WellformedReadFilter(header)));

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
                throw new UserException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }
        }

        return writer;
    }

    /**
     * Writes an appropriate VCF header, given our arguments, to the provided writer
     *
     * @param vcfWriter writer to which the header should be written
     */
    public void writeHeader( final VariantContextWriter vcfWriter ) {
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

        vcfWriter.writeHeader(new VCFHeader(headerInfo, sampleSet));
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
            final double[] genotypeLikelihoods = referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(activeRegionDetectionHackishSamplePloidy,sample.getValue().getBasePileup(), ref.getBase(), hcArgs.MIN_BASE_QUALTY_SCORE, averageHQSoftClips).getGenotypeLikelihoods();
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
            for ( final VariantContext vc : features.getValues(hcArgs.alleles) ) {
                if ( vc.isNotFiltered() ) {
                    givenAlleles.add(vc); // do something with these VCs during GGA mode
                }
            }

            // No alleles found in this region so nothing to do!
            if ( givenAlleles.isEmpty() ) {
                return referenceModelForNoVariation(region, true);
            }
        } else if( region.size() == 0 ) {
            // No reads here so nothing to do!
            return referenceModelForNoVariation(region, true);
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(region, givenAlleles);

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
        final Map<GATKRead, GATKRead> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());
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

        if ( haplotypeBAMWriter != null ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if ( hcArgs.disableOptimizations ) {
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            }
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(haplotypes, assemblyResult.getPaddedReferenceLoc(), haplotypes,
                                                             calledHaplotypeSet, readLikelihoods);
        }

        if( hcArgs.DEBUG ) {
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

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    private Map<GATKRead, GATKRead> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final Locatable paddedReferenceLoc) {
        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKRead, GATKRead> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKRead originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKRead realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead, realignedRead);
        }
        return result;
    }

    private boolean containsCalls(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        final List<VariantContext> calls = calledHaplotypes.getCalls();
        if (calls.isEmpty()) return false;
        for (final VariantContext call : calls)
            for (final Genotype genotype : call.getGenotypes())
                if (genotype.isCalled())
                    return true;
        return false;
    }

    /**
     * High-level function that runs the assembler on the given region's reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param region the region we should assemble
     * @param giveAlleles additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    private AssemblyResultSet assembleReads(final AssemblyRegion region, final List<VariantContext> giveAlleles) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeRegion(region); // handle overlapping fragments, clip adapter and low qual tails
        if( hcArgs.DEBUG ) {
            logger.info("Assembling " + region.getSpan() + " with " + region.size() + " reads:    (with overlap region = " + region.getExtendedSpan() + ")");
        }

        final byte[] fullReferenceWithPadding = region.getAssemblyRegionReference(referenceReader, REFERENCE_PADDING);
        final SimpleInterval paddedReferenceLoc = getPaddedReferenceLoc(region);
        final Haplotype referenceHaplotype = createReferenceHaplotype(region, paddedReferenceLoc);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        ReadErrorCorrector readErrorCorrector = null;
        if ( hcArgs.errorCorrectReads ) {
            readErrorCorrector = new ReadErrorCorrector(hcArgs.assemblerArgs.kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION, hcArgs.assemblerArgs.minObservationsForKmerToBeSolid, hcArgs.DEBUG, fullReferenceWithPadding);
        }

        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly(region, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, giveAlleles, readErrorCorrector, readsHeader);
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;
        } catch ( final Exception e ) {
            // Capture any exception that might be thrown, and write out the assembly failure BAM if requested
            if ( hcArgs.captureAssemblyFailureBAM ) {
                try ( final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File("assemblyFailure.bam"), null, readsHeader, false, false, false) ) {
                    for ( final GATKRead read : region.getReads() ) {
                        writer.addAlignment(read.convertToSAMRecord(readsHeader));
                    }
                }
            }
            throw e;
        }
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    private Haplotype createReferenceHaplotype(final AssemblyRegion region, final SimpleInterval paddedReferenceLoc) {
        return ReferenceConfidenceModel.createReferenceHaplotype(region, region.getAssemblyRegionReference(referenceReader), paddedReferenceLoc);
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
            final Haplotype refHaplotype = createReferenceHaplotype(region, paddedLoc);
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

        if ( haplotypeBAMWriter != null ) {
            haplotypeBAMWriter.close();
        }
    }

    private void finalizeRegion( final AssemblyRegion region ) {
        if ( region.isFinalized() ) {
            return;
        }

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKRead> readsToUse = new ArrayList<>(region.getReads().size());
        for( final GATKRead myRead : region.getReads() ) {
            GATKRead clippedRead;
            if ( hcArgs.errorCorrectReads ) {
                clippedRead = ReadClipper.hardClipLowQualEnds(myRead, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION);
            }
            else { // default case: clip low qual ends of reads
                clippedRead = ReadClipper.hardClipLowQualEnds(myRead, minTailQuality);
            }

            if ( hcArgs.dontUseSoftClippedBases || ! ReadUtils.hasWellDefinedFragmentSize(clippedRead) ) {
                // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
                clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
            }
            else {
                // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                // TODO -- reference haplotype start must be removed
                clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
            }

            clippedRead = clippedRead.isUnmapped() ? clippedRead : ReadClipper.hardClipAdaptorSequence(clippedRead);
            if ( ! clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0 ) {
                clippedRead = ReadClipper.hardClipToRegion( clippedRead, region.getExtendedSpan().getStart(), region.getExtendedSpan().getEnd() );
                if ( region.readOverlapsRegion(clippedRead) && clippedRead.getLength() > 0 ) {
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.
        // final List<GATKRead> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);
        Collections.sort(readsToUse, new ReadCoordinateComparator(readsHeader)); // TODO: sort may be unnecessary here

        // handle overlapping read pairs from the same fragment
        cleanOverlappingReadPairs(readsToUse);

        region.clearReads();
        region.addAll(readsToUse);
        region.setFinalized(true);
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

    private SimpleInterval getPaddedReferenceLoc( final AssemblyRegion region ) {
        final int padLeft = Math.max(region.getExtendedSpan().getStart() - REFERENCE_PADDING, 1);
        final int padRight = Math.min(region.getExtendedSpan().getEnd() + REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(region.getExtendedSpan().getContig()).getSequenceLength());
        return new SimpleInterval(region.getExtendedSpan().getContig(), padLeft, padRight);
    }

    private Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return splitReadsBySample(samplesList, reads);
    }

    private Map<String, List<GATKRead>> splitReadsBySample( final SampleList samplesList, final Collection<GATKRead> reads ) {
        final Map<String, List<GATKRead>> returnMap = new HashMap<>();
        final int sampleCount = samplesList.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            returnMap.put(samplesList.getSample(i), new ArrayList<>());
        }

        for ( final GATKRead read : reads ) {
            returnMap.get(ReadUtils.getSampleName(read, readsHeader)).add(read);
        }

        return returnMap;
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private void cleanOverlappingReadPairs(final List<GATKRead> reads) {
        for ( final List<GATKRead> perSampleReadList : splitReadsBySample(reads).values() ) {
            final FragmentCollection<GATKRead> fragmentCollection = FragmentCollection.create(perSampleReadList);
            for ( final List<GATKRead> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
            }
        }
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
