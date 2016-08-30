package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.MathArrays;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

import static java.lang.Math.pow;

/**
 * Call somatic SNPs and indels via local re-assembly of haplotypes
 *
 * <p>MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (<a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>) with the assembly-based machinery of HaplotypeCaller.</p>
 *
 * <p>The basic operation of MuTect2 proceeds similarly to that of the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>   </p>
 *
 * <h3>Differences from HaplotypeCaller</h3>
 * <p>While the HaplotypeCaller relies on a ploidy assumption (diploid by default) to inform its genotype likelihood and
 * variant quality calculations, MuTect2 allows for a varying allelic fraction for each variant, as is often seen in tumors with purity less
 * than 100%, multiple subclones, and/or copy number variation (either local or aneuploidy). MuTect2 also differs from the HaplotypeCaller in that it does apply some hard filters
 * to variants before producing output.</p>
 *
 * <h3>Usage examples</h3>
 * <p>These are example commands that show how to run MuTect2 for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the
 * Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Tumor/Normal variant calling</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar \
 *     -T MuTect2 \
 *     -R reference.fasta \
 *     -I tumor.bam \
 *     -I normal.bam \
 *     -tumor tumorSampleName \ // as in the BAM header
 *     -normal normalSampleName \ // as in the BAM header
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     [-L targets.interval_list] \
 *     -o output.vcf
 * </pre>
 *
 * <h4>Normal-only calling for panel of normals creation</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T MuTect2
 *     -R reference.fasta
 *     -I:tumor normal1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     --artifact_detection_mode \
 *     [-L targets.interval_list] \
 *     -o output.normal1.vcf
 * </pre>
 * <br />
 * For full PON creation, call each of your normals separately in artifact detection mode. Then use CombineVariants to
 * output only sites where a variant was seen in at least two samples:
 * <pre>
 * java -jar GenomeAnalysisTK.jar
 *     -T CombineVariants
 *     -R reference.fasta
 *     -V output.normal1.vcf -V output.normal2.vcf [-V output.normal2.vcf ...] \
 *     -minN 2 \
 *     --setKey "null" \
 *     --filteredAreUncalled \
 *     --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
 *     [-L targets.interval_list] \
 *     -o MuTect2_PON.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>MuTect2 currently only supports the calling of a single tumor-normal pair at a time</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Call somatic SNPs and indels via local re-assembly of haplotypes",
        oneLineSummary = "Call somatic SNPs and indels via local re-assembly of haplotypes",
        programGroup = VariantProgramGroup.class
)
public class MuTect2 extends AssemblyRegionWalker {
    private static final Logger logger = LogManager.getLogger(MuTect2.class);

    //TODO: move to GATKVCFConstants in gatk public
    public static final String TLOD_FWD_KEY =                       "TLOD_FWD";
    public static final String TLOD_REV_KEY =                       "TLOD_REV";
    public static final String TUMOR_SB_POWER_FWD_KEY =             "TUMOR_SB_POWER_FWD";
    public static final String TUMOR_SB_POWER_REV_KEY =             "TUMOR_SB_POWER_REV";
    public static final String TRIALLELIC_SITE_FILTER_NAME =                  "triallelic_site"; //M2
    public static final String STRAND_ARTIFACT_FILTER_NAME =                  "strand_artifact"; // M2
    public static final String CLUSTERED_READ_POSITION_FILTER_NAME =          "clustered_read_position"; // M2
    public static final String MEDIAN_LEFT_OFFSET_KEY =              "MEDIAN_LEFT_OFFSET";
    public static final String MEDIAN_RIGHT_OFFSET_KEY =             "MEDIAN_RIGHT_OFFSET";
    public static final String MAD_MEDIAN_LEFT_OFFSET_KEY =          "MAD_LEFT_OFFSET";
    public static final String MAD_MEDIAN_RIGHT_OFFSET_KEY =         "MAD_RIGHT_OFFSET";

    private SAMFileHeader header;
    protected SampleList samplesList;
    protected boolean printTCGAsampleHeader = false;
    private Set<String> normalReadGroupIDs;

    protected CachingIndexedFastaSequenceFile referenceReader;
    protected ReadThreadingAssembler assemblyEngine = null;
    protected ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;
    protected SomaticGenotypingEngine genotypingEngine = null;
    private HaplotypeBAMWriter haplotypeBAMWriter;
    final AssemblyRegionEvaluator mutect2AssemblyRegionEvaluator = new Mutect2AssemblyRegionEvaluator();

    private byte MIN_TAIL_QUALITY;
    private double log10GlobalReadMismappingRate;
    private static final int REFERENCE_PADDING = 500;
    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;
    private final static int maxReadsInRegionPerSample = 1000; // TODO -- should be an argument
    private final static int minReadsPerAlignmentStart = 5; // TODO -- should be an argument

    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @ArgumentCollection
    protected ReadThreadingAssemblerArgumentCollection RTAC = new ReadThreadingAssemblerArgumentCollection();

    @ArgumentCollection
    protected LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();

    @Argument(fullName = "debug_read_name", optional = true, doc="trace this read name through the calling process")
    protected String DEBUG_READ_NAME = null;

    @Hidden
    @Advanced
    @Argument(fullName = "MQ_filtering_level", shortName = "MQthreshold", optional = true, doc="Set an alternate MQ threshold for debugging")
    final public int MQthreshold = 20;

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     *  as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).
     *  Records that are filtered in the comp track will be ignored.
     *  Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Argument(fullName="comp", shortName = "comp", doc="comparison VCF file", optional = true)
    public List<FeatureInput<VariantContext>> comps = Collections.emptyList();

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     * //TODO: port TandemRepeatAnnotator and put it here
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", optional = true)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"DepthPerAlleleBySample", "BaseQualitySumPerAlleleBySample", "OxoGReadCounts"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Advanced
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", optional = true)
    protected List<String> annotationsToExclude = new ArrayList<>(Arrays.asList(new String[]{}));

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", optional = true)
    public List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{"StandardAnnotation", "StandardHCAnnotation"}));

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public String outputVCF = null;



    //TODO: HACK ALERT HACK ALERT HACK ALERT
    //TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam -I:normal normal.bam,
    //TODO: so for now we require the user to specify bams *both* as inputs, with -I tumor.bam -I normal.bam
    //TODO: *and* as sample names e.g. -tumor tumorSampleName -normal normalSampleName

    @Argument(fullName = "tumorSampleName", shortName = "tumor", doc = "BAM sample name of tumor", optional = false)
    protected String tumorSampleName = null;

    @Argument(fullName = "normalSampleName", shortName = "normal", doc = "BAM sample name of tumor", optional = true)
    protected String normalSampleName = null;

    //TODO: END OF HACK ALERT



    private VariantContextWriter vcfWriter;

    protected AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer();

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use read from this read group when making calls (but use all reads to build the assembly)", optional = true)
    protected String keepRG = null;

    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", optional = true)
    public byte MIN_BASE_QUALTY_SCORE = 10;

    public PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCRErrorModel.HOSTILE;

    @Hidden
    @Argument(fullName="errorCorrectReads", shortName="errorCorrectReads", doc = "Use an exploratory algorithm to error correct the kmers used during assembly.  May cause fundamental problems with the assembly graph itself", optional = true)
    protected boolean errorCorrectReads = false;

    @Hidden
    @Argument(fullName="captureAssemblyFailureBAM", shortName="captureAssemblyFailureBAM", doc="If specified, we will write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", optional = true)
    protected boolean captureAssemblyFailureBAM = false;

    @Advanced
    @Argument(fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", doc="If specified, we will not analyze soft clipped bases in the reads", optional = true)
    protected boolean dontUseSoftClippedBases = false;

    @Hidden
    @Argument(fullName="justDetermineAssemblyRegions", shortName="justDetermineAssemblyRegions", doc = "If specified, the HC won't actually do any assembly or calling, it'll just run the upfront active region determination code.  Useful for benchmarking and scalability testing", optional = true)
    protected boolean justDetermineAssemblyRegions = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName="doNotRunPhysicalPhasing", shortName="doNotRunPhysicalPhasing", doc="Disable physical phasing", optional = true)
    protected boolean doNotRunPhysicalPhasing = false;

    @Override
    protected int defaultReadShardSize() { return 5000; }

    @Override
    protected int defaultReadShardPadding() { return 100; }

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    public void onTraversalStart() {
        header = getHeaderForReads();
        vcfWriter = GATKVariantContextUtils.createVCFWriter(new File(outputVCF), header.getSequenceDictionary(), false);

        normalReadGroupIDs = header.getReadGroups().stream()
                .filter(rg -> rg.getSample().equals(normalSampleName))
                .map(rg -> rg.getId())
                .collect(Collectors.toSet());
        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(header)));
        if (!samplesList.asListOfSamples().contains(tumorSampleName)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given tumor" +
            " sample name " + tumorSampleName);
        } else if (normalSampleName != null && !samplesList.asListOfSamples().contains(normalSampleName)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given normal" +
                    " sample name " + normalSampleName);
        }

        //If the samples specified are exactly one normal and one tumor, use the TCGA VCF sample header format
        printTCGAsampleHeader = samplesList.numberOfSamples() == 2 && normalSampleName != null;

        final VariantAnnotatorEngine annotationEngine = initializeVCFOutput();

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(referenceArguments.getReferenceFile());
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(referenceArguments.getReferenceFile(), e);
        }

        // create and setup the assembler
        assemblyEngine = new ReadThreadingAssembler(RTAC.maxNumHaplotypesInPopulation, RTAC.kmerSizes, RTAC.dontIncreaseKmerSizesForCycles, RTAC.allowNonUniqueKmersInRef, RTAC.numPruningSamples);

        assemblyEngine.setErrorCorrectKmers(RTAC.errorCorrectKmers);
        assemblyEngine.setPruneFactor(RTAC.MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(MTAC.DEBUG);
        assemblyEngine.setDebugGraphTransformations(RTAC.debugGraphTransformations);
        assemblyEngine.setRecoverDanglingBranches(!RTAC.doNotRecoverDanglingBranches);
        assemblyEngine.setMinBaseQualityToUseInAssembly(MIN_BASE_QUALTY_SCORE);

        MIN_TAIL_QUALITY = (byte)(MIN_BASE_QUALTY_SCORE - 1);

        if ( RTAC.graphOutput != null ) assemblyEngine.setGraphWriter(new File(RTAC.graphOutput));

        // setup the likelihood calculation engine
        if ( LEAC.phredScaledGlobalReadMismappingRate < 0 ) LEAC.phredScaledGlobalReadMismappingRate = -1;

        // configure the global mismapping rate
        if ( LEAC.phredScaledGlobalReadMismappingRate < 0 ) {
            log10GlobalReadMismappingRate = - Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(LEAC.phredScaledGlobalReadMismappingRate);
            logger.info("Using global mismapping rate of " + LEAC.phredScaledGlobalReadMismappingRate + " => " + log10GlobalReadMismappingRate + " in log10 likelihood units");
        }

        likelihoodCalculationEngine = createLikelihoodCalculationEngine();

        genotypingEngine = new SomaticGenotypingEngine(samplesList, !doNotRunPhysicalPhasing, MTAC,
                tumorSampleName, normalSampleName, DEBUG_READ_NAME);
        genotypingEngine.setAnnotationEngine(annotationEngine);


        if ( MTAC.bamOutputPath != null ) {
            haplotypeBAMWriter = HaplotypeBAMWriter.create(MTAC.bamWriterType, new File(MTAC.bamOutputPath), header);
        }

        // why isn't this a constructor (instead of initialize)?  Since the method is package-friendly
        trimmer.initialize(MTAC.assemblyRegionTrimmerArgs, header.getSequenceDictionary(), MTAC.DEBUG,
                MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, false);
        final Set<String> sampleSet = samplesList.asSetOfSamples();

        if( MTAC.CONTAMINATION_FRACTION_FILE != null )
            MTAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(MTAC.CONTAMINATION_FRACTION_FILE, MTAC.CONTAMINATION_FRACTION, sampleSet, logger));

    }

    private VariantAnnotatorEngine initializeVCFOutput() {
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            annotationsToUse.add("ClusteredReadPosition");
        }
        final VariantAnnotatorEngine annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse,
                annotationsToUse,
                annotationsToExclude,
                MTAC.dbsnp.dbsnp,
                comps);

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

        //TODO: see comment in getM2HeaderLines about some M2 filters undefined in GATKVCFHEaderLines
        headerInfo.addAll(getM2HeaderLines());
        headerInfo.addAll(getSampleHeaderLines());

        // if printTCGAsampleHeader, we already checked for exactly 1 tumor and 1 normal in printTCGAsampleHeader assignment in initialize()
        final List<String> outputSampleNames = printTCGAsampleHeader ? Arrays.asList("TUMOR", "NORMAL") : samplesList.asListOfSamples();

        //TODO: this hack removes nulls due to the aforementioned problem with M2 header lines
        final Set<VCFHeaderLine> nonNullHEaderInfo = headerInfo.stream()
                .filter(line -> line != null).collect(Collectors.toSet());
        vcfWriter.writeHeader(new VCFHeader(nonNullHEaderInfo, outputSampleNames));

        return annotationEngine;
    }

    private Set<VCFHeaderLine> getM2HeaderLines(){
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NORMAL_LOD_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TUMOR_LOD_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.PANEL_OF_NORMALS_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.HAPLOTYPE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY));

        if (MTAC.ENABLE_STRAND_ARTIFACT_FILTER){
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(TLOD_FWD_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(TLOD_REV_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(TUMOR_SB_POWER_FWD_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(TUMOR_SB_POWER_REV_KEY));
        }

        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ALLELE_FRACTION_KEY));

        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.PON_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.ALT_ALLELE_IN_NORMAL_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.TUMOR_LOD_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME));


        // TODO: GATKVCFHeaderLines has a giant static block defining VCFHeaderLines
        // TODO: the following M2 filters do not yet appear there in GATK4
        // TODO: thus getFilterLine returns null
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(TRIALLELIC_SITE_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(STRAND_ARTIFACT_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(CLUSTERED_READ_POSITION_FILTER_NAME));


        if ( ! doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }
        return headerInfo;
    }

    private Set<VCFHeaderLine> getSampleHeaderLines(){
        final Set<VCFHeaderLine> sampleLines = new HashSet<>();
        if (printTCGAsampleHeader) {
            //NOTE: This will only list the first bam file for each tumor/normal sample if there is more than one
            final Map<String, String> normalSampleHeaderAttributes = new HashMap<>();
            normalSampleHeaderAttributes.put("ID", "NORMAL");
            normalSampleHeaderAttributes.put("SampleName", normalSampleName);
            final VCFSimpleHeaderLine normalSampleHeader = new VCFSimpleHeaderLine("SAMPLE", normalSampleHeaderAttributes);
            final Map<String, String> tumorSampleHeaderAttributes = new HashMap<>();
            tumorSampleHeaderAttributes.put("ID", "TUMOR");
            tumorSampleHeaderAttributes.put("SampleName", tumorSampleName);
            final VCFSimpleHeaderLine tumorSampleHeader = new VCFSimpleHeaderLine("SAMPLE", tumorSampleHeaderAttributes);

            sampleLines.add(normalSampleHeader);
            sampleLines.add(tumorSampleHeader);
        }
        return sampleLines;
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() { return mutect2AssemblyRegionEvaluator; }

    private final static List<VariantContext> NO_CALLS = Collections.emptyList();

    public List<VariantContext> callRegion(final AssemblyRegion originalAssemblyRegion, final FeatureContext featureContext ) {
        if ( justDetermineAssemblyRegions ) {
            return NO_CALLS;
        } else if( !originalAssemblyRegion.isActive() ) {
            return referenceModelForNoVariation(originalAssemblyRegion, true);
        } else if( originalAssemblyRegion.size() == 0 ) {
            return referenceModelForNoVariation(originalAssemblyRegion, true);
        }
        logReadInfo(DEBUG_READ_NAME, originalAssemblyRegion.getReads(), "Present in original active region");

        // create the assembly using just high quality reads (Q20 or higher).  We want to use lower
        // quality reads in the PairHMM (and especially in the normal) later, so we can't use a ReadFilter
        final AssemblyRegion assemblyAssemblyRegion = new AssemblyRegion(originalAssemblyRegion.getSpan(), originalAssemblyRegion.getSupportingStates(), originalAssemblyRegion.isActive(), originalAssemblyRegion.getExtension(), header);
        for (final GATKRead rec : originalAssemblyRegion.getReads()) {
            if (rec.getMappingQuality() >= MQthreshold ) {
                assemblyAssemblyRegion.add(rec);
            }
        }
        logReadInfo(DEBUG_READ_NAME, assemblyAssemblyRegion.getReads(), "Present in assembly active region");

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(assemblyAssemblyRegion, Collections.emptyList());
        final SortedSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(originalAssemblyRegion,allVariationEvents);
        if (!trimmingResult.isVariationPresent()) {
            return referenceModelForNoVariation(originalAssemblyRegion, false);
        }
        logReadInfo(DEBUG_READ_NAME, trimmingResult.getCallableRegion().getReads(), "Present in trimming result");

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        // it is conceivable that the region is active yet has no events upon assembling only the well-mapped reads
        if( ! assemblyResult.isVariationPresent() ) {
            return referenceModelForNoVariation(originalAssemblyRegion, false);
        }

        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping");

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalAssemblyRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);

        final Map<String, List<GATKRead>> perSampleFilteredReadList = splitReadsBySample(filteredReads);
        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping after filtering reads");

        // evaluate each sample's reads against all haplotypes

        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKRead>> reads = splitReadsBySample( regionForGenotyping.getReads() );
        for (final List<GATKRead> rec : reads.values()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after splitting assemblyResult by sample");
        }

        //TODO: this should be written as
        //TODO final Map<String, Integer> ARreads_origNormalMQ = regionForGenotyping.getReads().stream()
        //TODO        .collect(Collectors.toMap(GATKRead::getName, GATKRead::getMappingQuality));
        //TODO: but for some reason sometimes streamifying Mutect2 code causes a MalformedWalkerArgumentsException
        //TODO: this will probably not occur after the port to GATK4
        final HashMap<String, Integer> ARreads_origNormalMQ = new HashMap<>();
        for (final GATKRead read : regionForGenotyping.getReads()) {
            ARreads_origNormalMQ.put(read.getName(), read.getMappingQuality());
        }

        // modify MAPQ scores in normal to be high so that we don't do any base quality score capping
        for(final GATKRead rec : regionForGenotyping.getReads()) {
            if (isReadFromNormal(rec)) {
                rec.setMappingQuality(60);
            }
        }

        logger.debug("Computing read likelihoods with " + regionForGenotyping.getReads().size() + " reads against " + haplotypes.size() + " haplotypes across region " + assemblyResult.getRegionForGenotyping().toString());
        final ReadLikelihoods<Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);

        final Map<GATKRead,GATKRead> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());

        readLikelihoods.changeReads(readRealignments);
        for (final GATKRead rec : readRealignments.keySet()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after computing read likelihoods");
        }

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods,
                ARreads_origNormalMQ,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getSpan(),
                featureContext,
                header);

        if ( MTAC.bamOutputPath != null ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (MTAC.disableOptimizations)
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    assemblyResult.getPaddedReferenceLoc(),
                    haplotypes,
                    calledHaplotypeSet,
                    readLikelihoods);
        }

        if( MTAC.DEBUG ) { logger.info("----------------------------------------------------------------------------------"); }
        return annotateVCs(calledHaplotypes, featureContext);
    }

    //TODO: this is duplicated from HaplotypeCaller
    @Override
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        final List<VariantContext> callsInRegion = callRegion(region, featureContext);

        for ( final VariantContext call : callsInRegion ) {
            // Only include calls that start within the current read shard (as opposed to the padded regions around it).
            // This is critical to avoid duplicating events that span shard boundaries!
            if ( getCurrentReadShardBounds().contains(new SimpleInterval(call.getContig(), call.getStart(), call.getStart())) ) {
                vcfWriter.add(call);
            }
        }
    }

    private Set<String> calculateFilters(final FeatureContext featureContext, final VariantContext vc, final Map<String, Object> eventDistanceAttributes) {
        final Set<String> filters = new HashSet<>();

        final Integer eventCount = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY);
        final Integer maxEventDistance = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY);

        //TODO: verify that vc.getStart() is the right second argument to use here
        final Collection<VariantContext> panelOfNormalsVC = featureContext.getValues(MTAC.normalPanelFeatureInput, vc.getStart());

        //TODO: is this the correct criterion for filtering?
        //TODO: Takutp points out that we're filtering based on the PoN even if the alleles are different!
        //TODO: and what about genotypes?
        if (!panelOfNormalsVC.isEmpty()) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }

        // FIXME: how do we sum qscores here?
        // FIXME: parameterize thresholds
        // && sum of alt likelihood scores > 20

        // TODO: make the change to have only a single normal sample (but multiple tumors is ok...)
        int normalAltCounts = 0;
        double normalF = 0;
        int normalAltQualityScoreSum = 0;
        if (hasNormal()) {
            final Genotype normalGenotype = vc.getGenotype(normalSampleName);

            // NOTE: how do we get the non-ref depth here?
            normalAltCounts = normalGenotype.getAD()[1];
            normalF = (Double) normalGenotype.getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);

            final Object qss = normalGenotype.getExtendedAttribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY);
            if (qss != null) {
                normalAltQualityScoreSum = (Integer) ((Object[]) qss)[1];
            } else {
                logger.error("Null qss at " + vc.getStart());
            }
        }

        if ( (normalAltCounts > MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || normalF > MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION ) && normalAltQualityScoreSum > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
            filters.add(GATKVCFConstants.ALT_ALLELE_IN_NORMAL_FILTER_NAME);
        } else if ( eventCount > 1 && normalAltCounts >= 1) {
            filters.add(GATKVCFConstants.MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME);
        } else if (eventCount >= 3) {
            filters.add(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME);
        }

        // STR contractions, that is the deletion of one repeat unit of a short repeat (>1bp repeat unit)
        // such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we
        // hard filter them out by default
        if (vc.isIndel()) {
            final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                    .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
            final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa != null && rpa.length > 1 && ru.length() > 1) {
                final int refCount = rpa[0];
                final int altCount = rpa[1];

                if (refCount - altCount == 1) {
                    filters.add(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME);
                }
            }
        }

        // NOTE: what if there is a 3bp indel followed by a snp... we are comparing starts
        // so it would be thrown but it's really an adjacent event
        if ( eventCount >= 2 && maxEventDistance >= 3) {
            filters.add(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME);
        }

        // clustered read position filter
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER){
            final Double tumorFwdPosMedian = (Double) vc.getAttribute(MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMedian = (Double) vc.getAttribute(MEDIAN_RIGHT_OFFSET_KEY);
            final Double tumorFwdPosMAD = (Double) vc.getAttribute(MAD_MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMAD = (Double) vc.getAttribute(MAD_MEDIAN_RIGHT_OFFSET_KEY);
            //If the variant is near the read end (median threshold) and the positions are very similar (MAD threshold) then filter
            if ( (tumorFwdPosMedian != null && tumorFwdPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorFwdPosMAD != null && tumorFwdPosMAD <= MTAC.PIR_MAD_THRESHOLD) ||
                    (tumorRevPosMedian != null && tumorRevPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorRevPosMAD != null && tumorRevPosMAD <= MTAC.PIR_MAD_THRESHOLD))
                filters.add(CLUSTERED_READ_POSITION_FILTER_NAME);
        }
        // TODO: Move strand bias filter here

        return filters;
    }

    private final static byte REF_MODEL_DELETION_QUAL = (byte) 30;
    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @return genotype likelihoods of [AA,AB]
     */
    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadPileup pileup, final byte refBase, final byte minBaseQual, final double f) {
        final double[] genotypeLikelihoods = new double[2];
        final int AA = 0;
        final int AB=1;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {

                // TODO: why not use base qualities here?
                //double pobs = QualityUtils.qualToErrorProbLog10(qual);
                final double pobs = 1.0d - pow(10, (30 / -10.0));
                if( isNonRef(refBase, p)) {
                    genotypeLikelihoods[AB] += Math.log10(f*pobs + (1-f)*pobs/3.0d);
                    genotypeLikelihoods[AA] += Math.log10((1-pobs)/3);
                } else {
                    genotypeLikelihoods[AB] += Math.log10(f*(1-pobs)/3.0d + (1-f)*pobs);
                    genotypeLikelihoods[AA] += Math.log10(pobs);
                }
            }
        }

        return genotypeLikelihoods;
    }

    private boolean hasNormal() {
        return (normalSampleName != null);
    }

    //TODO: streamify in GATK4
    protected int getCountOfNonRefEvents(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        int i=0;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    i++;
                }
            }
        }
        return i;
    }

    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        final double f = calculateF(pileup, refBase, minBaseQual);
        return calcGenotypeLikelihoodsOfRefVsAny(pileup, refBase, minBaseQual, f);
    }

    private double calculateF(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        int refCount = 0, altCount = 0;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());

            // only consider deletions AND sites of sufficient quality
            if( p.isDeletion() || qual > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    altCount++;
                } else {
                    refCount++;
                }
            }
        }
        return (double) altCount / (refCount + altCount);
    }

    private boolean isNonRef(final byte refBase, final PileupElement p) {
        return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
    }

    int MIN_READ_LENGTH = 30; // private in superclass

    protected Set<GATKRead> filterNonPassingReads( final AssemblyRegion assemblyRegion) {
        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : assemblyRegion.getReads() ) {
            //TODO: Takuto points out that this is questionable.  Let's think hard abut it.
            // KCIBUL: only perform read quality filtering on tumor reads...
            if (isReadFromNormal(rec)) {
                if( rec.getLength() < MIN_READ_LENGTH ) {
                    readsToRemove.add(rec);
                }
            } else if( rec.getLength() < MIN_READ_LENGTH || rec.getMappingQuality() < MQthreshold || hasBadMate(rec) ||
                    (keepRG != null && !rec.getReadGroup().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        assemblyRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        return new PairHMMLikelihoodCalculationEngine( (byte)LEAC.gcpHMM, LEAC.pairHMM, log10GlobalReadMismappingRate, pcrErrorModel );
    }

    /**
     * FROM HC
     *
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    protected List<VariantContext> referenceModelForNoVariation(final AssemblyRegion region, final boolean needsToBeFinalized) {
        return NO_CALLS;
    }

    protected Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return HaplotypeCallerEngine.splitReadsBySample(samplesList, header, reads);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }

        if ( likelihoodCalculationEngine != null ) {
            likelihoodCalculationEngine.close();
        }
    }

    /**
     * High-level function that runs the assembler on the active region reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param assemblyRegion the region we should assemble
     * @param giveAlleles additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    protected AssemblyResultSet assembleReads(final AssemblyRegion assemblyRegion, final List<VariantContext> giveAlleles) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeAssemblyRegion(assemblyRegion); // handle overlapping fragments, clip adapter and low qual tails

        final byte[] fullReferenceWithPadding = assemblyRegion.getAssemblyRegionReference(referenceReader, REFERENCE_PADDING);
        final SimpleInterval paddedReferenceLoc = getPaddedLoc(assemblyRegion);
        final Haplotype referenceHaplotype = createReferenceHaplotype(assemblyRegion, paddedReferenceLoc);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        final ReadErrorCorrector readErrorCorrector = errorCorrectReads ? new ReadErrorCorrector(RTAC.kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION,
                RTAC.minObservationsForKmerToBeSolid, MTAC.DEBUG, fullReferenceWithPadding) : null;
        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly( assemblyRegion, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, giveAlleles,readErrorCorrector, header );
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;
        } catch ( final Exception e ) {
            //TODO: code duplication: copied from HaplotypeCallerEngine
            // Capture any exception that might be thrown, and write out the assembly failure BAM if requested
            if ( captureAssemblyFailureBAM ) {
                try ( final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File("assemblyFailure.bam"), null, header, false, false, false) ) {
                    for ( final GATKRead read : assemblyRegion.getReads() ) {
                        writer.addAlignment(read.convertToSAMRecord(header));
                    }
                }
            }
            throw e;
        }
    }

    //TODO: code duplication -- this is identical to code in HaplotypeCallerEngine
    private void finalizeAssemblyRegion( final AssemblyRegion assemblyRegion ) {
        if (assemblyRegion.isFinalized()) return;

        if( MTAC.DEBUG ) { logger.info("Assembling " + assemblyRegion.getSpan() + " with " + assemblyRegion.size() + " reads:    (with overlap region = " + assemblyRegion.getExtendedSpan() + ")"); }

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKRead> readsToUse = new ArrayList<>(assemblyRegion.getReads().size());
        for( final GATKRead myRead : assemblyRegion.getReads() ) {
            GATKRead clippedRead;
            if (errorCorrectReads)
                clippedRead = ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION );
            else  // default case: clip low qual ends of reads
                clippedRead= ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY );

            if ( dontUseSoftClippedBases || ! ReadUtils.hasWellDefinedFragmentSize(clippedRead) ) {
                // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
                clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
            } else {
                // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                // TODO -- reference haplotype start must be removed
                clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
            }

            clippedRead = ( clippedRead.isUnmapped()? clippedRead : ReadClipper.hardClipAdaptorSequence( clippedRead ) );
            if( !clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0 ) {
                clippedRead = ReadClipper.hardClipToRegion(clippedRead, assemblyRegion.getExtendedSpan().getStart(), assemblyRegion.getExtendedSpan().getEnd());
                if( assemblyRegion.readOverlapsRegion(clippedRead) && clippedRead.getLength() > 0 ) {
                    readsToUse.add(clippedRead);
                }
            }
        }

        Collections.sort(readsToUse, new ReadCoordinateComparator(header));
        //final List<GATKRead> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        assemblyRegion.clearReads();
        assemblyRegion.addAll(readsToUse);
        assemblyRegion.setFinalized(true);
    }

    private SimpleInterval getPaddedLoc( final AssemblyRegion assemblyRegion ) {
        final int padLeft = Math.max(assemblyRegion.getExtendedSpan().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(assemblyRegion.getExtendedSpan().getEnd()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(assemblyRegion.getExtendedSpan().getContig()).getSequenceLength());
        return new SimpleInterval(assemblyRegion.getContig(), padLeft, padRight);
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param assemblyRegion the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the SimpleInterval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    private Haplotype createReferenceHaplotype(final AssemblyRegion assemblyRegion, final SimpleInterval paddedReferenceLoc) {
        return ReferenceConfidenceModel.createReferenceHaplotype(assemblyRegion, assemblyRegion.getAssemblyRegionReference(referenceReader), paddedReferenceLoc);
    }

    public static void logReadInfo(final String readName, final Collection<GATKRead> records, final String message) {
        if (readName != null) {
            for (final GATKRead rec : records) {
                logReadInfo(readName, rec, message);
            }
        }
    }

    public static void logReadInfo(final String readName, final GATKRead rec, final String message) {
        if (readName != null && rec != null && readName.equals(rec.getName())) {
            logger.info("Found " + rec.toString() + " - " + message);
        }
    }

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    // TODO: migrate from HC -> HCUtils Class and share it!
    private Map<GATKRead,GATKRead> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final SimpleInterval paddedReferenceLoc) {

        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKRead,GATKRead> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKRead originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKRead realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead,realignedRead);
        }
        return result;
    }

    //TODO: is there a faster way to do this?
    private boolean isReadFromNormal(final GATKRead rec) {
        return normalSampleName != null && normalReadGroupIDs.contains(rec.getReadGroup());
    }

    final List<VariantContext> annotateVCs(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes, final FeatureContext featureContext) {
        final int eventCount = calledHaplotypes.getCalls().size();
        final Map<String, Object> eventDistanceAttributes = new HashMap<>();    //TODO: should be Map<String, Integer> -- see TODO below
        eventDistanceAttributes.put(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount);
        if (eventCount > 1) {
            final int lastPosition = calledHaplotypes.getCalls().get(0).getStart();
            final int[] eventDistances = new int[calledHaplotypes.getCalls().size() - 1];
            for (int n = 0; n < eventDistances.length; n++) {
                eventDistances[n] = Math.abs(calledHaplotypes.getCalls().get(n + 1).getStart() - lastPosition);
            }
            final int maxEventDistance = calledHaplotypes.getCalls().get(eventCount - 1).getStart() - calledHaplotypes.getCalls().get(0).getStart();
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY, MathUtils.arrayMin(eventDistances));
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY, maxEventDistance);
        } else {    //TODO: putting null is a hack -- we should remove this and update the integration test md5s
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY, null);
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY, null);
        }

        final List<VariantContext> annotatedCalls = new ArrayList<>();
        // can we do this with the Annotation classes instead?
        for (final VariantContext originalVC : calledHaplotypes.getCalls()) {
            final VariantContextBuilder vcb = new VariantContextBuilder(originalVC);

            final Map<String, Object> attributes = new HashMap<>(originalVC.getAttributes());
            attributes.putAll(eventDistanceAttributes);
            vcb.attributes(attributes);

            final Set<String> filters = new HashSet<>(originalVC.getFilters());

            final double tumorLod = originalVC.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, -1);
            if (tumorLod < MTAC.TUMOR_LOD_THRESHOLD) {
                filters.add(GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
            }

            // if we are in artifact detection mode, apply the thresholds for the LOD scores
            if (!MTAC.ARTIFACT_DETECTION_MODE) {
                filters.addAll(calculateFilters(featureContext, originalVC, eventDistanceAttributes));
            }

            vcb.filters(filters.isEmpty() ? VariantContext.PASSES_FILTERS : filters);

            if (printTCGAsampleHeader) {
                final Genotype tumorGenotype = new GenotypeBuilder(originalVC.getGenotype(tumorSampleName)).name("TUMOR").make();
                final Genotype normalGenotype = new GenotypeBuilder(originalVC.getGenotype(normalSampleName)).name("NORMAL").make();
                vcb.genotypes(Arrays.asList(tumorGenotype, normalGenotype));
            }
            annotatedCalls.add(vcb.make());
        }
        return annotatedCalls;
    }

    //TODO: See if we can share code with HaplotypeCaller's isActive
    private class Mutect2AssemblyRegionEvaluator implements  AssemblyRegionEvaluator {
        public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext featureContext) {
            if( context == null || context.getBasePileup().isEmpty() ) {
                return new ActivityProfileState(ref.getInterval(), 0.0);
            }

            final Map<String, AlignmentContext> splitContexts = context.splitContextBySampleName(header);
            final AlignmentContext tumorContext = splitContexts.get(tumorSampleName);
            final AlignmentContext normalContext = splitContexts.get(normalSampleName);

            // if there are no tumor reads... there is no activity!
            if (tumorContext == null) {
                return new ActivityProfileState(ref.getInterval(), 0);
            }

            final ReadPileup tumorPileup = tumorContext.getBasePileup().makeFilteredPileup(el -> el.getMappingQual() >= MQthreshold);
            final double[] tumorGLs = calcGenotypeLikelihoodsOfRefVsAny(tumorPileup, ref.getBase(), MIN_BASE_QUALTY_SCORE);
            final double tumorLod = tumorGLs[1] - tumorGLs[0];

            // NOTE: do I want to convert to a probability (or just keep this as a LOD score)

            // also at this point, we should throw out noisy sites (hence the nonRefInNormalCheck) but this is non-optimal
            double prob = 0;
            if (tumorLod > MTAC.INITIAL_TUMOR_LOD_THRESHOLD) {

                // TODO: should we even do this performance optimization?
                // in any case, we have to handle the case where there is no normal (and thus no normal context) which is
                // different than having a normal but having no reads (where we should not enter the active region)
                if (normalSampleName != null && normalContext != null) {
                    final int nonRefInNormal = getCountOfNonRefEvents(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE);

                    final double[] normalGLs = calcGenotypeLikelihoodsOfRefVsAny(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE, 0.5f);
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


    // If mates in a pair are mapping to different contigs, it is likely that at least one of them is in the wrong place.
    // Of course, it is also possible there is a structural variant, in which case we lose some sensitivity
    public static boolean hasBadMate(final GATKRead rec) {
        return (rec.isPaired() && !rec.mateIsUnmapped() && !rec.getAssignedContig().equals(rec.getMateContig()));
    }
}
