package org.broadinstitute.hellbender.tools.walkers.mutect;

import autovalue.shaded.org.apache.commons.lang.math.NumberUtils;
import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.MathUtils;
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
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static java.lang.Math.pow;

/**
 * Created by davidben on 9/15/16.
 */
public final class Mutect2Engine implements AssemblyRegionEvaluator {

    public static final String TUMOR_SB_POWER_FWD_KEY =             "TUMOR_SB_POWER_FWD";
    //TODO: move to GATKVCFConstants in gatk public
    public static final String TLOD_FWD_KEY =                       "TLOD_FWD";
    public static final String TLOD_REV_KEY =                       "TLOD_REV";
    public static final String TUMOR_SB_POWER_REV_KEY =             "TUMOR_SB_POWER_REV";
    public static final String TRIALLELIC_SITE_FILTER_NAME =                  "triallelic_site"; //M2
    public static final String STRAND_ARTIFACT_FILTER_NAME =                  "strand_artifact"; // M2
    public static final String CLUSTERED_READ_POSITION_FILTER_NAME =          "clustered_read_position"; // M2
    public static final String MEDIAN_LEFT_OFFSET_KEY =              "MEDIAN_LEFT_OFFSET";
    public static final String MEDIAN_RIGHT_OFFSET_KEY =             "MEDIAN_RIGHT_OFFSET";
    public static final String MAD_MEDIAN_LEFT_OFFSET_KEY =          "MAD_LEFT_OFFSET";
    public static final String MAD_MEDIAN_RIGHT_OFFSET_KEY =         "MAD_RIGHT_OFFSET";
    private static final Logger logger = LogManager.getLogger(Mutect2Engine.class);
    private final static List<VariantContext> NO_CALLS = Collections.emptyList();

    private M2ArgumentCollection MTAC;
    private SAMFileHeader header;

    private static final int MIN_READ_LENGTH = 30;

    private byte MIN_TAIL_QUALITY;

    private SampleList samplesList;
    private boolean printTCGAsampleHeader = false;

    private CachingIndexedFastaSequenceFile referenceReader;
    private ReadThreadingAssembler assemblyEngine;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private SomaticGenotypingEngine genotypingEngine;
    private Optional<HaplotypeBAMWriter> haplotypeBAMWriter;
    private VariantAnnotatorEngine annotationEngine;
    private final Mutect2FilteringEngine filteringEngine;

    private AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer();


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
        filteringEngine = new Mutect2FilteringEngine(MTAC);

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

        //If the samples specified are exactly one normal and one tumor, use the TCGA VCF sample header format
        printTCGAsampleHeader = samplesList.numberOfSamples() == 2 && MTAC.normalSampleName != null;

        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            MTAC.annotationsToUse.add("ClusteredReadPosition");
        }
        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(MTAC.annotationGroupsToUse,
                MTAC.annotationsToUse,
                MTAC.annotationsToExclude,
                MTAC.dbsnp.dbsnp,
                MTAC.comps);

        assemblyEngine = AssemblyBasedCallerUtils.createReadThreadingAssembler(MTAC);

        MIN_TAIL_QUALITY = (byte)(MTAC.minBaseQualityScore - 1);

        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs);

        genotypingEngine = new SomaticGenotypingEngine(samplesList, MTAC, MTAC.tumorSampleName, MTAC.normalSampleName, MTAC.DEBUG_READ_NAME);
        genotypingEngine.setAnnotationEngine(annotationEngine);
        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(MTAC, header);

        // why isn't this a constructor (instead of initialize)?  Since the method is package-friendly
        trimmer.initialize(MTAC.assemblyRegionTrimmerArgs, header.getSequenceDictionary(), MTAC.debug,
                MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, false);

        if( MTAC.CONTAMINATION_FRACTION_FILE != null ) {
            MTAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(MTAC.CONTAMINATION_FRACTION_FILE, MTAC.CONTAMINATION_FRACTION, samplesList.asSetOfSamples(), logger));
        }
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

        //TODO: see comment in getM2HeaderLines about some M2 filters undefined in GATKVCFHEaderLines
        headerInfo.addAll(getM2HeaderLines());
        headerInfo.addAll(getSampleHeaderLines());

        // if printTCGAsampleHeader, we already checked for exactly 1 tumor and 1 normal in printTCGAsampleHeader assignment in initialize()
        final List<String> outputSampleNames = printTCGAsampleHeader ? Arrays.asList("TUMOR", "NORMAL") : samplesList.asListOfSamples();

        //TODO: this hack removes nulls due to the aforementioned problem with M2 header lines
        final Set<VCFHeaderLine> nonNullHeaderInfo = headerInfo.stream()
                .filter(line -> line != null).collect(Collectors.toSet());
        final VCFHeader vcfHeader = new VCFHeader(nonNullHeaderInfo, outputSampleNames);
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        vcfWriter.writeHeader(vcfHeader);
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


        if ( ! MTAC.doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }
        return headerInfo;
    }

    private Set<VCFHeaderLine> getSampleHeaderLines(){
        final Set<VCFHeaderLine> sampleLines = new HashSet<>();
        if (printTCGAsampleHeader) {
            //NOTE: This will only list the first bam file for each tumor/normal sample if there is more than one
            final Map<String, String> normalSampleHeaderAttributes = ImmutableMap.of("ID", "NORMAL", "SampleName", MTAC.normalSampleName);
            final Map<String, String> tumorSampleHeaderAttributes = ImmutableMap.of("ID", "TUMOR", "SampleName", MTAC.tumorSampleName);
            sampleLines.add(new VCFSimpleHeaderLine("SAMPLE", normalSampleHeaderAttributes));
            sampleLines.add(new VCFSimpleHeaderLine("SAMPLE", tumorSampleHeaderAttributes));
        }
        return sampleLines;
    }

    public List<VariantContext> callRegion(final AssemblyRegion originalAssemblyRegion, final FeatureContext featureContext ) {
        if ( MTAC.justDetermineActiveRegions || !originalAssemblyRegion.isActive() || originalAssemblyRegion.size() == 0 ) {
            return NO_CALLS;
        }
        logReadInfo(MTAC.DEBUG_READ_NAME, originalAssemblyRegion.getReads(), "Present in original active region");

        final AssemblyRegion assemblyActiveRegion = AssemblyBasedCallerUtils.assemblyRegionWithWellMappedReads(originalAssemblyRegion, MTAC.MIN_MAPPING_QUALITY_SCORE, header);

        logReadInfo(MTAC.DEBUG_READ_NAME, assemblyActiveRegion.getReads(), "Present in assembly active region");

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyActiveRegion, Collections.emptyList(), MTAC, header, samplesList, logger, referenceReader, assemblyEngine);
        final SortedSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(originalAssemblyRegion,allVariationEvents);
        if (!trimmingResult.isVariationPresent()) {
            return NO_CALLS;
        }
        logReadInfo(MTAC.DEBUG_READ_NAME, trimmingResult.getCallableRegion().getReads(), "Present in trimming result");

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        // it is conceivable that the region is active yet has no events upon assembling only the well-mapped reads
        if( ! assemblyResult.isVariationPresent() ) {
            return NO_CALLS;
        }

        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        logReadInfo(MTAC.DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping");

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalAssemblyRegion?
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);

        // TODO: this quantity and all downstream uses of it seems like it can be obtained from
        // TODO: ReadLikelihoods<Allele>::sampelReads
        final Map<String, List<GATKRead>> perSampleFilteredReadList = splitReadsBySample(filteredReads);
        logReadInfo(MTAC.DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping after filtering reads");

        final Map<String,List<GATKRead>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        //TODO: this obtains the old behavior where the second of a read pair was the one recorded
        final Map<String, Integer> ARreads_origNormalMQ = regionForGenotyping.getReads().stream()
                .collect(Collectors.toMap(GATKRead::getName, GATKRead::getMappingQuality, (read, mate) -> mate));

        // modify MAPQ scores in normal to be high so that we don't do any base quality score capping
        regionForGenotyping.getReads().stream().filter(this::isReadFromNormal).forEach(rec -> rec.setMappingQuality(60));

        final ReadLikelihoods<Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods,
                ARreads_origNormalMQ,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getSpan(),
                featureContext,
                header);

        writeBamOutput(assemblyResult, readLikelihoods, calledHaplotypes);

        if( MTAC.debug) { logger.info("----------------------------------------------------------------------------------"); }
        return applyAnnotationsAndFilters(calledHaplotypes, featureContext);
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

    private Set<String> calculateFilters(final FeatureContext featureContext, final VariantContext vc, final Map<String, Object> eventDistanceAttributes) {
        final Set<String> filters = new HashSet<>();

        final Integer eventCount = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY);
        final Integer maxEventDistance = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY);

        filteringEngine.applyPanelOfNormalsFilter(vc, filters, featureContext);

        // TODO: make the change to have only a single normal sample (but multiple tumors is ok...)
        if (hasNormal()) {
            final Genotype normalGenotype = vc.getGenotype(MTAC.normalSampleName);

            // NOTE: how do we get the non-ref depth here?
            final int normalAltCounts = normalGenotype.getAD()[1];
            final double normalF = (Double) normalGenotype.getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);
            int normalAltQualityScoreSum = 0;

            final Object qss = normalGenotype.getExtendedAttribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY);
            if (qss != null) {
                normalAltQualityScoreSum = (Integer) ((Object[]) qss)[1];
            } else {
                logger.error("Null qss at " + vc.getStart());
            }

            if ( (normalAltCounts > MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || normalF > MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION ) && normalAltQualityScoreSum > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
                filters.add(GATKVCFConstants.ALT_ALLELE_IN_NORMAL_FILTER_NAME);
            } else if ( eventCount > 1 && normalAltCounts >= 1) {
                filters.add(GATKVCFConstants.MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME);
            }
        }


        else if (eventCount >= 3) {
            filters.add(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME);
        }

        filteringEngine.applySTRFilter(vc, filters);

        // NOTE: what if there is a 3bp indel followed by a snp... we are comparing starts
        // so it would be thrown but it's really an adjacent event
        if ( eventCount >= 2 && maxEventDistance >= 3) {
            filters.add(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME);
        }

        filteringEngine.applyClusteredReadPositionFilter(vc, filters);

        // TODO: Move strand bias filter here

        return filters;
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

    protected Set<GATKRead> filterNonPassingReads( final AssemblyRegion assemblyRegion) {
        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : assemblyRegion.getReads() ) {
            //TODO: Takuto points out that this is questionable.  Let's think hard abut it.
            // KCIBUL: only perform read quality filtering on tumor reads...
            if (isReadFromNormal(rec)) {
                if( rec.getLength() < MIN_READ_LENGTH ) {
                    readsToRemove.add(rec);
                }
            } else if( rec.getLength() < MIN_READ_LENGTH || rec.getMappingQuality() < MTAC.MIN_MAPPING_QUALITY_SCORE || hasBadMate(rec) ||
                    (MTAC.keepRG != null && !rec.getReadGroup().equals(MTAC.keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        assemblyRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    protected Map<String, List<GATKRead>> splitReadsBySample( final Collection<GATKRead> reads ) {
        return AssemblyBasedCallerUtils.splitReadsBySample(samplesList, header, reads);
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

    private boolean isReadFromNormal(final GATKRead rec) {
        return MTAC.normalSampleName != null && ReadUtils.getSampleName(rec, header).equals(MTAC.normalSampleName);
    }

    final List<VariantContext> applyAnnotationsAndFilters(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes, final FeatureContext featureContext) {
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
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY, NumberUtils.min(eventDistances));
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY, maxEventDistance);
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
                final Genotype tumorGenotype = new GenotypeBuilder(originalVC.getGenotype(MTAC.tumorSampleName)).name("TUMOR").make();
                final Genotype normalGenotype = new GenotypeBuilder(originalVC.getGenotype(MTAC.normalSampleName)).name("NORMAL").make();
                vcb.genotypes(Arrays.asList(tumorGenotype, normalGenotype));
            }
            annotatedCalls.add(vcb.make());
        }
        return annotatedCalls;
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


    // If mates in a pair are mapping to different contigs, it is likely that at least one of them is in the wrong place.
    // Of course, it is also possible there is a structural variant, in which case we lose some sensitivity
    public static boolean hasBadMate(final GATKRead rec) {
        return (rec.isPaired() && !rec.mateIsUnmapped() && !rec.getAssignedContig().equals(rec.getMateContig()));
    }

}
