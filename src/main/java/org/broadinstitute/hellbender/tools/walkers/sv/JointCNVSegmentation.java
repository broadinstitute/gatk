package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.lang3.tuple.MutablePair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.samples.Sex;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@BetaFeature
@CommandLineProgramProperties(
        summary = "Gathers single-sample segmented gCNV VCFs, harmonizes breakpoints, and outputs a cohort VCF with genotypes.",
        oneLineSummary = "Combined single-sample segmented gCNV VCFs.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class JointCNVSegmentation extends MultiVariantWalkerGroupedOnStart {

    private SortedSet<String> samples;
    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary dictionary;
    private SVDepthOnlyCallDefragmenter defragmenter;
    private SVClusterEngine clusterEngine;
    private List<GenomeLoc> callIntervals;
    private String currentContig;
    private SampleDB sampleDB;
    private boolean doDefragmentation;

    protected final OneShotLogger oneShotLogger = new OneShotLogger(logger);

    public static final String MIN_QUALITY_LONG_NAME = "minimum-qs-score";
    public static final String MODEL_CALL_INTERVALS = "model-call-intervals";
    public static final String BREAKPOINT_SUMMARY_STRATEGY = "breakpoint-summary-strategy";

    @Argument(fullName = MIN_QUALITY_LONG_NAME, doc = "Minimum QS score to combine a variant segment")
    private int minQS = 20;

    @Argument(fullName = MODEL_CALL_INTERVALS, doc = "Intervals used for gCNV calls.  Should be preprocessed and filtered to line up with model calls. Required for exomes.")
    private File modelCallIntervalList;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY, doc = "Strategy to use for choosing a representative value for a breakpoint cluster.")
    private SVClusterEngine.BreakpointSummaryStrategy breakpointSummaryStrategy = SVClusterEngine.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined output file", optional=false)
    private File outputFile;

    @Argument(
            doc = "Reference copy-number on autosomal intervals.",
            fullName = PostprocessGermlineCNVCalls.AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int refAutosomalCopyNumber = 2;

    @Argument(
            doc = "Contigs to treat as allosomal (i.e. choose their reference copy-number allele according to " +
                    "the sample karyotype).",
            fullName = PostprocessGermlineCNVCalls.ALLOSOMAL_CONTIG_LONG_NAME,
            optional = true
    )
    private List<String> allosomalContigList = Arrays.asList("X","Y","chrX","chrY");

    /**
     * See https://software.broadinstitute.org/gatk/documentation/article.php?id=7696 for more details on the PED
     * format. Note that each -ped argument can be tagged with NO_FAMILY_ID, NO_PARENTS, NO_SEX, NO_PHENOTYPE to
     * tell the GATK PED parser that the corresponding fields are missing from the ped file.
     *
     */
    @Argument(fullName=StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName=StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for samples", optional=true)
    private File pedigreeFile = null;

    @Override
    public boolean doDictionaryCrossValidation() {
        return false;
    }

    //require a reference to do dictionary validation since there may be too many samples for cross-validating
    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        sampleDB = initializeSampleDB();

        if (sampleDB == null) {
            logger.warn("No pedigree file supplied for sex genotype calls. Ploidy will be inferred from segments VCF genotype ploidy.");
        }

        final Collection<String> samples = getSamplesForVariants();
        if (samples != null) {
            for (String sample : samples) {
                if (sampleDB.getSample(sample) == null) {
                    logger.warn("Sample " + sample + " is missing from the supplied pedigree file. Ploidy will be inferred from segments VCF genotype ploidy.");
                }
            }
        }


        dictionary = getBestAvailableSequenceDictionary();
        //dictionary will not be null because this tool requiresReference()

        final GenomeLocParser parser = new GenomeLocParser(this.dictionary);

        if (modelCallIntervalList == null) {
            callIntervals = null;
        } else {
        final List<GenomeLoc> inputCoverageIntervals = IntervalUtils.featureFileToIntervals(parser, modelCallIntervalList.getAbsolutePath());
        final List<GenomeLoc> inputTraversalIntervals = IntervalUtils.genomeLocsFromLocatables(parser,getTraversalIntervals());
            callIntervals = IntervalUtils.mergeListsBySetOperator(inputCoverageIntervals, inputTraversalIntervals, IntervalSetRule.INTERSECTION);
        }

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary, 0.8, callIntervals);
        clusterEngine = new SVClusterEngine(dictionary, true, breakpointSummaryStrategy);

        vcfWriter = getVCFWriter();
    }

    //TODO: this is the third copy of this method
    /**
     * Entry-point function to initialize the samples database from input data
     */
    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);  //strict will warn about missing samples
        if (pedigreeFile != null) {
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        }
        return sampleDBBuilder.getFinalSampleDB();
    }

    private VariantContextWriter getVCFWriter() {
        samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));

        VariantContextWriter writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /**
     * @param variantContexts  VariantContexts from driving variants with matching start position
     *                         NOTE: This will never be empty
     * @param referenceContext ReferenceContext object covering the reference of the longest spanning VariantContext
     * @param readsContexts
     */
    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        if (currentContig == null) {
            currentContig = variantContexts.get(0).getContig(); //variantContexts should have identical start, so choose 0th arbitrarily
        } else if (!variantContexts.get(0).getContig().equals(currentContig)) {
            processClusters();
            currentContig = variantContexts.get(0).getContig();
        }
        for (final VariantContext vc : variantContexts) {
            if (vc.getGenotypes().size() != 1) {
                oneShotLogger.warn("Multi-sample VCFs found, which are assumed to be pre-clustered. Skipping defragmentation.");
                doDefragmentation = false;
            } else {
                doDefragmentation = true;
            }
            final SVCallRecord record = SVCallRecordUtils.createDepthOnlyFromGCNVWithOriginalGenotypes(vc, minQS);
            if (record != null) {
                if (doDefragmentation) {
                    defragmenter.add(new SVCallRecordWithEvidence(record, Collections.emptyList(), Collections.emptyList(), Collections.emptyList(), null));
                } else {
                    clusterEngine.add(new SVCallRecordWithEvidence(record, Collections.emptyList(), Collections.emptyList(), Collections.emptyList(), null));
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        processClusters();
        return null;
    }

    private void processClusters() {
        if (!defragmenter.isEmpty()) {
            final List<SVCallRecord> defragmentedCalls = defragmenter.getOutput();
            defragmentedCalls.stream().forEachOrdered(clusterEngine::add);
        }
        //Jack and Isaac cluster first and then defragment
        final List<SVCallRecord> clusteredCalls = clusterEngine.getOutput();
        write(clusteredCalls);
    }

    private void write(final List<SVCallRecord> calls) {
        final ReferenceSequenceFile reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        final List<VariantContext> sortedCalls = calls.stream()
                .sorted(Comparator.comparing(c -> new SimpleInterval(c.getContigA(), c.getPositionA(), c.getPositionB()), //VCs have to be sorted by end as well
                        IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(record -> buildVariantContext(record, reference))
                .collect(Collectors.toList());
        Iterator<VariantContext> it = sortedCalls.iterator();
        ArrayList<VariantContext> overlappingVCs = new ArrayList<>();
        if (!it.hasNext()) {
            return;
        }
        VariantContext prev = it.next();
        overlappingVCs.add(prev);
        int clusterEnd = prev.getEnd();
        String clusterContig = prev.getContig();
        //gather groups of overlapping VCs and update the genotype copy numbers appropriately
        while (it.hasNext()) {
            final VariantContext curr = it.next();
            if (curr.getStart() < clusterEnd && curr.getContig().equals(clusterContig)) {
                overlappingVCs.add(curr);
                if (curr.getEnd() > clusterEnd) {
                    clusterEnd = curr.getEnd();
                }
            } else {
                final List<VariantContext> resolvedVCs = resolveVariantContexts(overlappingVCs);
                for (final VariantContext vc : resolvedVCs) { vcfWriter.add(vc); }
                overlappingVCs = new ArrayList<>();
                overlappingVCs.add(curr);
                clusterEnd = curr.getEnd();
                clusterContig = curr.getContig();
            }
        }
        //write out the last set of overlapping VCs
        final List<VariantContext> resolvedVCs = resolveVariantContexts(overlappingVCs);
        for (final VariantContext vc : resolvedVCs) { vcfWriter.add(vc); }
    }

    /**
     * Ensure genotype calls are consistent for overlapping variant contexts
     * Note that we assume that a sample will not occur twice with the same copy number because it should have been defragmented
     * @param overlappingVCs
     * @return
     */
    private List<VariantContext> resolveVariantContexts(final List<VariantContext> overlappingVCs) {
        Utils.nonNull(overlappingVCs);
        final List<VariantContext> resolvedVCs = new ArrayList<>();
        final Iterator<VariantContext> it = overlappingVCs.iterator();

        final Map<String, MutablePair<Integer, Integer>> sampleCopyNumbers = new LinkedHashMap<>();  //sampleName, copyNumber, endPos -- it's safe to just use position because if the VCs overlap then they must be on the same contig
        while (it.hasNext()) {
            final VariantContext curr = it.next();
            for (final Genotype g : curr.getGenotypes()) {
                //if this sample is in the table and we have a new variant for this sample, update the table
                final String s = g.getSampleName();
                if (sampleCopyNumbers.containsKey(s) && sampleCopyNumbers.get(s).getLeft() != Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString())) {
                    sampleCopyNumbers.get(s).setRight(Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()));
                    sampleCopyNumbers.get(s).setLeft(curr.getEnd());
                }
            }
            resolvedVCs.add(updateGenotypes(curr, sampleCopyNumbers));
            //update copy number table for subsequent VCs using variants genotypes from input VCs
            for (final Genotype g : curr.getGenotypes()) {
                if (g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.CN)) {
                    sampleCopyNumbers.put(g.getSampleName(), new MutablePair<>(Integer.parseInt(g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN).toString()), curr.getAttributeAsInt(VCFConstants.END_KEY, curr.getStart())));
                }
            }
        }
        return resolvedVCs;
    }

    /**
     *
     * @param vc VariantContext with just variant samples
     * @param sampleCopyNumbers may be modified to remove terminated variants
     * @return new VariantContext with AC and AF
     */
    private VariantContext updateGenotypes(final VariantContext vc, final Map<String, MutablePair<Integer, Integer>> sampleCopyNumbers) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final List<Genotype> newGenotypes = new ArrayList<>();
        final Allele vcRefAllele = vc.getReference();
        final Map<Allele, Long> alleleCountMap = new HashMap<>();
        alleleCountMap.put(GATKSVVCFConstants.DEL_ALLELE, 0L);
        alleleCountMap.put(GATKSVVCFConstants.DUP_ALLELE, 0L);
        int alleleNumber = 0;
        for (final String sample : samples) {
            final int samplePloidy;
            //"square off" the genotype matrix by adding homRef calls
            if (!sampleCopyNumbers.containsKey(sample) && !vc.hasGenotype(sample)) {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                samplePloidy = getSamplePloidy(sample, vc.getContig(), null);
                genotypeBuilder.alleles(GATKVariantContextUtils.makePloidyLengthAlleleList(samplePloidy, vcRefAllele));
                genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, samplePloidy);
                newGenotypes.add(genotypeBuilder.make());
            } else {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                final Genotype g = vc.getGenotype(sample); //may be null
                samplePloidy = getSamplePloidy(sample, vc.getContig(), g);
                final int copyNumber = g != null ? Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, samplePloidy).toString()) :
                        samplePloidy;
                final List<Allele> alleles;
                if (doDefragmentation || g == null) {  //if it's a multi-sample VCF we can trust the input genotypes
                    alleles = GATKSVVariantContextUtils.makeGenotypeAlleles(copyNumber, samplePloidy, vcRefAllele);
                } else {
                    alleles = g.getAlleles();
                    for (Allele a : alleles) {
                        if (a.isReference()) {
                            a = vcRefAllele;
                        }
                    }
                }
                genotypeBuilder.alleles(alleles);
                //check for genotype in VC because we don't want to count overlapping events (in sampleCopyNumbers map) towards AC
                if (vc.hasGenotype(sample) && alleles.contains(GATKSVVCFConstants.DEL_ALLELE)) {
                    final Long temp = alleleCountMap.get(GATKSVVCFConstants.DEL_ALLELE);
                    alleleCountMap.put(GATKSVVCFConstants.DEL_ALLELE, temp + alleles.stream().filter(Allele::isNonReference).count());
                } else if (vc.hasGenotype(sample) && (copyNumber > samplePloidy)) {
                    final Long temp = alleleCountMap.get(GATKSVVCFConstants.DUP_ALLELE);
                    alleleCountMap.put(GATKSVVCFConstants.DUP_ALLELE, temp + 1); //best we can do for dupes is carrier frequency
                }
                genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, copyNumber);
                if (sampleCopyNumbers.containsKey(sample)) {
                    if (sampleCopyNumbers.get(sample).getRight() > vc.getStart()) {
                        genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, sampleCopyNumbers.get(sample).getLeft());
                    }
                }
                newGenotypes.add(genotypeBuilder.make());
            }
            alleleNumber += samplePloidy;
        }
        builder.genotypes(newGenotypes);
        if (alleleNumber > 0) {
            if (vc.getAlternateAlleles().size() == 1) {
                final long AC;
                if (vc.getAlternateAllele(0).equals(GATKSVVCFConstants.DUP_ALLELE)) {
                    AC = alleleCountMap.get(GATKSVVCFConstants.DUP_ALLELE);
                } else {
                    AC = alleleCountMap.get(GATKSVVCFConstants.DEL_ALLELE);
                }
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, AC)
                        .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, Double.valueOf(AC) / alleleNumber)
                        .attribute(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
            } else {
                final List<Long> alleleCounts = new ArrayList<>();
                final List<Double> alleleFreqs = new ArrayList<>();
                //if we merged and del and a dupe from different callsets, then make sure the VC has both alleles
                //if (alleleCountMap.get(GATKSVVCFConstants.DEL_ALLELE) > 0 && alleleCountMap.get(GATKSVVCFConstants.DUP_ALLELE) > 0) {
                //    builder.alleles(vc.getReference().getDisplayString(), GATKSVVCFConstants.DUP_ALLELE.getDisplayString(), GATKSVVCFConstants.DEL_ALLELE.getDisplayString());
                //}
                for (final Allele a : builder.getAlleles()) {
                    if (a.isReference()) {
                        continue;
                    }
                    alleleCounts.add(alleleCountMap.containsKey(a) ? alleleCountMap.get(a) : 0L);
                    alleleFreqs.add(alleleCountMap.containsKey(a) ? Double.valueOf(alleleCountMap.get(a)) / alleleNumber : 0L);
                }
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, alleleCounts)
                    .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs)
                    .attribute(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
            }
        }
        return builder.make();
    }

    /**
     *
     * @param sampleName
     * @param contig
     * @param g may be null
     * @return
     */
    private int getSamplePloidy(final String sampleName, final String contig, final Genotype g) {
        if (!allosomalContigList.contains(contig)) {
            return refAutosomalCopyNumber;
        }
        int samplePloidy = 1;
        if (sampleDB == null || sampleDB.getSample(sampleName) == null) {
            if (g != null) {
                samplePloidy = g.getPloidy();
            } else {
                oneShotLogger.warn("Samples missing from pedigree assumed to have ploidy 1 on allosomes.");
            }

        }
        else if (sampleDB != null && sampleDB.getSample(sampleName) != null) {
            final Sex sampleSex = sampleDB.getSample(sampleName).getSex();
            if (contig.equals("X") || contig.equals("chrX")) {
                if (sampleSex.equals(Sex.FEMALE)) {
                    samplePloidy = 2;
                } else if (sampleSex.equals(Sex.MALE)) {
                    samplePloidy = 1;
                }
            } else if (contig.equals("Y") || contig.equals("chrY")) {
                if (sampleSex.equals(Sex.FEMALE)) {
                    samplePloidy = 0;
                } else if (sampleSex.equals(Sex.MALE)) {
                    samplePloidy = 1;
                }
            }
        }
        return samplePloidy;
    }

    public VariantContext buildVariantContext(final SVCallRecord call, final ReferenceSequenceFile reference) {
        Utils.nonNull(call);
        Utils.nonNull(reference);
        final List<Allele> outputAlleles = new ArrayList<>();
        final Allele refAllele = Allele.create(ReferenceUtils.getRefBaseAtPosition(reference, call.getContigA(), call.getPositionA()), true);
        outputAlleles.add(refAllele);
        if (!call.getType().equals(StructuralVariantType.CNV)) {
            outputAlleles.add(Allele.create("<" + call.getType().name() + ">", false));
        } else {
            outputAlleles.add(GATKSVVCFConstants.DEL_ALLELE);
            outputAlleles.add(GATKSVVCFConstants.DUP_ALLELE);
        }

        final VariantContextBuilder builder = new VariantContextBuilder("", call.getContigA(), call.getPositionA(), call.getPositionB(),
                outputAlleles);
        builder.attribute(VCFConstants.END_KEY, call.getPositionB());
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        if (call.getType().equals(StructuralVariantType.CNV)) {
            builder.attribute(VCFConstants.SVTYPE, "MCNV");
        } else {
            builder.attribute(VCFConstants.SVTYPE, call.getType());
        }
        final List<Genotype> genotypes = new ArrayList<>();
        for (final Genotype g : call.getGenotypes()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g);
            //update reference alleles
            List<Allele> newGenotypeAlleles = new ArrayList<>();
            for (Allele a : g.getAlleles()) {
                if (a.isReference()) {
                    newGenotypeAlleles.add(refAllele);
                } else {
                    newGenotypeAlleles.add(a);
                }
            }
            genotypeBuilder.alleles(newGenotypeAlleles);
             if (g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.CN)) {
                genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN));
            }
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder.make();
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
