package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.reflections.Reflections;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Perform "quick and dirty" joint genotyping on one or more samples pre-called with HaplotypeCaller
 *
 * <p>
 * This tool is designed to perform joint genotyping on multiple samples pre-called with HaplotypeCaller to produce a
 * multi-sample callset in a super extra highly scalable manner. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller
 * with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
 *
 * <h3>Input</h3>
 * <p>
 * A GenomicsDB containing the samples to joint-genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A final VCF in which all samples have been jointly genotyped.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Perform joint genotyping on a set of GVCFs stored in a GenomicsDB</h4>
 * <pre>
 * gatk --javaOptions "-Xmx4g" GnarlyGenotyper \
 *   -R reference.fasta \
 *   -V gendb://genomicsdb \
 *   --onlyOutputCallsStartingInIntervals \
 *   -O output.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p><ul><li>Only GenomicsDB instances can be used as input for this tool.</li>
 * <li>To generate all the annotations necessary for VQSR, input variants must include the QUALapprox, VarDP and MQ_DP
 * </ul></p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool assumes all diploid genotypes.</p>
 *
 */
@CommandLineProgramProperties(summary = "Perform \"quick and dirty\" joint genotyping on one or more samples pre-called with HaplotypeCaller",
        oneLineSummary = "Perform \"quick and dirty\" joint genotyping on one or more samples pre-called with HaplotypeCaller",
        programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
public final class GnarlyGenotyper extends CombineGVCFs {

    private static final OneShotLogger warning = new OneShotLogger(GnarlyGenotyper.class);

    private static final boolean SUMMARIZE_PLs = false;  //for very large numbers of samples, save on space and hail import time by summarizing PLs with genotype quality metrics

    public static final int PIPELINE_MAX_ALT_COUNT = 6;

    private double INDEL_QUAL_THRESHOLD;
    private double SNP_QUAL_THRESHOLD;

    private static final int ASSUMED_PLOIDY = GATKVariantContextUtils.DEFAULT_PLOIDY;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private ReferenceConfidenceVariantContextMerger merger;

    private final RMSMappingQuality mqCalculator = RMSMappingQuality.getInstance();

    // cache the ploidy 2 PL array sizes for increasing numbers of alts up to the maximum of PIPELINE_MAX_ALT_COUNT
    private final int[] likelihoodSizeCache = new int[PIPELINE_MAX_ALT_COUNT + 1];
    private final static ArrayList<GenotypeLikelihoodCalculator> glcCache = new ArrayList<>();

    private final Set<Class<? extends InfoFieldAnnotation>> allASAnnotations = new HashSet<>();


    @Argument(fullName = "output-database-name", shortName = "output-db",
            doc="File to which the sites-only annotation database derived from these input samples should be written", optional=true)
    private String outputDbName = null;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    @Argument(fullName = "keep-all-sites", shortName = "keep-all",
            doc="Retain low quality and non-variant sites, applying appropriate filters", optional=true)
    private boolean keepAllSites = false;

    /**
     * This option can only be activated if intervals are specified.
     */
    @Advanced
    @Argument(fullName = GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME,
            doc="Restrict variant output to sites that start within provided intervals",
            optional=true)
    private boolean onlyOutputCallsStartingInIntervals = false;

    @Argument(fullName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            shortName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            doc = "Boolean flag to read in all data in between intervals.  Improves performance reading from GenomicsDB " +
                    "using large lists of +intervals, as in exome sequencing, especially if GVCF data only exists for " +
                    "specified intervals.")
    private boolean mergeInputIntervals = false;

    @Hidden
    @Argument(fullName = "strip-allele-specific-annotations", shortName = "strip-as", doc = "Remove raw AS values and don't calculate finalized values")
    private boolean stripASAnnotations = false;


    private VariantContextWriter vcfWriter;
    private VariantContextWriter annotationDBwriter = null;

    /** these are used when {@link #onlyOutputCallsStartingInIntervals) is true */
    private List<SimpleInterval> intervals;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * Get the largest interval per contig that contains the intervals specified on the command line.
     * @param getIntervals intervals to be transformed
     * @param sequenceDictionary used to validate intervals
     * @return a list of one interval per contig spanning the input intervals after processing and validation
     */
    @Override
    protected List<SimpleInterval> transformTraversalIntervals(final List<SimpleInterval> getIntervals, final SAMSequenceDictionary sequenceDictionary) {
        if (mergeInputIntervals) {
            return IntervalUtils.getSpanningIntervals(getIntervals, sequenceDictionary);
        } else {
            return getIntervals;
        }
    }

    @Override
    public boolean useVariantAnnotations() {
        return true;
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, AS_StandardAnnotation.class);
    }

    @Override
    protected GenomicsDBOptions getGenomicsDBOptions() {
        if (genomicsDBOptions == null) {
            genomicsDBOptions = new GenomicsDBOptions(referenceArguments.getReferencePath(), true, PIPELINE_MAX_ALT_COUNT);
        }
        return genomicsDBOptions;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        final VCFHeader inputVCFHeader = getHeaderForVariants();

        if(onlyOutputCallsStartingInIntervals) {
            if( !intervalArgumentCollection.intervalsSpecified()) {
                throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required if --" + GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME + " was specified.");
            }
        }
        intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()) :
                Collections.emptyList();

        setupVCFWriter(inputVCFHeader, getSamplesForVariants());

        //we don't apply the prior to the QUAL approx in ReblockGVCF, so do it here
        INDEL_QUAL_THRESHOLD = genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING - 10*Math.log10(genotypeArgs.indelHeterozygosity);
        SNP_QUAL_THRESHOLD = genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING - 10*Math.log10(genotypeArgs.snpHeterozygosity);

        if (!SUMMARIZE_PLs) {
            GenotypeLikelihoodCalculators GLCprovider = new GenotypeLikelihoodCalculators();

            //initialize PL size cache -- HTSJDK cache only goes up to 4 alts, but I need 6
            for (final int numAlleles : IntStream.rangeClosed(1, PIPELINE_MAX_ALT_COUNT + 1).boxed().collect(Collectors.toList())) {
                likelihoodSizeCache[numAlleles - 1] = GenotypeLikelihoods.numLikelihoods(numAlleles, ASSUMED_PLOIDY);
                glcCache.add(numAlleles - 1, GLCprovider.getInstance(ASSUMED_PLOIDY, numAlleles));
            }
        }

        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, false);

        merger = new ReferenceConfidenceVariantContextMerger(annotationEngine, getHeaderForVariants(), false);

        Reflections reflections = new Reflections("org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific");
        //allASAnnotations.addAll(reflections.getSubTypesOf(InfoFieldAnnotation.class));  //we don't want AS_InbreedingCoeff
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_StrandBiasTest.class));
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_RankSumTest.class));
        allASAnnotations.add(AS_RMSMappingQuality.class);
        allASAnnotations.add(AS_QualByDepth.class);
    }

    private void setupVCFWriter(VCFHeader inputVCFHeader, Set<String> samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCFWriter.GVCF_BLOCK));

        //add header for new filter
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.MONOMORPHIC_FILTER_NAME));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        // add headers for annotations added by this tool
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.FISHER_STRAND_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.RMS_MAPPING_QUALITY_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = createVCFWriter(outputFile);
        if (outputDbName != null) {
            annotationDBwriter = createVCFWriter(new File(outputDbName));
        }

        final VCFHeader dbHeader = new VCFHeader(headerLines);
        if (SUMMARIZE_PLs) {
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.GENOTYPE_QUALITY_BY_ALLELE_BALANCE));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.GENOTYPE_QUALITY_BY_ALT_CONFIDENCE));
        }
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(samples));
        vcfWriter.writeHeader(vcfHeader);
        if (outputDbName != null) {
            annotationDBwriter.writeHeader(dbHeader);
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext ref) {
        final List<VariantContext> mergedVCs = new ArrayList<>();

        // If we need to stop at an intermediate site since the last apply, do so (caused by gvcfBlocks, contexts ending, etc...)
        if (!variantContextsOverlappingCurrentMerge.isEmpty()) {
            Locatable last = prevPos != null && prevPos.getContig().equals(variantContextsOverlappingCurrentMerge.get(0).getContig()) ? prevPos : variantContextsOverlappingCurrentMerge.get(0);
            // If on a different contig, close out all the queued states on the current contig
            int end = last.getContig().equals(ref.getWindow().getContig())
                    ? ref.getInterval().getStart() - 1
                    : variantContextsOverlappingCurrentMerge.stream().mapToInt(VariantContext::getEnd).max().getAsInt();

            createIntermediateVariants(new SimpleInterval(last.getContig(), last.getStart(), end), mergedVCs);
        }

        mergeWithNewVCs(variantContexts, ref);

        // Update the stored reference if it has a later stop position than the current stored reference
        if ((storedReferenceContext == null) ||
                (!ref.getWindow().contigsMatch(storedReferenceContext.getWindow())) ||
                (storedReferenceContext.getWindow().getEnd() < ref.getWindow().getEnd())) {
            storedReferenceContext = ref;
        }

        for (final VariantContext vc : mergedVCs) {
            doTheThings(vc);
        }
    }

    private void doTheThings(final VariantContext mergedVC) {
        //GenomicsDB merged all the annotations, but we still need to finalize MQ and QD annotations
        //builder gets the finalized annotations and dbBuilder gets the raw annotations for the database
        VariantContextBuilder builder = new VariantContextBuilder(mqCalculator.finalizeRawMQ(mergedVC));

        SimpleInterval variantStart = new SimpleInterval(mergedVC.getContig(), mergedVC.getStart(), mergedVC.getStart());
        //return early if there's no non-symbolic ALT since GDB already did the merging
        if ( !mergedVC.isVariant() || !GenotypeGVCFs.isProperlyPolymorphic(mergedVC)
                || mergedVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) == 0
                || (onlyOutputCallsStartingInIntervals && !intervals.stream().anyMatch(interval -> interval.contains(variantStart)))) {
            if (keepAllSites) {
                builder.filter(GATKVCFConstants.MONOMORPHIC_FILTER_NAME);
                vcfWriter.add(builder.make());
            }
            return;
        }

        //return early if variant doesn't meet QUAL threshold
        if (!mergedVC.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY)) {
            warning.warn("Variant will not be output because it is missing the " + GATKVCFConstants.RAW_QUAL_APPROX_KEY + "key assigned by the ReblockGVCFs tool -- if the input did come from ReblockGVCFs, check the GenomicsDB vidmap.json annotation info");
        }
        final double QUALapprox = mergedVC.getAttributeAsDouble(GATKVCFConstants.RAW_QUAL_APPROX_KEY, 0.0);
        //TODO: do we want to apply the indel prior to mixed sites?
        final boolean isIndel = mergedVC.getReference().length() > 1 || mergedVC.getAlternateAlleles().stream().anyMatch(allele -> allele.length() > 1);
        final double sitePrior = isIndel ? HomoSapiensConstants.INDEL_HETEROZYGOSITY : HomoSapiensConstants.SNP_HETEROZYGOSITY;
        if((isIndel && QUALapprox < INDEL_QUAL_THRESHOLD) || (!isIndel && QUALapprox < SNP_QUAL_THRESHOLD)) {
            if (keepAllSites) {
                builder.filter(GATKVCFConstants.MONOMORPHIC_FILTER_NAME);
                vcfWriter.add(builder.make());
            }
            return;
        }

        //GenomicsDB merged all the annotations, but we still need to finalize MQ and QD annotations
        //vcfBuilder gets the finalized annotations and annotationDBBuilder gets the raw annotations for the database
        final VariantContextBuilder vcfBuilder = new VariantContextBuilder(mqCalculator.finalizeRawMQ(mergedVC));
        final VariantContextBuilder annotationDBBuilder = new VariantContextBuilder(mergedVC);

        final int variantDP = mergedVC.getAttributeAsInt(GATKVCFConstants.VARIANT_DEPTH_KEY, 0);
        double QD = QUALapprox / (double)variantDP;
        vcfBuilder.attribute(GATKVCFConstants.QUAL_BY_DEPTH_KEY, QD).log10PError(QUALapprox/-10.0-Math.log10(sitePrior));

        int[] SBsum = {0,0,0,0};

        final List<Allele> targetAlleles;
        final boolean removeNonRef;
        if (mergedVC.getAlleles().contains(Allele.NON_REF_ALLELE)) { //Hail combine output doesn't give NON_REFs
            targetAlleles = mergedVC.getAlleles().subList(0, mergedVC.getAlleles().size() - 1);
            removeNonRef = true;
        }
        else {
            targetAlleles = mergedVC.getAlleles();
            removeNonRef = false;
        }

        final Map<Allele, Integer> alleleCountMap = new HashMap<>();
        //initialize the count map
        for (final Allele a : targetAlleles) {
            alleleCountMap.put(a, 0);
        }

        //Get AC and SB annotations
        //remove the NON_REF allele and update genotypes if necessary
        final GenotypesContext calledGenotypes = iterateOnGenotypes(mergedVC, targetAlleles, alleleCountMap, SBsum, removeNonRef, SUMMARIZE_PLs);
        Integer numCalledAlleles = 0;
        if (mergedVC.hasGenotypes()) {
            for (final Allele a : targetAlleles) {
                numCalledAlleles += alleleCountMap.get(a);
            }
            final List<Integer> targetAlleleCounts = new ArrayList<>();
            final List<Double> targetAlleleFreqs = new ArrayList<>();
            for (final Allele a : targetAlleles) {
                if (!a.isReference()) {
                    targetAlleleCounts.add(alleleCountMap.get(a));
                    targetAlleleFreqs.add((double) alleleCountMap.get(a) / numCalledAlleles);
                }
            }
            vcfBuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, targetAlleleCounts.size() == 1 ? targetAlleleCounts.get(0) : targetAlleleCounts);
            vcfBuilder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, targetAlleleFreqs.size() == 1 ? targetAlleleFreqs.get(0) : targetAlleleFreqs);
            vcfBuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, numCalledAlleles);

            annotationDBBuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, targetAlleleCounts.size() == 1 ? targetAlleleCounts.get(0) : targetAlleleCounts);
            annotationDBBuilder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, targetAlleleFreqs.size() == 1 ? targetAlleleFreqs.get(0) : targetAlleleFreqs);
            annotationDBBuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, numCalledAlleles);
        } else {
            if (mergedVC.hasAttribute(GATKVCFConstants.SB_TABLE_KEY)) {
                SBsum = GATKProtectedVariantContextUtils.getAttributeAsIntArray(mergedVC, GATKVCFConstants.SB_TABLE_KEY, () -> null, 0);
            }
            annotationDBBuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, mergedVC.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
            annotationDBBuilder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, mergedVC.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
            annotationDBBuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, mergedVC.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
        }

        final List<Integer> gtCounts = Arrays.stream(((String)mergedVC.getAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY)).split(",")).mapToInt(Integer::parseInt).boxed().collect(Collectors.toList());
        final int refCount = numCalledAlleles/2 - gtCounts.get(1) - gtCounts.get(2);
        gtCounts.set(0, refCount);
        Pair<Integer, Double> eh = ExcessHet.calculateEH(mergedVC, new GenotypeCounts(gtCounts.get(0), gtCounts.get(1), gtCounts.get(2)), numCalledAlleles/2);
        vcfBuilder.attribute(GATKVCFConstants.EXCESS_HET_KEY, String.format("%.4f", eh.getRight()));
        vcfBuilder.attribute(GATKVCFConstants.FISHER_STRAND_KEY, FisherStrand.makeValueObjectForAnnotation(FisherStrand.pValueForContingencyTable(StrandBiasTest.decodeSBBS(SBsum))));
        vcfBuilder.attribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY, StrandOddsRatio.formattedValue(StrandOddsRatio.calculateSOR(StrandBiasTest.decodeSBBS(SBsum))));
        annotationDBBuilder.attribute(GATKVCFConstants.SB_TABLE_KEY, SBsum);
        //TODO: annotationDBBuilder add all the raw AS annotations

        vcfBuilder.genotypes(calledGenotypes);
        annotationDBBuilder.noGenotypes();
        vcfBuilder.alleles(targetAlleles);

        for (final Class c : allASAnnotations) {
            try {
                InfoFieldAnnotation annotation = (InfoFieldAnnotation) c.newInstance();
                if (annotation instanceof AS_StandardAnnotation && annotation instanceof ReducibleAnnotation) {
                    ReducibleAnnotation ann = (ReducibleAnnotation) annotation;
                    if (mergedVC.hasAttribute(ann.getRawKeyName())) {
                        if (!stripASAnnotations) {
                            final Map<String, Object> finalValue = ann.finalizeRawData(vcfBuilder.make(), mergedVC);
                            finalValue.forEach((key, value) -> vcfBuilder.attribute(key, value));
                            annotationDBBuilder.attribute(ann.getRawKeyName(), mergedVC.getAttribute(ann.getRawKeyName()));
                        }
                    }
                }
            }
            catch (Exception e) {
                throw new IllegalStateException("Something went wrong: ", e);
            }
        }
        //since AS_FS and AS_SOR share the same raw key, we have to wait to remove raw keys until all the finalized values are added
        for (final Class c : allASAnnotations) {
            try {
                InfoFieldAnnotation annotation = (InfoFieldAnnotation) c.newInstance();
                if (annotation instanceof AS_StandardAnnotation  && annotation instanceof ReducibleAnnotation) {
                    ReducibleAnnotation ann = (ReducibleAnnotation) annotation;
                    if (mergedVC.hasAttribute(ann.getRawKeyName())) {
                        vcfBuilder.rmAttribute(ann.getRawKeyName());
                    }
                }
            }
            catch (Exception e) {
                throw new IllegalStateException("Something went wrong: ", e);
            }
        }
        if (mergedVC.hasAttribute(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY) && mergedVC.hasAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY)) {
            List<Integer> dps = Arrays.asList(mergedVC.getAttributeAsString(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, "")
                    .split("|")).stream().map(Integer::parseInt).collect(Collectors.toList());
            
        }

        annotationDBBuilder.alleles(targetAlleles);



        VariantContext result = vcfBuilder.make();
        if (annotationDBwriter != null) {
            annotationDBwriter.add(annotationDBBuilder.make());  //we don't seem to have a sites-only option anymore, so do it manually
        }

        //variantStart = new SimpleInterval(result.getContig(), result.getStart(), result.getStart());
        if (!onlyOutputCallsStartingInIntervals || intervals.stream().anyMatch(interval -> interval.contains(variantStart))) {
            vcfWriter.add(result);
        }

    }

    //assume input genotypes are diploid

    /**
     * Remove the NON_REF allele from the genotypes, updating PLs, ADs, and GT calls
     * @param vc the input variant with NON_REF
     * @return a GenotypesContext
     */
    private GenotypesContext iterateOnGenotypes(final VariantContext vc, final List<Allele> targetAlleles,
                                                final Map<Allele,Integer> targetAlleleCounts, final int[] SBsum,
                                                final boolean nonRefReturned, final boolean summarizePLs) {
        final List<Allele> inputAllelesWithNonRef = vc.getAlleles();
        if(nonRefReturned && !inputAllelesWithNonRef.get(inputAllelesWithNonRef.size()-1).equals(Allele.NON_REF_ALLELE)) {
            throw new IllegalStateException("This tool assumes that the NON_REF allele is listed last, as in HaplotypeCaller GVCF output,"
            + " but that was not the case at position " + vc.getContig() + ":" + vc.getStart() + ".");
        }
        final GenotypesContext mergedGenotypes = GenotypesContext.create();

        int newPLsize = -1;
        if (!summarizePLs) {
            final int maximumAlleleCount = inputAllelesWithNonRef.size();
            final int numConcreteAlts = maximumAlleleCount - 1; //-1 for NON_REF
            if (maximumAlleleCount <= PIPELINE_MAX_ALT_COUNT) {
                newPLsize = likelihoodSizeCache[numConcreteAlts - 1]; //-1 for zero-indexed array
            } else {
                newPLsize = GenotypeLikelihoods.numLikelihoods(maximumAlleleCount, ASSUMED_PLOIDY);
            }
        }

        for ( final Genotype g : vc.getGenotypes() ) {
            final String name = g.getSampleName();
            if(g.getPloidy() != ASSUMED_PLOIDY && !isGDBnoCall(g)) {
                throw new UserException.BadInput("This tool assumes diploid genotypes, but sample " + name + " has ploidy "
                        + g.getPloidy() + " at position " + vc.getContig() + ":" + vc.getStart() + ".");
            }
            final Genotype calledGT;
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g);
            genotypeBuilder.name(name);
            if (isGDBnoCall(g) || g.getAllele(0).equals(Allele.NON_REF_ALLELE) || g.getAllele(1).equals(Allele.NON_REF_ALLELE)) {
                genotypeBuilder.alleles(GATKVariantContextUtils.noCallAlleles(ASSUMED_PLOIDY));
            }
            else if (nonRefReturned) {
                if (g.hasAD()) {
                    final int[] AD = trimADs(g, targetAlleles.size());
                    genotypeBuilder.AD(AD);
                }
                else if (g.countAllele(Allele.NON_REF_ALLELE) > 0) {
                    genotypeBuilder.alleles(GATKVariantContextUtils.noCallAlleles(ASSUMED_PLOIDY)).noGQ();
                }
            }
            if (g.hasPL()) {
                if (summarizePLs) {
                    summarizePLs(genotypeBuilder, g, vc);
                } else {
                    final int[] PLs = trimPLs(g, newPLsize);
                    genotypeBuilder.PL(PLs);
                    genotypeBuilder.GQ(MathUtils.secondSmallestMinusSmallest(PLs, 0));
                    //If GenomicsDB returns no-call genotypes like CombineGVCFs, then we need to actually find the GT from PLs
                    makeGenotypeCall(genotypeBuilder, GenotypeLikelihoods.fromPLs(PLs).getAsVector(), targetAlleles);
                }
            }
            final Map<String, Object> attrs = new HashMap<>(g.getExtendedAttributes());
            attrs.remove(GATKVCFConstants.MIN_DP_FORMAT_KEY);
            //attrs.remove(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
            calledGT = genotypeBuilder.attributes(attrs).make();
            mergedGenotypes.add(calledGT);

            if (g.hasAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                try {
                    @SuppressWarnings("unchecked")
                    final List<Integer> sbbsList = (ArrayList<Integer>) g.getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
                    MathUtils.addToArrayInPlace(SBsum, Ints.toArray(sbbsList));
                }
                catch (ClassCastException e) {
                    throw new IllegalStateException("The GnarlyGenotyper tool assumes that input variants come from " +
                            "GenomicsDB and have SB FORMAT fields that have already been parsed into ArrayLists.");
                }
            }

            //running total for AC values
            for (int i = 0; i < ASSUMED_PLOIDY; i++) {
                Allele a = calledGT.getAllele(i);
                int count = targetAlleleCounts.containsKey(a) ? targetAlleleCounts.get(a) : 0;
                if (!a.equals(Allele.NO_CALL)) {
                    targetAlleleCounts.put(a,count+1);
                }
            }
        }
        return mergedGenotypes;
    }

    /**
     *
     * @param gb
     * @param genotypeLikelihoods
     * @param allelesToUse
     */
    private static void makeGenotypeCall(final GenotypeBuilder gb,
                                        final double[] genotypeLikelihoods,
                                        final List<Allele> allelesToUse) {

        if ( genotypeLikelihoods == null || !GATKVariantContextUtils.isInformative(genotypeLikelihoods) ) {
            gb.alleles(GATKVariantContextUtils.noCallAlleles(ASSUMED_PLOIDY)).noGQ();
        } else {
            final int maxLikelihoodIndex = MathUtils.maxElementIndex(genotypeLikelihoods);
            final GenotypeLikelihoodCalculator glCalc = glcCache.get(allelesToUse.size());
            final GenotypeAlleleCounts alleleCounts = glCalc.genotypeAlleleCountsAt(maxLikelihoodIndex);

            gb.alleles(alleleCounts.asAlleleList(allelesToUse));
            final int numAltAlleles = allelesToUse.size() - 1;
            if ( numAltAlleles > 0 ) {
                gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(maxLikelihoodIndex, genotypeLikelihoods));
            }
        }
    }

    /**
     * Save space in the VCF output by omitting the PLs and summarizing their info in ABGQ and ALTGQ
     * ABGQ is the next best (Phred-scaled) genotype likelihood over genotypes with the called alleles and
     * different allele counts (e.g. het vs homRef or homVar)
     * ALTGQ is the next best (Phred-scaled) genotype likelihood if one of the called alts is dropped from the VC
     * (e.g. a 0/2 het might become a 0/1 het if the 2 allele is removed)
     * @param gb a builder to be modified with ABGQ and ALTGQ
     * @param g the original genotype (should not have NON_REF called)
     * @param vc the original VariantContext including NON_REF
     */
    static void summarizePLs(final GenotypeBuilder gb,
                                    final Genotype g,
                                    final VariantContext vc) {
        final List<Allele> calledAlleles = g.getAlleles();
        final List<Integer> calledAllelePLPositions = getPLindicesForAlleles(vc, calledAlleles);

        final int[] PLs = g.getPL();
        //ABGQ is for GTs where both alleleIndex1 and alleleIndex2 are in calledAllelePLPositions
        //ALTGQ is for GTs where not both alleleIndex1 and alleleIndex2 are in calledAllelePLPositions
        int ABGQ = Integer.MAX_VALUE;
        int ALTGQ = Integer.MAX_VALUE;

        if (g.isHet()) {
            for (int i : calledAllelePLPositions) {
                if (PLs[i] == 0) {
                    continue;
                }
                if (PLs[i] < ABGQ) {
                    ABGQ = PLs[i];
                }
            }
        }
        //ABGQ can be any position that has the homozygous allele
        else {
            for (int i = 0; i < PLs.length; i++) {
                boolean match1 = false;
                boolean match2 = false;
                if (PLs[i] == 0) {
                    continue;
                }
                //all this is matching alleles based on their index in vc.getAlleles()
                GenotypeLikelihoods.GenotypeLikelihoodsAllelePair PLalleleAltArrayIndexes = GenotypeLikelihoods.getAllelePair(i); //this call assumes ASSUMED_PLOIDY is 2 (diploid)
                if (calledAllelePLPositions.contains(PLalleleAltArrayIndexes.alleleIndex1)) {
                    match1 = true;
                }
                if (calledAllelePLPositions.contains(PLalleleAltArrayIndexes.alleleIndex2)) {
                    match2 = true;
                }
                if (match1 || match2) {
                    if (PLs[i] < ABGQ) {
                        ABGQ = PLs[i];
                    }
                }
            }
            if (g.isHomRef()) {
                ALTGQ = ABGQ;
            }
        }

        if (!g.isHomRef()) {
            final Set<Allele> comparisonAlleles = new HashSet<>(vc.getAlleles());
            List<Integer> comparisonAllelePLPositions;
            if (!g.getAllele(0).isReference()) {
                comparisonAlleles.remove(g.getAllele(0));
                comparisonAllelePLPositions = getPLindicesForAlleles(vc, new ArrayList<>(comparisonAlleles));
                for (final int i : comparisonAllelePLPositions) {
                    if (PLs[i] < ALTGQ) {
                        ALTGQ = PLs[i];
                    }
                }
                comparisonAlleles.add(g.getAllele(0));
            }
            if (!g.getAllele(1).isReference()) {
                comparisonAlleles.remove(g.getAllele(1));
                comparisonAllelePLPositions = getPLindicesForAlleles(vc, new ArrayList<>(comparisonAlleles));
                for (final int i : comparisonAllelePLPositions) {
                    if (PLs[i] < ALTGQ) {
                        ALTGQ = PLs[i];
                    }
                }
            }
        }

        gb.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, PLs[0]);
        gb.attribute(GATKVCFConstants.GENOTYPE_QUALITY_BY_ALLELE_BALANCE, ABGQ);
        gb.attribute(GATKVCFConstants.GENOTYPE_QUALITY_BY_ALT_CONFIDENCE, ALTGQ);
        gb.noPL();
    }

    /**
     * Some versions of GenomicsDB report no-calls as `0` or `.`
     * @param g
     * @return
     */
    private boolean isGDBnoCall(Genotype g) {
        return g.getPloidy() == 1 && (g.getAllele(0).isReference() || g.getAllele(0).isNoCall());
    }

    /**
     * @param g  genotype from a VC including the NON_REF for which to update the PLs
     * @param newPLsize  number of PL entries for alleles without NON_REF
     * @return updated PLs
     */
    private static int[] trimPLs(final Genotype g, final int newPLsize) {
        final int[] oldPLs = g.getPL();
        final int[] newPLs = new int[newPLsize];
        System.arraycopy(oldPLs, 0, newPLs, 0, newPLsize);
        return newPLs;
    }

    /**
     * Trim the AD array to fit the set of alleles without NON_REF
     * Reads supporting the NON_REF will be dropped
     * @param g genotype from a VC including the NON_REF for which to update the ADs
     * @param newAlleleNumber number of alleles not including NON_REF
     * @return updated ADs
     */
    private static int[] trimADs(final Genotype g, final int newAlleleNumber) {
        final int[] oldADs = g.getAD();
        final int[] newADs = new int[newAlleleNumber];
        System.arraycopy(oldADs, 0, newADs, 0, newAlleleNumber);
        return newADs;
    }

     /**
     *
     * @param vc
     * @param calledAlleles should be size 2
     * @return variable-length list of PL positions for genotypes including {@code calledAlleles}
     * e.g. {0,1,2} for a REF/ALT0 call, {0,3,5} for a REF/ALT2 call, {0} for a REF/REF call, {2} for a ALT0/ALT0 call
     */
    private static List<Integer> getPLindicesForAlleles(final VariantContext vc, final List<Allele> calledAlleles) {
        final List<Integer> calledAllelePLPositions = new ArrayList<>();
        for (final Allele a : calledAlleles) {
            final int[] x = vc.getGLIndicesOfAlternateAllele(a);
            calledAllelePLPositions.addAll(Arrays.stream(x).boxed().collect(Collectors.toList()));
        }
        return calledAllelePLPositions.stream().distinct().collect(Collectors.toList());

    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
        if (annotationDBwriter != null) {
            annotationDBwriter.close();
        }
    }
}