package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
 *
 * <p>
 * This tool is designed to perform joint genotyping on a single input, which may contain one or many samples. In any
 * case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or
 * `-ERC BP_RESOLUTION`.
 *
 *
 * <h3>Input</h3>
 * <p>
 * The GATK4 GenotypeGVCFs tool can take only one input track.  Options are 1) a single single-sample GVCF 2) a single
 * multi-sample GVCF created by CombineGVCFs or 3) a GenomicsDB workspace created by GenomicsDBImport.
 * A sample-level GVCF is produced by HaplotypeCaller with the `-ERC GVCF` setting.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A final VCF in which all samples have been jointly genotyped.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Perform joint genotyping on a singular sample by providing a single-sample GVCF or on a cohort by providing a combined multi-sample GVCF</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" GenotypeGVCFs \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V input.g.vcf.gz \
 *   -O output.vcf.gz
 * </pre>
 *
 * <h4>Perform joint genotyping on GenomicsDB workspace created with GenomicsDBImport</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" GenotypeGVCFs \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V gendb://my_database \
 *   -O output.vcf.gz \
 *   --tmp-dir=/path/to/large/tmp
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 *   <li>Only GVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call GVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</li>
 *   <li>Cannot take multiple GVCF files in one command.</li>
 *   <li>The amount of temporary disk storage required by GenomicsDBImport may exceed what is available in the default location: `/tmp`. The command line argument `--tmp-dir` can be used to specify an alternate temperary storage location with sufficient space.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
@CommandLineProgramProperties(summary = "Perform joint genotyping on a single-sample GVCF from HaplotypeCaller or a multi-sample GVCF from CombineGVCFs or GenomicsDBImport",
        oneLineSummary = "Perform joint genotyping on one or more samples pre-called with HaplotypeCaller",
        programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
public final class GenotypeGVCFs extends VariantLocusWalker {

    public static final String PHASED_HOM_VAR_STRING = "1|1";
    public static final String ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME = "only-output-calls-starting-in-intervals";
    public static final String ALL_SITES_LONG_NAME = "include-non-variant-sites";
    public static final String ALL_SITES_SHORT_NAME = "all-sites";
    public static final String KEEP_COMBINED_LONG_NAME = "keep-combined-raw-annotations";
    public static final String KEEP_COMBINED_SHORT_NAME = "keep-combined";
    private static final String GVCF_BLOCK = "GVCFBlock";
    private VCFHeader outputHeader;


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private File outputFile;

    @Argument(fullName=ALL_SITES_LONG_NAME, shortName=ALL_SITES_SHORT_NAME, doc="Include loci found to be non-variant after genotyping", optional=true)
    private boolean includeNonVariants = false;

    /**
     * Import all data between specified intervals.   Improves performance using large lists of intervals, as in exome
     * sequencing, especially if GVCF data only exists for specified intervals.  Use with
     * --only-output-calls-starting-in-intervals if input GVCFs contain calls outside the specified intervals.
     */
    @Argument(fullName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            shortName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            doc = "Boolean flag to import all data in between intervals.")
    private boolean mergeInputIntervals = false;

    /**
     * "Genotype" somatic GVCFs, outputting genotypes according to confidently called alt alleles, which may lead to inconsistent ploidy
     * Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.
     */
    @Argument(fullName= CombineGVCFs.SOMATIC_INPUT_LONG_NAME, doc = "Finalize input GVCF according to somatic (i.e. Mutect2) TLODs (BETA feature)")
    protected boolean somaticInput = false;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     */
    @Argument(fullName=M2ArgumentCollection.EMISSION_LOD_LONG_NAME, shortName = M2ArgumentCollection.EMISSION_LOG_SHORT_NAME,
    doc = "LOD threshold to emit variant to VCF.")
    protected double tlodThreshold = 3.5;  //allow for some lower quality variants


    /**
     * Margin of error in allele fraction to consider a somatic variant homoplasmic, i.e. if there is less than a 0.1% reference allele fraction, those reads are likely errors
     */
    @Argument(fullName=CombineGVCFs.ALLELE_FRACTION_DELTA_LONG_NAME, doc = "Margin of error in allele fraction to consider a somatic variant homoplasmic")
    protected double afTolerance = 1e-3;  //based on Q30 as a "good" base quality score

    /**
     * If specified, keep the combined raw annotations (e.g. AS_SB_TABLE) after genotyping.  This is applicable to Allele-Specific annotations
     */
    @Argument(fullName=KEEP_COMBINED_LONG_NAME, shortName = KEEP_COMBINED_SHORT_NAME, doc = "If specified, keep the combined raw annotations")
    protected boolean keepCombined = false;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * This option can only be activated if intervals are specified.
     */
    @Advanced
    @Argument(fullName= ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME,
            doc="Restrict variant output to sites that start within provided intervals",
            optional=true)
    private boolean onlyOutputCallsStartingInIntervals = false;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set
     * when appropriate. Note that dbSNP is not used in any way for the genotyping calculations themselves.
     */
    @ArgumentCollection
    private final DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private GenotypingEngine<?> genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private ReferenceConfidenceVariantContextMerger merger;

    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();

    // INFO Header names that require alt alleles
    final LinkedHashSet<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();

    private VariantContextWriter vcfWriter;

    /** these are used when {@link #onlyOutputCallsStartingInIntervals) is true */
    private List<SimpleInterval> intervals;

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
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class);
    }

    @Override
    public void onTraversalStart() {
        if (somaticInput) {
            logger.warn("Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.");
        }

        if (!includeNonVariants) {
            changeTraversalModeToByVariant();
        }

        final VCFHeader inputVCFHeader = getHeaderForVariants();

        if(onlyOutputCallsStartingInIntervals) {
            if( !hasUserSuppliedIntervals()) {
                throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required if --" + ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME + " was specified.");
            }
        }

        intervals = hasUserSuppliedIntervals() ? intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()) :
                Collections.emptyList();

        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples()); //todo should this be getSampleNamesInOrder?

        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, keepCombined);

        // Request INFO field annotations inheriting from RankSumTest and RMSAnnotation added to remove list
        for ( final InfoFieldAnnotation annotation :  annotationEngine.getInfoAnnotations() ) {
            if ( annotation instanceof RankSumTest ||
                    annotation instanceof AS_RMSMappingQuality ||
                    annotation instanceof RMSMappingQuality) {
                final List<String> keyNames = annotation.getKeyNames();
                if ( !keyNames.isEmpty() ) {
                    infoFieldAnnotationKeyNamesToRemove.add(keyNames.get(0));
                }
            }
        }

        // We only want the engine to generate the AS_QUAL key if we are using AlleleSpecific annotations.
        genotypingEngine = new MinimalGenotypingEngine(createUAC(), samples, new GeneralPloidyFailOverAFCalculatorProvider(genotypeArgs), annotationEngine.isRequestedReducibleRawKey(GATKVCFConstants.AS_QUAL_KEY));

        merger = new ReferenceConfidenceVariantContextMerger(annotationEngine, getHeaderForVariants(), somaticInput);

        if ( includeNonVariants ) {
            // Save INFO header names that require alt alleles
            for ( final VCFHeaderLine headerLine : inputVCFHeader.getMetaDataInInputOrder() ) {
                if (headerLine instanceof VCFInfoHeaderLine ) {
                    if (((VCFInfoHeaderLine) headerLine).getCountType() == VCFHeaderLineCount.A) {
                        infoHeaderAltAllelesLineNames.add(((VCFInfoHeaderLine) headerLine).getID());
                    }
                }
            }
        }

        setupVCFWriter(inputVCFHeader, samples);
    }

    private static boolean annotationShouldBeSkippedForHomRefSites(VariantAnnotation annotation) {
        return annotation instanceof RankSumTest || annotation instanceof RMSMappingQuality || annotation instanceof AS_RMSMappingQuality;
    }

    private void setupVCFWriter(final VCFHeader inputVCFHeader, final SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if (keepCombined) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_QUAL_KEY));
        }
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        outputHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        vcfWriter.writeHeader(outputHeader);
    }

    @Override
    public void apply(final Locatable loc, List<VariantContext> variants, ReadsContext reads, ReferenceContext ref, FeatureContext features) {

        final List<VariantContext> variantsToProcess = getVariantSubsetToProcess(loc, variants);

        ref.setWindow(10, 10); //TODO this matches the gatk3 behavior but may be unnecessary
        final VariantContext mergedVC = merger.merge(variantsToProcess, loc, includeNonVariants ? ref.getBase() : null, !includeNonVariants, false);
        final VariantContext regenotypedVC = somaticInput ? regenotypeSomaticVC(mergedVC, ref, features, includeNonVariants) :
                regenotypeVC(mergedVC, ref, features, includeNonVariants);
        if (regenotypedVC != null) {
            final SimpleInterval variantStart = new SimpleInterval(regenotypedVC.getContig(), regenotypedVC.getStart(), regenotypedVC.getStart());
            if (!GATKVariantContextUtils.isSpanningDeletionOnly(regenotypedVC) &&
                    (!onlyOutputCallsStartingInIntervals || intervals.stream().anyMatch(interval -> interval.contains (variantStart)))) {
                vcfWriter.add(regenotypedVC);
            }
        }
    }

    // If includeNonVariants is set, we're using group-by-locus traversal. To match GATK3 GenotypeGVCFs,
    // see if there is a variant in the overlapping group that starts exactly at the locus start position, and if so
    // prioritize and process only that variant. Otherwise process all of the overlapping variants.
    private List<VariantContext> getVariantSubsetToProcess(final Locatable loc, List<VariantContext> preProcessedVariants) {
        if (includeNonVariants) {
            final List<VariantContext> matchingStart =
                    preProcessedVariants.stream().filter(vc -> vc.getStart() == loc.getStart()).collect(Collectors.toList());
            if (matchingStart.size() == 0) {
                return preProcessedVariants;
            }
            else if (matchingStart.size() == 1) {
                return matchingStart;
            }
            // since this tool only accepts a single input source, there should never be
            // more than one variant at a given starting locus
            throw new IllegalStateException(
                    String.format(
                            "Variant input contains more than one variant starting at location: %s",
                            new SimpleInterval(matchingStart.get(0))));
        } else {
            return preProcessedVariants;
        }
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    private VariantContext regenotypeVC(final VariantContext originalVC, final ReferenceContext ref, final FeatureContext features, boolean includeNonVariants) {
        Utils.nonNull(originalVC);

        final VariantContext result;

        if ( originalVC.isVariant()  && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
            // only re-genotype polymorphic sites
            final VariantContext regenotypedVC = calculateGenotypes(originalVC);
            if (regenotypedVC == null || (!isProperlyPolymorphic(regenotypedVC) && !includeNonVariants)) {
                return null;
            }
            if (isProperlyPolymorphic(regenotypedVC) || includeNonVariants) {
                // Note that reversetrimAlleles must be performed after the annotations are finalized because the reducible annotation data maps
                // were generated and keyed on the un reverseTrimmed alleles from the starting VariantContexts. Thus reversing the order will make
                // it difficult to recover the data mapping due to the keyed alleles no longer being present in the variant context.
                final VariantContext withGenotypingAnnotations = addGenotypingAnnotations(originalVC.getAttributes(), regenotypedVC);
                final VariantContext withAnnotations = annotationEngine.finalizeAnnotations(withGenotypingAnnotations, originalVC);
                final int[] relevantIndices = regenotypedVC.getAlleles().stream().mapToInt(a -> originalVC.getAlleles().indexOf(a)).toArray();
                final VariantContext trimmed = GATKVariantContextUtils.reverseTrimAlleles(withAnnotations);
                final GenotypesContext updatedGTs = subsetAlleleSpecificFormatFields(outputHeader, trimmed.getGenotypes(), relevantIndices);
                result = new VariantContextBuilder(trimmed).genotypes(updatedGTs).make();
            } else if (includeNonVariants) {
                result = originalVC;
            } else {
                return null;
            }
        } else {
            result = originalVC;
        }


        // if it turned monomorphic then we either need to ignore or fix such sites
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        if (result.isPolymorphicInSamples()) {
            // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
            final VariantContext reannotated = annotationEngine.annotateContext(result, features, ref, null, a -> true);
            return new VariantContextBuilder(reannotated).genotypes(cleanupGenotypeAnnotations(reannotated, false)).make();
        } else if (includeNonVariants) {
            // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
            VariantContext reannotated = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            reannotated = annotationEngine.annotateContext(reannotated, features, ref, null, GenotypeGVCFs::annotationShouldBeSkippedForHomRefSites);
            return removeNonRefAlleles(reannotated);
        } else {
            return null;
        }
    }

    private VariantContext calculateGenotypes(VariantContext vc){
        /*
         * Query the VariantContext for the appropriate model.  If type == MIXED, one would want to use model = BOTH.
         * However GenotypingEngine.getAlleleFrequencyPriors throws an exception if you give it anything but a SNP or INDEL model.
         */
        final GenotypeLikelihoodsCalculationModel model = vc.getType() == VariantContext.Type.INDEL
                ? GenotypeLikelihoodsCalculationModel.INDEL
                : GenotypeLikelihoodsCalculationModel.SNP;
        return genotypingEngine.calculateGenotypes(vc, model, null);
    }

    /**
     * Remove NON-REF alleles from the variant context
     *
     * @param vc   the variant context
     * @return variant context with the NON-REF alleles removed if multiallelic or replaced with NO-CALL alleles if biallelic
     */
    private VariantContext removeNonRefAlleles(final VariantContext vc) {

        // If NON_REF is the only alt allele, ignore this site
        final List<Allele> newAlleles = new ArrayList<>();
        // Only keep alleles that are not NON-REF
        for ( final Allele allele : vc.getAlleles() ) {
            if ( !allele.equals(Allele.NON_REF_ALLELE) ) {
                newAlleles.add(allele);
            }
        }

        // If no alt allele, then remove INFO fields that require alt alleles
        if ( newAlleles.size() == 1 ) {
            final VariantContextBuilder builder = new VariantContextBuilder(vc).alleles(newAlleles);
            for ( final String name : infoHeaderAltAllelesLineNames ) {
                builder.rmAttributes(Arrays.asList(name));
            }
            return builder.make();
        } else {
            return vc;
        }
    }

    private GenotypesContext subsetAlleleSpecificFormatFields(final VCFHeader outputHeader, final GenotypesContext originalGs, final int[] relevantIndices) {
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());
        for (final Genotype g : originalGs) {
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            final Set<String> keys = g.getExtendedAttributes().keySet();
            for (final String key : keys) {
                final VCFFormatHeaderLine headerLine = outputHeader.getFormatHeaderLine(key);
                final Object attribute;
                if (headerLine.getCountType().equals(VCFHeaderLineCount.INTEGER) && headerLine.getCount() == 1) {
                    attribute = g.getAnyAttribute(key);
                }
                else {
                    attribute = ReferenceConfidenceVariantContextMerger.generateAnnotationValueVector(headerLine.getCountType(),
                            GATKProtectedVariantContextUtils.attributeToList(g.getAnyAttribute(key)), relevantIndices);
                }
                gb.attribute(key, attribute);
            }
            newGTs.add(gb.make());
        }
        return newGTs;
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    private VariantContext regenotypeSomaticVC(final VariantContext originalVC, final ReferenceContext ref, final FeatureContext features, boolean includeNonVariants) {
        Utils.nonNull(originalVC);

        final VariantContext result;
        if ( originalVC.isVariant()  && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
            result = callSomaticGenotypes(originalVC);
        } else if (includeNonVariants) {
            result = originalVC;
        } else {
            result = null;
        }
        return result;
    }


    /**
     * Drop low quality alleles and call genotypes
     * CombineGVCFs will convert calls to no-call (of varying ploidy, as is the case in somatic)
     *
     * @param vc input VariantContext with no-called genotypes
     * @return a VC with called genotypes and low quality alleles removed, may be null
     */
    private VariantContext callSomaticGenotypes(final VariantContext vc) {
        final List<Genotype> newGenotypes = new ArrayList<>();
        final GenotypesContext genotypes = vc.getGenotypes();
        final double[] perAlleleLikelihoodSums = new double[vc.getAlleles().size()];  //needs the ref for the subsetting utils

        for(final Genotype g : genotypes) {
            GenotypeBuilder gb = new GenotypeBuilder(g);
            final double[] tlodArray = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
            final double[] variantAFArray = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 0.0);
            double variantAFtotal = 0;
            final List<Allele> calledAlleles = new ArrayList<>();
            for(int i = 0; i < vc.getAlleles().size()-1; i++) {
                variantAFtotal += variantAFArray[i];
                if (tlodArray[i] > tlodThreshold) {
                    calledAlleles.add(vc.getAlternateAllele(i));
                    perAlleleLikelihoodSums[i+1] += tlodArray[i];
                }
            }
            //hack for weird Mutect2 ploidy -- if the variant is non-homoplasmic, call the reference allele too
            if(variantAFtotal < 1-afTolerance && (!g.hasAD() || g.getAD()[0] > 0)) {

                calledAlleles.add(0, vc.getReference());
            }
            //"ploidy" gets set according to the size of the alleles List in the Genotype
            gb.alleles(calledAlleles);
            newGenotypes.add(gb.make());
        }

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final VariantContext regenotypedVC = builder.genotypes(newGenotypes).make();

        final int maxAltAlleles = ((UnifiedArgumentCollection)genotypingEngine.getConfiguration()).genotypeArgs.MAX_ALTERNATE_ALLELES;
        List<Allele> allelesToKeep;

        //we need to make sure all alleles pass the tlodThreshold
        allelesToKeep = new ArrayList<>(perAlleleLikelihoodSums.length-1);
        allelesToKeep.add(vc.getReference());
        for (int i = 1; i < perAlleleLikelihoodSums.length; i++) {
            if (perAlleleLikelihoodSums[i] > tlodThreshold) {
                allelesToKeep.add(vc.getAlternateAllele(i-1));
            }
        }

        if (regenotypedVC.getAlternateAlleles().size() > maxAltAlleles) {
            allelesToKeep = AlleleSubsettingUtils.filterToMaxNumberOfAltAllelesBasedOnScores(maxAltAlleles, allelesToKeep, perAlleleLikelihoodSums);
        }

        if (allelesToKeep.size() == 1) {
            return null;
        }

        //if we didn't drop alleles then we're done!
        if (allelesToKeep.size() == regenotypedVC.getAlleles().size()) {
            return regenotypedVC;
        }

        final int[] relevantIndices = allelesToKeep.stream().mapToInt(a -> regenotypedVC.getAlleles().indexOf(a)).toArray();

        //do another pass over genotypes to drop the alleles that aren't called
        final GenotypesContext reducedGenotypes = AlleleSubsettingUtils.subsetSomaticAlleles(outputHeader, regenotypedVC.getGenotypes(), allelesToKeep, relevantIndices);
        final VariantContext subsetVC = builder.alleles(allelesToKeep).genotypes(reducedGenotypes).make();
        final VariantContext trimmedVC = GATKVariantContextUtils.trimAlleles(subsetVC, true, true);
        if (isProperlyPolymorphic(trimmedVC)) {
            return trimmedVC;
        }
        else {
            return null;
        }
    }


    /**
     * Determines whether the provided VariantContext has real alternate alleles.
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    public static boolean isProperlyPolymorphic(final VariantContext vc) {
        //obvious cases
        if (vc == null || vc.getAlternateAlleles().isEmpty()) {
            return false;
        } else if (vc.isBiallelic()) {
            return !(GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0)) || vc.isSymbolic());
        } else if (GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0)) && vc.getAlternateAllele(1).equals(Allele.NON_REF_ALLELE)){
            return false;
        } else {
            return true;
        }
    }

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations and any other annotations generated by the genotyper.
        final Map<String, Object> attrs = new LinkedHashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY)) {
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
        }
        if (newVC.hasAttribute(GATKVCFConstants.AS_QUAL_KEY)) {
            attrs.put(GATKVCFConstants.AS_QUAL_KEY, newVC.getAttribute(GATKVCFConstants.AS_QUAL_KEY));
        }
        return new VariantContextBuilder(newVC).attributes(attrs).make();
    }


    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     * 5. move GQ to RGQ if the site is monomorphic
     *
     * @param vc            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    @VisibleForTesting
    static List<Genotype> cleanupGenotypeAnnotations(final VariantContext vc, final boolean createRefGTs) {
        final GenotypesContext oldGTs = vc.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) ) {
                depth = parseInt(oldGT.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY));
                builder.DP(depth);
                attrs.remove(GATKVCFConstants.MIN_DP_FORMAT_KEY);
            }

            attrs.remove(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, PHASED_HOM_VAR_STRING);
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && vc.isVariant() ) {
                final int[] AD = new int[vc.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                // move the GQ to RGQ
                if (oldGT.hasGQ()) {
                    builder.noGQ();
                    attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
                }

                //keep 0 depth samples and 0 GQ samples as no-call
                if (depth > 0 && oldGT.hasGQ() && oldGT.getGQ() > 0) {
                    final List<Allele> refAlleles = Collections.nCopies(oldGT.getPloidy(), vc.getReference());
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    private static int parseInt(Object attribute){
        if( attribute instanceof String) {
            return Integer.parseInt((String)attribute);
        } else if ( attribute instanceof Number){
            return ((Number) attribute).intValue();
        } else {
            throw new IllegalArgumentException("Expected a Number or a String but found something else.");
        }
    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection(genotypeArgs);

        //whether to emit non-variant sites is not contained in genotypeArgs and must be passed to uac separately
        //Note: GATK3 uses OutputMode.EMIT_ALL_CONFIDENT_SITES when includeNonVariants is requested
        //GATK4 uses EMIT_ALL_SITES to ensure LowQual sites are emitted.
        uac.outputMode = includeNonVariants ? OutputMode.EMIT_ALL_SITES : OutputMode.EMIT_VARIANTS_ONLY;
        return uac;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
