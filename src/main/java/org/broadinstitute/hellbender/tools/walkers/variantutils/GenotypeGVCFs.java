package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Perform joint genotyping on gVCF files produced by HaplotypeCaller
 *
 * <p>
 * GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant discovery
 * (see Best Practices documentation for more details) using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller, or result from combining such gVCF files using CombineGVCFs. This tool performs the multi-sample
 * joint aggregation step and merges the records together in a sophisticated manner: at each position of the input
 * gVCFs, this tool will combine all spanning records, produce correct genotype likelihoods, re-genotype the newly
 * merged record, and then re-annotate it.</p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more HaplotypeCaller gVCFs to genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined, genotyped VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T GenotypeGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *   --variant sample2.g.vcf \
 *   -o output.vcf
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
public final class GenotypeGVCFs extends VariantWalker {

    /**
     * The gVCF files to merge together
     */
    @Argument(fullName="variant", shortName = "V", doc="One or more input gVCF files")
    public List<RodBindingCollection<VariantContext>> variantCollections;
    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    @Argument(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName="includeNonVariantSites", shortName="allSites", doc="Include loci found to be non-variant after genotyping", optional=true)
    public boolean INCLUDE_NON_VARIANTS = false;

    /**
     * Uniquify all sample names (intended for use with multiple inputs for the same sample)
     */
    @Hidden
    @Advanced
    @Argument(fullName="uniquifySamples", shortName="uniquifySamples", doc="Assume duplicate samples are present and uniquify all names with '.variant' and file number index")
    public boolean uniquifySamples = false;

   @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Which annotations to recompute for the combined output VCF file.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to recompute.  The single value 'none' removes the default annotations", required=false)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(
            new String[]{InbreedingCoeff.class.getSimpleName(),
                    FisherStrand.class.getSimpleName(),
                    QualByDepth.class.getSimpleName(),
                    ChromosomeCounts.class.getSimpleName(),
                    StrandOddsRatio.class.getSimpleName()}));

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", optional=true)
    protected List<String> annotationGroupsToUse = new ArrayList<>();


    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    // the genotyping engine
    private GenotypingEngine genotypingEngine;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    @Override
    public void onTraversalStart() {
        boolean inputsAreTagged = false;

        // collect the actual rod bindings into a list for use later
        for ( final RodBindingCollection<VariantContext> variantCollection : variantCollections ) {
            variants.addAll(variantCollection.getRodBindings());
            if (uniquifySamples) {
                for (final RodBinding<VariantContext> rb : variantCollection.getRodBindings()) {
                    //are inputs passed in with -V:fileTag ?
                    if (!rb.getTags().isEmpty()) inputsAreTagged = true;
                }
            }
        }
        //RodBinding tags are used in sample uniquification
        if (inputsAreTagged)
            logger.warn("Output uniquified VCF may not be suitable for input to CombineSampleData because input VCF(s) contain tags.");

        final Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(toolkit, variants);

        final GATKVariantContextUtils.GenotypeMergeType mergeType;
        if(uniquifySamples) {
            mergeType = GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY;
        }
        else
            mergeType = GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE;

        final SampleList samples = new IndexedSampleList(SampleUtils.getSampleList(vcfRods, mergeType));
        // create the genotyping engine
        genotypingEngine = new UnifiedGenotypingEngine(createUAC(), samples, toolkit.getGenomeLocParser(), GeneralPloidyFailOverAFCalculatorProvider.createThreadSafeProvider(toolkit, genotypeArgs, logger),
                toolkit.getArguments().BAQMode);
        // create the annotation engine
        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationGroupsToUse), annotationsToUse, Collections.<String>emptyList(), this, toolkit);

        // take care of the VCF headers
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp != null && dbsnp.dbsnp.isBound() )
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);

        final Set<String> sampleNameSet = SampleListUtils.asSet(samples);
        final VCFHeader vcfHeader = new VCFHeader(headerLines, sampleNameSet);
        vcfWriter.writeHeader(vcfHeader);

        //now that we have all the VCF headers, initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)
        annotationEngine.invokeAnnotationInitializationMethods(headerLines);

        logger.info("Notice that the -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this is automatically determined by the input variant files");
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ){
        final Locatable loc = variant;
        final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.merge(featureContext.getPrioritizedValue(variants, loc), loc, INCLUDE_NON_VARIANTS ? referenceContext.getBase() : null, true, uniquifySamples, annotationEngine);
        if ( combinedVC == null ) {
            return;
        }
        final VariantContext vc = regenotypeVC(featureContext, referenceContext, combinedVC);
        if (vc != null) {
            vcfWriter.add(vc);
        }
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     *
     * @param featureContext        the ref tracker
     * @param ref            the ref context
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    protected VariantContext regenotypeVC(final FeatureContext featureContext, final ReferenceContext ref, final VariantContext originalVC) {
        if ( originalVC == null ) throw new IllegalArgumentException("originalVC cannot be null");

        VariantContext rawResult = originalVC;

        // only re-genotype polymorphic sites
        if ( rawResult.isVariant() ) {
            VariantContext regenotypedVC = genotypingEngine.calculateGenotypes(rawResult);
            if ( ! isProperlyPolymorphic(regenotypedVC) ) {
                if (!INCLUDE_NON_VARIANTS)
                    return null;
            }
            else {
                rawResult = addGenotypingAnnotations(rawResult.getAttributes(), regenotypedVC);
            }
        }

        //At this point we should already have DP and AD annotated
        VariantContext result = annotationEngine.finalizeAnnotations(rawResult, originalVC);
        //do trimming after allele-specific annotation reduction or the mapping is difficult
        result = GATKVariantContextUtils.reverseTrimAlleles(result);

        // if it turned monomorphic then we either need to ignore or fix such sites
        boolean createRefGTs = false;
        if ( result.isMonomorphicInSamples() ) {
            if ( !INCLUDE_NON_VARIANTS )
                return null;
            createRefGTs = true;
        }

        // Re-annotate and fix/remove some of the original annotations.
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        if ( createRefGTs ) {
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            result = annotationEngine.annotateContext(tracker, ref, null, result);
        } else {
            result = annotationEngine.annotateContext(tracker, ref, null, result);
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, false)).make();
        }

        return result;
    }

    /**
     * Determines whether the provided VariantContext has real alternate alleles
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    private boolean isProperlyPolymorphic(final VariantContext vc) {
        return ( vc != null &&
                !vc.getAlternateAlleles().isEmpty() &&
                (!vc.isBiallelic() ||
                        (!vc.getAlternateAllele(0).equals(Allele.SPAN_DEL) &&
                                !vc.getAlternateAllele(0).equals(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED))
                )
        );
    }

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations
        final Map<String, Object> attrs = new HashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY))
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));

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
     * @param VC            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    private List<Genotype> cleanupGenotypeAnnotations(final VariantContext VC, final boolean createRefGTs) {
        final GenotypesContext oldGTs = VC.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute("MIN_DP") ) {
                depth = Integer.parseInt((String)oldGT.getAnyAttribute("MIN_DP"));
                builder.DP(depth);
                attrs.remove("MIN_DP");
            }

            // move the GQ to RGQ
            if ( createRefGTs && oldGT.hasGQ() ) {
                builder.noGQ();
                attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
            }

            // remove SB
            attrs.remove("SB");

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "1|1");
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && VC.isVariant() ) {
                final int[] AD = new int[VC.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                final int ploidy = oldGT.getPloidy();
                final List<Allele> refAlleles = Collections.nCopies(ploidy,VC.getReference());

                //keep 0 depth samples and 0 GQ samples as no-call
                if (depth > 0 && oldGT.hasGQ() && oldGT.getGQ() > 0) {
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = genotypeArgs.clone();
        return uac;
    }

    @Override
    public Object onTraversalDone() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
        return null;
    }
}
