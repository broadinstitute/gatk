package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Perform joint genotyping on gVCF files produced by HaplotypeCaller
 *
 * <p>
 * GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant discovery
 * (see Best Practices documentation for more details) using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller, or result from combining such gVCF files using CombineGVCFs. This tool will produce correct genotype
 * likelihoods, re-genotype the newly merged record, and then re-annotate it.</p>
 *
 * <h3>Input</h3>
 * <p>
 * One HaplotypeCaller gVCF to genotype
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined, genotyped VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk-launch GenotypeGVCFs \
 *   -R reference.fasta \
 *   -V sample1.g.vcf \
 *   -O output.vcf
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
@CommandLineProgramProperties(summary = "genotype a gvcf file to produce a vcf", oneLineSummary = "genotype a gvcf file", programGroup = VariantProgramGroup.class)
@DocumentedFeature
public final class GenotypeGVCFs extends VariantWalker {

    public static final String PHASED_HOM_VAR_STRING = "1|1";
    private static final String GVCF_BLOCK = "GVCFBlock";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private File outputFile;

    //@Argument(fullName="includeNonVariantSites", shortName="allSites", doc="Include loci found to be non-variant after genotyping", optional=true)
    //TODO This option is currently not supported.
    private boolean includeNonVariants = false;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Which annotations to recompute for the combined output VCF file.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to recompute.", optional=true)
    private List<String> annotationsToUse = new ArrayList<>();


    @Advanced
    @Argument(fullName="annotationsToExclude", shortName="AX", doc="One or more specific annotations to exclude from recomputation.", optional=true)
    private List<String> annotationsToExclude = new ArrayList<>();

    /**
     * Which groups of annotations to add to the output VCF file.
     * Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Advanced
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", optional=true)
    private List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{StandardAnnotation.class.getSimpleName()}));

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    private final DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private GenotypingEngine<?> genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private ReferenceConfidenceVariantContextMerger merger;

    private VariantContextWriter vcfWriter;

    //todo: remove this when the reducible annotation framework is in place
    //this is ok because RMS_MAPPING is stateless, but it needs to be instantiated in order to call its
    //overridden instance methods
    private static final RMSMappingQuality RMS_MAPPING_QUALITY = new RMSMappingQuality();

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples()); //todo should this be getSampleNamesInOrder?

        genotypingEngine = new MinimalGenotypingEngine(createUAC(), samples, new GeneralPloidyFailOverAFCalculatorProvider(genotypeArgs));

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbsnp.dbsnp, Collections.emptyList());

        merger = new ReferenceConfidenceVariantContextMerger();

        setupVCFWriter(inputVCFHeader, samples);
    }

    private static boolean annotationShouldBeSkippedForHomRefSites(VariantAnnotation annotation) {
        return annotation instanceof RankSumTest || annotation instanceof RMSMappingQuality || annotation instanceof AS_RMSMappingQuality;
    }

    private void setupVCFWriter(VCFHeader inputVCFHeader, SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        ref.setWindow(10, 10); //TODO this matches the gatk3 behavior but may be unnecessary
        final VariantContext mergedVC = merger.merge(Collections.singletonList(variant), variant, includeNonVariants ? ref.getBase() : null, true, false);
        final VariantContext regenotypedVC = regenotypeVC(mergedVC, ref, features, includeNonVariants);
        if (regenotypedVC != null) {
            vcfWriter.add(regenotypedVC);
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
            if (isProperlyPolymorphic(regenotypedVC)) {
                final VariantContext allelesTrimmed = GATKVariantContextUtils.reverseTrimAlleles(regenotypedVC);
                final VariantContext withAnnotations = addGenotypingAnnotations(originalVC.getAttributes(), allelesTrimmed);
                //TODO: remove this when proper support for reducible annotations is added
                result = RMS_MAPPING_QUALITY.finalizeRawMQ(withAnnotations);
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
            final VariantContext reannotated = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            return annotationEngine.annotateContext(reannotated, features, ref, null, GenotypeGVCFs::annotationShouldBeSkippedForHomRefSites);
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
     * Determines whether the provided VariantContext has real alternate alleles.
     *
     * There is a bit of a hack to handle the <NON-REF> case because it is not defined in htsjdk.Allele
     * We check for this as a biallelic symbolic allele.
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    @VisibleForTesting
    static boolean isProperlyPolymorphic(final VariantContext vc) {
        //obvious cases
        if (vc == null || vc.getAlternateAlleles().isEmpty()) {
            return false;
        } else if (vc.isBiallelic()) {
            return !(isSpanningDeletion(vc.getAlternateAllele(0)) || vc.isSymbolic());
        } else {
            return true;
        }
    }

    @VisibleForTesting
    static boolean isSpanningDeletion(final Allele allele){
        return allele.equals(Allele.SPAN_DEL) || allele.equals(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
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
        final Map<String, Object> attrs = new LinkedHashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY)) {
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
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
        return uac;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
