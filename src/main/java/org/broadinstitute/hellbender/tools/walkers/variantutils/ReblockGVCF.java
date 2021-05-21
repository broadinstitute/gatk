package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Condense homRef blocks in a single-sample GVCF
 *
 * <p>
 * ReblockGVCF compresses a GVCF by merging hom-ref blocks that were produced using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller according to new GQ band parameters.  A joint callset produced with GVCFs reprocessed by ReblockGVCF will have
 * lower precision for hom-ref genotype qualities at variant sites, but the input data footprint can be greatly reduced
 * if the default GQ band parameters are used.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A HaplotypeCaller-produced GVCF to reblock
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A smaller GVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk ReblockGVCF \
 *   -R reference.fasta \
 *   -V sample1.g.vcf \
 *   -O sample1.reblocked.g.vcf
 * </pre>
 *
 * Invocation as for use with GnarlyGenotyper in the "Biggest Practices"
 * <pre>
 *  gatk ReblockGVCF \
 *    -R reference.fasta \
 *    -V sample1.g.vcf \
 *    -drop-low-quals \
 *    -rgq-threshold 10 \
 *    -do-qual-approx \
 *    -O sample1.reblocked.g.vcf
 *  * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only single-sample GVCF files produced by HaplotypeCaller can be used as input for this tool.</p>
 * <h3>Special note on ploidy</h3>
 * <p>This tool assumes diploid genotypes.</p>
 *
 */
@BetaFeature
@CommandLineProgramProperties(summary = "Compress a single-sample GVCF from HaplotypeCaller by merging homRef blocks using new GQ band parameters",
        oneLineSummary = "Condenses homRef blocks in a single-sample GVCF",
        programGroup = OtherProgramGroup.class,
        omitFromCommandLine = true)
@DocumentedFeature
public final class ReblockGVCF extends VariantWalker {

    private static final OneShotLogger logger = new OneShotLogger(ReblockGVCF.class);

    private final static int PLOIDY_TWO = 2;  //assume diploid genotypes

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    private GATKPath outputFile;

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Output the band lower bound for each GQ block instead of the min GQ -- for better compression
     */
    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS, doc = "Output the band lower bound for each GQ block regardless of the data it represents", optional = true)
    private boolean floorBlocks = false;

    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.GQ_BAND_LONG_NAME, shortName=HaplotypeCallerArgumentCollection.GQ_BAND_SHORT_NAME,
            doc="Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)")
    public List<Integer> GVCFGQBands = new ArrayList<>();
    {
        GVCFGQBands.add(20); GVCFGQBands.add(100);
    }

    @Advanced
    @Argument(fullName="drop-low-quals", shortName="drop-low-quals", doc="Exclude variants and homRef blocks that are GQ0 from the reblocked GVCF to save space; drop low quality/uncalled alleles")
    protected boolean dropLowQuals = false;

    @Advanced
    @Argument(fullName="rgq-threshold-to-no-call", shortName="rgq-threshold", doc="Reference genotype quality (PL[0]) value below which variant sites will be converted to GQ0 homRef calls")
    protected double rgqThreshold = 0.0;

    @Advanced
    @Argument(fullName="do-qual-score-approximation", shortName="do-qual-approx", doc="Add necessary INFO field annotation to perform QUAL approximation downstream; required for GnarlyGenotyper")
    protected boolean doQualApprox = false;

    @Advanced
    @Argument(fullName="allow-missing-hom-ref-data", doc="Fill in homozygous reference genotypes with no PLs and no GQ with PL=[0,0,0].  Necessary for input from Regeneron's WeCall variant caller.")
    protected boolean allowMissingHomRefData = false;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;
    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = Arrays.asList(GVCFWriter.GVCF_BLOCK, GATKVCFConstants.HAPLOTYPE_SCORE_KEY,
            GATKVCFConstants.INBREEDING_COEFFICIENT_KEY, GATKVCFConstants.MLE_ALLELE_COUNT_KEY,
            GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, GATKVCFConstants.EXCESS_HET_KEY);

    private VariantContextWriter vcfWriter;

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, AS_StandardAnnotation.class);
    }

    @Override
    public void onTraversalStart() {
        VCFHeader inputHeader = getHeaderForVariants();
        if (inputHeader.getGenotypeSamples().size() > 1) {
            throw new UserException.BadInput("ReblockGVCF is a single sample tool, but the input GVCF has more than 1 sample.");
        }
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeaders);
        // Remove GCVFBlocks, legacy headers, and annotations that aren't informative for single samples
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCFWriter.GVCF_BLOCK) ||
                (vcfHeaderLine.getKey().equals("INFO")) && ((VCFInfoHeaderLine)vcfHeaderLine).getID().equals(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED) ||  //remove old (maybe wrong type) and add new with deprecated note
                (vcfHeaderLine.getKey().equals("INFO")) && infoFieldAnnotationKeyNamesToRemove.contains(((VCFInfoHeaderLine)vcfHeaderLine).getID()));

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        genotypingEngine = createGenotypingEngine(new IndexedSampleList(inputHeader.getGenotypeSamples()));
        createAnnotationEngine();

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED));  //NOTE: this is deprecated, but keep until we reprocess all GVCFs
        if (inputHeader.hasInfoLine(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED));
        }

        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        VariantContextWriter writer = createVCFWriter(outputFile);

        try {
            vcfWriter = new GVCFWriter(writer, new ArrayList<Number>(GVCFGQBands), PLOIDY_TWO, floorBlocks);
        } catch ( IllegalArgumentException e ) {
            throw new IllegalArgumentException("GQBands are malformed: " + e.getMessage(), e);
        }
        vcfWriter.writeHeader(new VCFHeader(headerLines, inputHeader.getGenotypeSamples()));

        if (genotypeArgs.samplePloidy != PLOIDY_TWO) {
            throw new UserException.BadInput("The -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this is tool assumes a diploid sample");
        }
    }

    private HaplotypeCallerGenotypingEngine createGenotypingEngine(SampleList samples) {
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        // create the genotyping engine
        hcArgs.standardArgs.outputMode = OutputMode.EMIT_ALL_ACTIVE_SITES;
        hcArgs.standardArgs.annotateAllSitesWithPLs = true;
        hcArgs.standardArgs.genotypeArgs = genotypeArgs.clone();
        hcArgs.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;   //this is important to force emission of all alleles at a multiallelic site
        return new HaplotypeCallerGenotypingEngine(hcArgs, samples, true, false);

    }

    @VisibleForTesting
    protected void createAnnotationEngine() {
        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, false);
    }

    // get VariantContexts from input gVCFs and regenotype
    @Override
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        final VariantContext newVC = regenotypeVC(variant);
        if (newVC != null) {
            vcfWriter.add(newVC);
        }
    }

    /**
     * Re-genotype (and re-annotate) a VariantContext
     * Note that the GVCF write takes care of the actual homRef block merging based on {@code GVCFGQBands}
     *
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
     private VariantContext regenotypeVC(final VariantContext originalVC) {
        VariantContext result = originalVC;

        //Pass back ref-conf homRef sites/blocks to be combined by the GVCFWriter
        if (isHomRefBlock(result)) {
            return filterHomRefBlock(result);
        }

        //don't need to calculate quals for sites with no data whatsoever or sites already genotyped homRef,
        // but if STAND_CALL_CONF > 0 we need to drop low quality alleles and regenotype
         //Note that spanning deletion star alleles will be considered low quality
        if (result.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 && !isHomRefCall(result) && dropLowQuals) {
            result = genotypingEngine.calculateGenotypes(originalVC);
        }

        if (result == null) {
            return null;
        }

        //variants with PL[0] less than threshold get turned to homRef with PL=[0,0,0], shouldn't get INFO attributes
        //make sure we can call het variants with GQ >= rgqThreshold in joint calling downstream
        if(shouldBeReblocked(result)) {
            return lowQualVariantToGQ0HomRef(result, originalVC);
        }
        //high quality variant
        else {
            return cleanUpHighQualityVariant(result, originalVC);
        }
    }

    private boolean isHomRefBlock(final VariantContext result) {
        return result.getAlternateAlleles().contains(Allele.NON_REF_ALLELE) && result.hasAttribute(VCFConstants.END_KEY);
    }

    /**
     * determine if VC is a homRef "call", i.e. an annotated variant with non-symbolic alt alleles and homRef genotypes
     * we treat these differently from het/homVar calls or homRef blocks
     * @param result VariantContext to process
     * @return true if VC is a 0/0 call and not a homRef block
     */
    private boolean isHomRefCall(final VariantContext result) {
        final Genotype genotype = result.getGenotype(0);
        return genotype.isHomRef() && result.getLog10PError() != VariantContext.NO_LOG10_PERROR;
    }

    private VariantContext filterHomRefBlock(final VariantContext result) {
        final Genotype genotype = result.getGenotype(0);
        if (dropLowQuals && (!genotype.hasGQ() || genotype.getGQ() < rgqThreshold || genotype.getGQ() == 0)) {
            return null;
        }
        else if (genotype.isCalled() && genotype.isHomRef()) {
            if (!genotype.hasPL()) {
                if (genotype.hasGQ()) {
                    logger.warn("PL is missing for hom ref genotype at at least one position for sample " + genotype.getSampleName() + ": " + result.getContig() + ":" + result.getStart() +
                            ".  Using GQ to determine quality.");
                    final int gq = genotype.getGQ();
                    final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                    final VariantContextBuilder vcBuilder = new VariantContextBuilder(result);
                    vcBuilder.genotypes(gBuilder.GQ(gq).make());
                    return vcBuilder.make();
                } else {
                    final String message = "Homozygous reference genotypes must contain GQ or PL. Both are missing for hom ref genotype at "
                            + result.getContig() + ":" + result.getStart() + " for sample " + genotype.getSampleName() + ".";
                    if (allowMissingHomRefData) {
                        logger.warn(message);
                        final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                        final VariantContextBuilder vcBuilder = new VariantContextBuilder(result);
                        vcBuilder.genotypes(gBuilder.GQ(0).PL(new int[]{0,0,0}).make());
                        return vcBuilder.make();
                    } else {
                        throw new UserException.BadInput(message);
                    }
                }
            }
            return result;
        }
        else if (!genotype.isCalled() && genotype.hasPL() && genotype.getPL()[0] == 0) {
            return result;
        }
        else {
            return null;
        }
    }

    @VisibleForTesting
    protected boolean shouldBeReblocked(final VariantContext result) {
        final Genotype genotype = result.getGenotype(0);
        return !genotype.isCalled() || (genotype.hasPL() && genotype.getPL()[0] < rgqThreshold) || genotype.isHomRef();
    }

    /**
     * "reblock" a variant by converting its genotype to homRef, changing PLs, adding reblock END tags and other attributes
     * @param result  a variant already determined to be low quality
     * @param originalVC the variant context with the original, full set of alleles
     * @return
     */
    @VisibleForTesting
    public VariantContext lowQualVariantToGQ0HomRef(final VariantContext result, final VariantContext originalVC) {
        if(dropLowQuals && (!isHomRefCall(result) || !result.getGenotype(0).isCalled())) {
            return null;
        }
        if (!result.getGenotype(0).isCalled()) {
            return new VariantContextBuilder(result).attributes(null).make();
        }

        final Map<String, Object> attrMap = new HashMap<>();
        final GenotypeBuilder gb = changeCallToGQ0HomRef(result, attrMap);

        //there are some cases where there are low quality variants with homRef calls AND alt alleles!
        //TODO: the best thing would be to take the most likely alt's likelihoods
        if (isHomRefCall(originalVC)) {
            final Genotype genotype = result.getGenotype(0);
            final int[] idxVector = originalVC.getGLIndicesOfAlternateAllele(Allele.NON_REF_ALLELE);   //this is always length 3
            if (genotype.hasPL()) {
                final int[] multiallelicPLs = genotype.getPL();
                final int[] newPLs = new int[3];
                newPLs[0] = multiallelicPLs[idxVector[0]];
                newPLs[1] = multiallelicPLs[idxVector[1]];
                newPLs[2] = multiallelicPLs[idxVector[2]];
                gb.PL(newPLs);
            }
            if (genotype.hasAD()) {
                int depth = (int) MathUtils.sum(genotype.getAD());
                gb.DP(depth);
                gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, depth);
            }
        }

        VariantContextBuilder builder = new VariantContextBuilder(result);
        final Genotype newG = gb.make();
        return builder.alleles(Arrays.asList(newG.getAlleles().get(0), Allele.NON_REF_ALLELE)).unfiltered().log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap).genotypes(newG).make(); //genotyping engine will add lowQual filter, so strip it off
    }

    /**
     * Note that this modifies {@code attrMap} as a side effect
     * @param result a VC to be converted to a GQ0 homRef call
     * @param attrMap the new VC attribute map, to update the END tag as necessary
     * @return a GenotypeBuilder to make a 0/0 call with PLs=[0,0,0]
     */
    @VisibleForTesting
    protected GenotypeBuilder changeCallToGQ0HomRef(final VariantContext result, final Map<String, Object> attrMap) {
        Genotype genotype = result.getGenotype(0);
        Allele newRef = result.getReference();
        GenotypeBuilder gb = new GenotypeBuilder(genotype);
        //NB: If we're dropping a deletion allele, then we need to trim the reference and add an END tag with the vc stop position
        if (result.getReference().length() > 1) {
            attrMap.put(VCFConstants.END_KEY, result.getEnd());
            newRef = Allele.create(newRef.getBases()[0], true);
            gb.alleles(Collections.nCopies(PLOIDY_TWO, newRef));
        }
        //if GT is not homRef, correct it
        if (!isHomRefCall(result)) {
            gb.PL(new int[3]);  //3 for diploid PLs, automatically initializes to zero
            gb.GQ(0).noAD().alleles(Collections.nCopies(PLOIDY_TWO, newRef)).noAttributes();
        }
        return gb;
    }

    @VisibleForTesting
    protected VariantContext cleanUpHighQualityVariant(final VariantContext result, final VariantContext originalVC) {
        Map<String, Object> attrMap = new HashMap<>();
        Map<String, Object> origMap = originalVC.getAttributes();
        //copy over info annotations
        for(final InfoFieldAnnotation annotation : annotationEngine.getInfoAnnotations()) {
            for (final String key : annotation.getKeyNames()) {
                if (infoFieldAnnotationKeyNamesToRemove.contains(key)) {
                    continue;
                }
                if (origMap.containsKey(key)) {
                    attrMap.put(key, origMap.get(key));
                }
            }
            if (annotation instanceof ReducibleAnnotation) {
                for (final String rawKey : ((ReducibleAnnotation)annotation).getRawKeyNames()) {
                    if (infoFieldAnnotationKeyNamesToRemove.contains(rawKey)) {
                        continue;
                    }
                    if (origMap.containsKey(rawKey)) {
                        attrMap.put(rawKey, origMap.get(rawKey));
                    }
                }
            }
        }
        final Genotype genotype = result.getGenotype(0);
        VariantContextBuilder builder = new VariantContextBuilder(result);  //QUAL from result is carried through
        builder.attributes(attrMap);

        boolean allelesNeedSubsetting = result.getNAlleles() != originalVC.getNAlleles();
        List<Allele> allelesToDrop = new ArrayList<>();
        if (dropLowQuals) {
            //drop alleles that aren't called iff we're dropping low quality variants (mostly because this can introduce GVCF gaps if deletion alleles are dropped)
            for (final Allele currAlt : result.getAlternateAlleles()) {
                boolean foundMatch = false;
                for (final Allele gtAllele : genotype.getAlleles()) {
                    if (gtAllele.equals(currAlt, false)) {
                        foundMatch = true;
                        break;
                    }
                    if (gtAllele.equals(Allele.NON_REF_ALLELE)) {
                        if (dropLowQuals) { //don't regenotype, just drop it -- this is a GQ 0 case if ever I saw one
                            return null;
                        } else {
                            GenotypeBuilder gb = changeCallToGQ0HomRef(result, attrMap);
                            return builder.alleles(Arrays.asList(result.getReference(), Allele.NON_REF_ALLELE)).unfiltered().log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap).genotypes(gb.make()).make();
                        }
                    }
                }
                if (!foundMatch && !currAlt.isSymbolic()) {
                    allelesNeedSubsetting = true;
                    allelesToDrop.add(currAlt);
                }
            }
        }

        //remove any AD reads for the non-ref
        int nonRefInd = result.getAlleleIndex(Allele.NON_REF_ALLELE);
        boolean genotypesWereModified = false;
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        GenotypesContext newGenotypes = result.getGenotypes();
        Genotype g = genotype;
        if(g.hasAD()) {
            int[] ad = g.getAD();
            if (ad.length >= nonRefInd && ad[nonRefInd] > 0) { //only initialize a builder if we have to
                genotypesWereModified = true;
                GenotypeBuilder gb = new GenotypeBuilder(g);
                ad[nonRefInd] = 0;
                gb.AD(ad).DP((int) MathUtils.sum(ad));
                genotypesArray.add(gb.make());
                newGenotypes = GenotypesContext.create(genotypesArray);
            }
        }
        else {
            genotypesArray.add(g);
        }

        //all VCs should get new RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, but preserve deprecated keys if present
        if (!originalVC.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY)) {
            //we're going to approximate depth for MQ calculation with the site-level DP (should be informative and uninformative reads), which is pretty safe because it will only differ if reads are missing MQ
            final Integer rawMqValue = originalVC.hasAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED) ?
                    (int)Math.round(originalVC.getAttributeAsDouble(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0.0)) :
                    (int)Math.round(originalVC.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, 60.0) *
                            originalVC.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, 60.0) *
                            originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            attrMap.put(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,
                    String.format("%d,%d", rawMqValue, originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0)));
            if (originalVC.hasAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
                attrMap.put(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, originalVC.getAttributeAsDouble(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0)); //NOTE: this annotation is deprecated, but keep it here so we don't have to reprocess gnomAD v3 GVCFs again
                attrMap.put(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED, originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0)); //NOTE: this annotation is deprecated, but keep it here so we don't have to reprocess gnomAD v3 GVCFs again
            }
        }


        int[] relevantIndices = new int[result.getNAlleles()];
        List<Allele> newAlleleSet = new ArrayList<>();
        for(final Allele a : result.getAlleles()) {
            newAlleleSet.add(a);
        }
        if(allelesNeedSubsetting) {
            newAlleleSet.removeAll(allelesToDrop);
            builder.alleles(newAlleleSet);
            final GenotypesContext gc;
            if(!genotypesWereModified) {
                gc = AlleleSubsettingUtils.subsetAlleles(result.getGenotypes(), PLOIDY_TWO, result.getAlleles(), newAlleleSet, null, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0), false);
                builder.genotypes(gc);
            }
            else {
                gc = AlleleSubsettingUtils.subsetAlleles(newGenotypes, PLOIDY_TWO, result.getAlleles(), newAlleleSet, null, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0), false);
                builder.genotypes(gc);
            }
            g = gc.get(0);
            relevantIndices = newAlleleSet.stream().mapToInt(a -> originalVC.getAlleles().indexOf(a)).toArray();
        }

        //copy over info annotations
        for(final InfoFieldAnnotation annotation : annotationEngine.getInfoAnnotations()) {
            for (final String key : annotation.getKeyNames()) {
                if (infoFieldAnnotationKeyNamesToRemove.contains(key)) {
                    continue;
                }
                if (origMap.containsKey(key)) {
                    attrMap.put(key, origMap.get(key));
                }
            }
            if (annotation instanceof ReducibleAnnotation) {
                for (final String rawKey : ((ReducibleAnnotation)annotation).getRawKeyNames()) {
                    if (infoFieldAnnotationKeyNamesToRemove.contains(rawKey)) {
                        continue;
                    }
                    if (origMap.containsKey(rawKey)) {
                        if (allelesNeedSubsetting && AnnotationUtils.isAlleleSpecific(annotation)) {
                            List<String> alleleSpecificValues = AnnotationUtils.getAlleleLengthListOfString(originalVC.getAttributeAsString(rawKey, null));
                            final List<?> subsetList = alleleSpecificValues.size() > 0 ? AlleleSubsettingUtils.remapRLengthList(alleleSpecificValues, relevantIndices, "")
                                    : Collections.nCopies(relevantIndices.length, "");
                            attrMap.put(rawKey, AnnotationUtils.encodeAnyASListWithRawDelim(subsetList));
                        } else {
                            attrMap.put(rawKey, origMap.get(rawKey));
                        }
                    }
                }
            }
        }
        //do QUAL calcs after we potentially drop alleles
        if (doQualApprox && g.hasPL()) {
            attrMap.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, g.getPL()[0]);
            int varDP = QualByDepth.getDepth(result.getGenotypes(), null);
            if (varDP == 0) {  //prevent QD=Infinity case
                varDP = originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 1); //if there's no VarDP and no DP, just prevent Infs/NaNs and QD will get capped later
            }
            attrMap.put(GATKVCFConstants.VARIANT_DEPTH_KEY, varDP);
            if (annotationEngine.getInfoAnnotations().stream()
                    .anyMatch(infoFieldAnnotation -> infoFieldAnnotation.getClass().getSimpleName().equals("AS_QualByDepth"))) {
                     final List<String> quals = new ArrayList<>();
                for (final Allele alt : newAlleleSet) {
                    if (alt.isReference()) {
                        //GDB expects an empty field for ref
                        continue;
                    }
                    if (alt.equals(Allele.NON_REF_ALLELE)) {
                        quals.add("0");
                        continue;
                    }
                    final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(result.getGenotypes(),
                            HomoSapiensConstants.DEFAULT_PLOIDY, result.getAlleles(), Arrays.asList(result.getReference(), alt), null,
                            GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, result.getAttributeAsInt(VCFConstants.DEPTH_KEY,0), false);
                    if (gc.get(0).hasPL()) {
                        quals.add(Integer.toString(gc.get(0).getPL()[0]));
                    } else {  //AlleleSubsettingUtils can no-call genotypes with super duper low GQs
                        quals.add("0");
                    }
                }
                attrMap.put(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM +String.join(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, quals));
                List<Integer> as_varDP = AS_QualByDepth.getAlleleDepths(AlleleSubsettingUtils.subsetAlleles(result.getGenotypes(),
                        HomoSapiensConstants.DEFAULT_PLOIDY, result.getAlleles(), newAlleleSet, null,
                        GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, result.getAttributeAsInt(VCFConstants.DEPTH_KEY,0), false));
                if (as_varDP != null) {
                    attrMap.put(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, as_varDP.stream().map( n -> Integer.toString(n)).collect(Collectors.joining(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM)));
                }
            }
        } else {  //manually copy annotations that might be from reblocking and aren't part of AnnotationEngine
            if (result.hasAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY)) {
                attrMap.put(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, result.getAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
            }
            if (result.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY)) {
                attrMap.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, result.getAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
            }

        }
        attrMap.put(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, g.getAlleles().stream().anyMatch(Allele::isReference) ? Arrays.asList(0,1,0) : Arrays.asList(0,0,1)); //ExcessHet currently uses rounded/integer genotype counts, so do the same here
        builder.attributes(attrMap);

        if (allelesNeedSubsetting) {
            //only trim if we're subsetting alleles, and we only subset if we're allowed to drop sites, as per the -drop-low-quals arg
            return GATKVariantContextUtils.reverseTrimAlleles(builder.attributes(attrMap).unfiltered().make());
        }
        return builder.attributes(attrMap).genotypes(newGenotypes).unfiltered().make();
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
