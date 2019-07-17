package org.broadinstitute.hellbender.tools.evoquer;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.reflections.Reflections;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A stand-alone version of the Gnarly Genotyper, extracted and modified from
 * Laura's (@ldgauthier) branch, ldg_getVQSRinput.
 *
 * (Thanks, Laura!)
 *
 * This is currently a static class and it probably shouldn't be.
 * I made it this way to move quickly with the intention of fixing it later.
 * If you don't like this, come shriek at me and maybe I will fix it.
 *
 * Created by jonn on 4/23/19.
 */
final class GnarlyGenotyperEngine {

    private static final Logger logger = LogManager.getLogger(GnarlyGenotyperEngine.class);

    //==================================================================================================================
    // Private Static Members:

    private static final double  INDEL_QUAL_THRESHOLD   = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.INDEL_HETEROZYGOSITY);
    private static final double  SNP_QUAL_THRESHOLD     = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.SNP_HETEROZYGOSITY);
    private static final int     ASSUMED_PLOIDY         = GATKVariantContextUtils.DEFAULT_PLOIDY;
    private static final int     PIPELINE_MAX_ALT_COUNT = 50;

    private static final boolean SUMMARIZE_PLs          = false;  //for very large numbers of samples, save on space and hail import time by summarizing PLs with genotype quality metrics
    private static final boolean stripASAnnotations     = false;

    private static final RMSMappingQuality mqCalculator = RMSMappingQuality.getInstance();

    // cache the ploidy 2 PL array sizes for increasing numbers of alts up to the maximum of PIPELINE_MAX_ALT_COUNT
    private static final int[]                                   likelihoodSizeCache = new int[ PIPELINE_MAX_ALT_COUNT + 1 ];
    private static final ArrayList<GenotypeLikelihoodCalculator> glcCache            = new ArrayList<>();

    static {
        final GenotypeLikelihoodCalculators GLCprovider = new GenotypeLikelihoodCalculators();

        //initialize PL size cache -- HTSJDK cache only goes up to 4 alts, but I need 6
        for (final int numAlleles : IntStream.rangeClosed(1, PIPELINE_MAX_ALT_COUNT + 1).boxed().collect(Collectors.toList())) {
            likelihoodSizeCache[numAlleles - 1] = GenotypeLikelihoods.numLikelihoods(numAlleles, ASSUMED_PLOIDY);
            glcCache.add(numAlleles - 1, GLCprovider.getInstance(ASSUMED_PLOIDY, numAlleles));
        }
    }

    //==================================================================================================================
    // Public Static Methods:

    /**
     * Finalizes the genotypes for the given intermediate {@link VariantContext} object, as created by
     * {@link EvoquerEngine#evokeInterval(SimpleInterval)}.
     *
     * This intermediate object is equivalent to one line of a GVCF for all samples.
     *
     * @param variant The intermediate {@link VariantContext} to process into a final {@link VariantContext} with a final genotype.
     * @return A {@link VariantContext} with a determined genotype across all samples or {@code null} if the given variant did not meet reporting thresholds.
     */
    // TODO: Can we remove this SuppressWarnings({"unchecked", "rawtypes"}) block?
    @SuppressWarnings({"unchecked", "rawtypes"})
    static VariantContext finalizeGenotype(final VariantContext variant) {
        //GenomicsDB merged all the annotations, but we still need to finalize MQ and QD annotations
        //builder gets the finalized annotations and dbBuilder gets the raw annotations for the database
        VariantContextBuilder builder = new VariantContextBuilder(mqCalculator.finalizeRawMQ(variant));

        SimpleInterval variantStart = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getStart());
        //return early if there's no non-symbolic ALT since GDB already did the merging
        if ( !variant.isVariant() || !GenotypeGVCFs.isProperlyPolymorphic(variant)
                || variant.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) == 0 ) {

            return null;
        }

        //return early if variant doesn't meet QUAL threshold
        if (!variant.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY) && !variant.hasAttribute(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY)) {
            logger.warn("Variant will not be output because it is missing the " + GATKVCFConstants.RAW_QUAL_APPROX_KEY + "/" + GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY + " key assigned by the ReblockGVCFs tool -- if the input did come from ReblockGVCFs, check the GenomicsDB vidmap.json annotation info");
            return null;
        }

        double QUALapprox = 0.0;
        if ( variant.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY) ) {
            QUALapprox = variant.getAttributeAsDouble(GATKVCFConstants.RAW_QUAL_APPROX_KEY, 0.0);
        } else {
            final List<Integer> as_quals = Arrays.stream(variant.getAttributeAsString(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, "").split(AnnotationUtils.AS_SPLIT_REGEX)).mapToInt(Integer::parseInt).boxed().collect(Collectors.toList());
            for ( final int as_qual : as_quals ) {
                QUALapprox += as_qual;  // FOR NOW!!! (Laura said this was ok as a stopgap)
            }
        }

        //TODO: do we want to apply the indel prior to mixed sites?
        final boolean isIndel = variant.getReference().length() > 1 || variant.getAlternateAlleles().stream().anyMatch(allele -> allele.length() > 1);
        final double sitePrior = isIndel ? HomoSapiensConstants.INDEL_HETEROZYGOSITY : HomoSapiensConstants.SNP_HETEROZYGOSITY;
        if((isIndel && QUALapprox < INDEL_QUAL_THRESHOLD) || (!isIndel && QUALapprox < SNP_QUAL_THRESHOLD)) {
            return null;
        }

        //GenomicsDB merged all the annotations, but we still need to finalize MQ and QD annotations
        //vcfBuilder gets the finalized annotations and annotationDBBuilder gets the raw annotations for the database
        final VariantContextBuilder vcfBuilder = new VariantContextBuilder(mqCalculator.finalizeRawMQ(variant));
        final VariantContextBuilder annotationDBBuilder = new VariantContextBuilder(variant);

        // final int variantDP = variant.getAttributeAsInt(GATKVCFConstants.VARIANT_DEPTH_KEY, 0);
        // double QD = QUALapprox / (double)variantDP;
        // vcfBuilder.attribute(GATKVCFConstants.QUAL_BY_DEPTH_KEY, QD).log10PError(QUALapprox/-10.0-Math.log10(sitePrior));
        vcfBuilder.rmAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY);
        //TODO: AS_QUALapprox needs non-ref removed, delimiters fixed

        //TODO: fix weird reflection logging?
        Reflections reflections = new Reflections("org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific");
        final Set<Class<? extends InfoFieldAnnotation>> allASAnnotations = reflections.getSubTypesOf(InfoFieldAnnotation.class);
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_StrandBiasTest.class));
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_RankSumTest.class));

        int[] SBsum = {0,0,0,0};

        final List<Allele> targetAlleles;
        final boolean removeNonRef;
        if (variant.getAlleles().contains(Allele.NON_REF_ALLELE)) { //Hail combine output doesn't give NON_REFs
            targetAlleles = variant.getAlleles().subList(0, variant.getAlleles().size() - 1);
            removeNonRef = true;
        }
        else {
            targetAlleles = variant.getAlleles();
            removeNonRef = false;
        }

        final Map<Allele, Integer> alleleCountMap = new HashMap<>();
        //initialize the count map
        for (final Allele a : targetAlleles) {
            alleleCountMap.put(a, 0);
        }

        //Get AC and SB annotations
        //remove the NON_REF allele and update genotypes if necessary
        final int[] rawGenotypeCounts = new int[3];
        final GenotypesContext calledGenotypes = iterateOnGenotypes(variant, targetAlleles, alleleCountMap, SBsum, removeNonRef, SUMMARIZE_PLs, variant.hasAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY) ? null : rawGenotypeCounts);
        Integer numCalledAlleles = 0;
        if (variant.hasGenotypes()) {
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
            if (variant.hasAttribute(GATKVCFConstants.SB_TABLE_KEY)) {
                SBsum = GATKProtectedVariantContextUtils.getAttributeAsIntArray(variant, GATKVCFConstants.SB_TABLE_KEY, () -> null, 0);
            }
            annotationDBBuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, variant.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
            annotationDBBuilder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, variant.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
            annotationDBBuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, variant.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
        }
        final List<Integer> gtCounts;
        if (variant.hasAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY)) {
            gtCounts = variant.getAttributeAsIntList(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, 0);
        } else {
            gtCounts = Arrays.stream(rawGenotypeCounts).boxed().collect(Collectors.toList());
        }
        final int refCount = Math.max(numCalledAlleles/2 - gtCounts.get(1) - gtCounts.get(2), 0);
        //homRefs don't get counted properly because ref blocks aren't annotated
        gtCounts.set(0, refCount);
        Pair<Integer, Double> eh = ExcessHet.calculateEH(variant, new GenotypeCounts(gtCounts.get(0), gtCounts.get(1), gtCounts.get(2)), numCalledAlleles/2);
        vcfBuilder.attribute(GATKVCFConstants.EXCESS_HET_KEY, String.format("%.4f", eh.getRight()));
        vcfBuilder.rmAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY);
        // vcfBuilder.attribute(GATKVCFConstants.FISHER_STRAND_KEY, FisherStrand.makeValueObjectForAnnotation(FisherStrand.pValueForContingencyTable(StrandBiasTest.decodeSBBS(SBsum))));
        vcfBuilder.attribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY, StrandOddsRatio.formattedValue(StrandOddsRatio.calculateSOR(StrandBiasTest.decodeSBBS(SBsum))));
        annotationDBBuilder.attribute(GATKVCFConstants.SB_TABLE_KEY, SBsum);

        vcfBuilder.genotypes(calledGenotypes);
        annotationDBBuilder.noGenotypes();
        vcfBuilder.alleles(targetAlleles);

        for (final Class c : allASAnnotations) {
            try {
                InfoFieldAnnotation annotation = (InfoFieldAnnotation) c.newInstance();
                if (annotation instanceof AS_StandardAnnotation && annotation instanceof ReducibleAnnotation) {
                    ReducibleAnnotation ann = (ReducibleAnnotation) annotation;
                    if (variant.hasAttribute(ann.getRawKeyName())) {
                        if (!stripASAnnotations) {
                            final Map<String, Object> finalValue = ann.finalizeRawData(vcfBuilder.make(), variant);
                            finalValue.forEach((key, value) -> vcfBuilder.attribute(key, value));
                            annotationDBBuilder.attribute(ann.getRawKeyName(), variant.getAttribute(ann.getRawKeyName()));
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
                    if (variant.hasAttribute(ann.getRawKeyName())) {
                        vcfBuilder.rmAttribute(ann.getRawKeyName());
                    }
                }
            }
            catch (Exception e) {
                throw new IllegalStateException("Something went wrong: ", e);
            }
        }
        if (variant.hasAttribute(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY) && variant.hasAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY)) {
            List<Integer> dps = Arrays.asList(variant.getAttributeAsString(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, "")
                    .split(AnnotationUtils.AS_SPLIT_REGEX)).stream().map(Integer::parseInt).collect(Collectors.toList());
            //TODO: do something with this value!
        }

        annotationDBBuilder.alleles(targetAlleles);

        return vcfBuilder.make();
   }


    //==================================================================================================================
    // Private Instance Methods:

    //assume input genotypes are diploid

    /**
     * Remove the NON_REF allele from the genotypes, updating PLs, ADs, and GT calls
     * @param vc the input variant with NON_REF
     * @return a GenotypesContext
     */
    private static GenotypesContext iterateOnGenotypes(final VariantContext vc, final List<Allele> targetAlleles,
                                                final Map<Allele,Integer> targetAlleleCounts, final int[] SBsum,
                                                final boolean nonRefReturned, final boolean summarizePLs,
                                                final int[] rawGenotypeCounts) {
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

            //re-tally genotype counts if they are missing from the original VC
            if (rawGenotypeCounts != null) {
                int altCount = (int)g.getAlleles().stream().filter(a -> !a.isReference()).count();
                rawGenotypeCounts[altCount]++;
            }
        }
        return mergedGenotypes;
    }


    /**
     * Some versions of GenomicsDB report no-calls as `0` or `.`
     * @param g
     * @return
     */
    private static boolean isGDBnoCall(Genotype g) {
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
}
