package org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.reflections.Reflections;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Guts of the GnarlyGenotyper
 *
 *
 */

public final class GnarlyGenotyperEngine {

    private final static double INDEL_QUAL_THRESHOLD = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.INDEL_HETEROZYGOSITY);
    private final static double SNP_QUAL_THRESHOLD = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.SNP_HETEROZYGOSITY);

    private static final int ASSUMED_PLOIDY = GATKVariantContextUtils.DEFAULT_PLOIDY;

    private static final RMSMappingQuality mqCalculator = RMSMappingQuality.getInstance();

    // cache the ploidy 2 PL array sizes for increasing numbers of alts up to the maximum of maxAltAllelesToOutput
    private int[] likelihoodSizeCache;
    private final ArrayList<GenotypeLikelihoodCalculator> glcCache = new ArrayList<>();
    private Set<Class<? extends InfoFieldAnnotation>> allASAnnotations;

    private final int maxAltAllelesToOutput;
    private final boolean summarizePls;  //for very large numbers of samples, save on space and hail import time by summarizing PLs with genotype quality metrics
    private final boolean keepAllSites;
    private final boolean stripASAnnotations;


    public GnarlyGenotyperEngine(final boolean keepAllSites, final int maxAltAllelesToOutput, final boolean summarizePls, final boolean stripASAnnotations) {
        this.maxAltAllelesToOutput = maxAltAllelesToOutput;
        this.summarizePls = summarizePls;
        this.keepAllSites = keepAllSites;
        this.stripASAnnotations = stripASAnnotations;

        if (!summarizePls) {
            final GenotypeLikelihoodCalculators GLCprovider = new GenotypeLikelihoodCalculators();

            //initialize PL size cache -- HTSJDK cache only goes up to 4 alts, but I need 6
            likelihoodSizeCache = new int[maxAltAllelesToOutput + 1 + 1]; //+1 for ref and +1 so index == numAlleles
            glcCache.add(null); //add a null at index zero because zero alleles (incl. ref) makes no sense
            for (final int numAlleles : IntStream.rangeClosed(1, maxAltAllelesToOutput + 1).boxed().collect(Collectors.toList())) {
                likelihoodSizeCache[numAlleles] = GenotypeLikelihoods.numLikelihoods(numAlleles, ASSUMED_PLOIDY);
                //GL calculator cache is indexed by the total number of alleles, including ref
                glcCache.add(numAlleles, GLCprovider.getInstance(ASSUMED_PLOIDY, numAlleles));
            }
        }

        //TODO: fix weird reflection logging?
        final Reflections reflections = new Reflections("org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific");
        allASAnnotations = reflections.getSubTypesOf(InfoFieldAnnotation.class);
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_StrandBiasTest.class));
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_RankSumTest.class));
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    public VariantContext finalizeGenotype(final VariantContext variant) {
        return finalizeGenotype(variant, null);
    }

    /**
     * Calculates final annotation value for VC and calls genotypes, if necessary
     * @param variant   expected to have <NON_REF>
     * @param annotationDBBuilder   may be null, otherwise modified to retain the combined form of the annotations
     * @return may be null if variant doesn't exceed the quality threshold
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    public VariantContext finalizeGenotype(final VariantContext variant, final VariantContextBuilder annotationDBBuilder) {
        //GenomicsDB or Evoquer merged all the annotations, but we still need to finalize MQ and QD annotations
        //return a VC with the finalized annotations and dbBuilder gets the raw annotations for the database

        // Prefer non-AS QUALapprox, but if we only have AS QUALapprox sum the values
        final double QUALapprox;
        if (variant.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY)) {
            QUALapprox = variant.getAttributeAsInt(GATKVCFConstants.RAW_QUAL_APPROX_KEY, 0);
        }
        else if (variant.hasAttribute(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY)) {
            List<Integer> alleleSpecificQualList = AS_QualByDepth.parseQualList(variant);
            QUALapprox = alleleSpecificQualList.stream()
                    .mapToInt(Integer::intValue)
                    .sum();
        }
        else {
            QUALapprox = 0;
        }

        //Don't apply the indel prior to mixed sites if there's a SNP, but don't count a '*' as a SNP
        final boolean hasSnpAllele = variant.getAlternateAlleles().stream().anyMatch(allele -> allele != Allele.SPAN_DEL && allele.length() == variant.getReference().length());
        final boolean isIndel = !hasSnpAllele;
        final double sitePrior = isIndel ? HomoSapiensConstants.INDEL_HETEROZYGOSITY : HomoSapiensConstants.SNP_HETEROZYGOSITY;
        System.out.println("KCIBUL -- in the qualapprox w/ " + QUALapprox + " for isIndel " + isIndel + " and SNP:" + SNP_QUAL_THRESHOLD + " and INDEL:" + INDEL_QUAL_THRESHOLD);
        if((isIndel && QUALapprox < INDEL_QUAL_THRESHOLD) || (!isIndel && QUALapprox < SNP_QUAL_THRESHOLD)) {
            if (keepAllSites) {
                final VariantContextBuilder builder = new VariantContextBuilder(mqCalculator.finalizeRawMQ(variant));
                //we don't need the rest of the annotations since this is filtered and won't go to VQSR
                builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
                builder.attribute(GATKVCFConstants.AC_ADJUSTED_KEY, 0);
                return builder.make();
            }
            return null;
        }

        //GenomicsDB merged all the annotations, but we still need to finalize MQ and QD annotations
        //vcfBuilder gets the finalized annotations and annotationDBBuilder (if present) gets the raw annotations for the database
        final VariantContext vcWithMQ = mqCalculator.finalizeRawMQ(variant);
        final VariantContextBuilder vcfBuilder = new VariantContextBuilder(vcWithMQ);

        //Because AS_StrandBias annotations both use and return the raw key
        final Map<String, Object> annotationsToBeModified = new HashMap<>(vcWithMQ.getAttributes());
        for (final Class c : allASAnnotations) {
            try {
                final InfoFieldAnnotation annotation = (InfoFieldAnnotation) c.getDeclaredConstructor().newInstance();
                if (annotation instanceof AS_StandardAnnotation && annotation instanceof ReducibleAnnotation) {
                    final ReducibleAnnotation ann = (ReducibleAnnotation) annotation;
                    if (variant.hasAttribute(ann.getPrimaryRawKey())) {
                        if (!stripASAnnotations) {
                            //here we still have the non-ref
                            final Map<String, Object> finalValue = ann.finalizeRawData(vcfBuilder.make(), variant);
                            finalValue.forEach((key, value) -> annotationsToBeModified.put(key, value));
                            if (annotationDBBuilder != null) {
                                annotationDBBuilder.attribute(ann.getPrimaryRawKey(), variant.getAttribute(ann.getPrimaryRawKey()));
                            }
                        }
                    }
                }
            }
            catch (final Exception e) {
                throw new IllegalStateException("Something went wrong: ", e);
            }
        }
        vcfBuilder.attributes(annotationsToBeModified);

        // tolerate lack of VarDP annotation
        if ( variant.hasAttribute(GATKVCFConstants.VARIANT_DEPTH_KEY) ) {
            final int variantDP = variant.getAttributeAsInt(GATKVCFConstants.VARIANT_DEPTH_KEY, 0);
            final double QD = QUALapprox / (double) variantDP;
            vcfBuilder.attribute(GATKVCFConstants.QUAL_BY_DEPTH_KEY, QD).log10PError(QUALapprox / -10.0 - Math.log10(sitePrior));
        }

        if (!keepAllSites) {
            System.out.println("KCIBUL -- in the QD Check");

            vcfBuilder.rmAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY);
        }

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
        final GenotypesContext calledGenotypes = iterateOnGenotypes(variant, targetAlleles, alleleCountMap, SBsum, removeNonRef, summarizePls, variant.hasAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY) ? null : rawGenotypeCounts);
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

            if (annotationDBBuilder != null) {
                annotationDBBuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, targetAlleleCounts.size() == 1 ? targetAlleleCounts.get(0) : targetAlleleCounts);
                annotationDBBuilder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, targetAlleleFreqs.size() == 1 ? targetAlleleFreqs.get(0) : targetAlleleFreqs);
                annotationDBBuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, numCalledAlleles);
            }
        } else {
            if (variant.hasAttribute(GATKVCFConstants.SB_TABLE_KEY)) {
                SBsum = VariantContextGetters.getAttributeAsIntArray(variant, GATKVCFConstants.SB_TABLE_KEY, () -> null, 0);
            }
            if (annotationDBBuilder != null) {
                annotationDBBuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, variant.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
                annotationDBBuilder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, variant.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
                annotationDBBuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, variant.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
            }
        }
        if (variant.hasAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY) || variant.hasGenotypes()) {
            final List<Integer> gtCounts;
            if (variant.hasAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY)) {
                gtCounts = variant.getAttributeAsIntList(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, 0);
            } else {
                gtCounts = Arrays.stream(rawGenotypeCounts).boxed().collect(Collectors.toList());
            }
            final int refCount = Math.max(numCalledAlleles / 2 - gtCounts.get(1) - gtCounts.get(2), 0);
            //homRefs don't get counted properly because ref blocks aren't annotated
            gtCounts.set(0, refCount);
            final Pair<Integer, Double> eh = ExcessHet.calculateEH(new GenotypeCounts(gtCounts.get(0), gtCounts.get(1), gtCounts.get(2)), numCalledAlleles / 2);
            vcfBuilder.attribute(GATKVCFConstants.EXCESS_HET_KEY, String.format("%.4f", eh.getRight()));
            vcfBuilder.rmAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY);
            if (annotationDBBuilder != null) {
                //replace db copy with updated hom ref count
                annotationDBBuilder.attribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, StringUtils.join(gtCounts, AnnotationUtils.LIST_DELIMITER));
            }
        }

        vcfBuilder.attribute(GATKVCFConstants.FISHER_STRAND_KEY, FisherStrand.makeValueObjectForAnnotation(FisherStrand.pValueForContingencyTable(StrandBiasTest.decodeSBBS(SBsum))));
        vcfBuilder.attribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY, StrandOddsRatio.formattedValue(StrandOddsRatio.calculateSOR(StrandBiasTest.decodeSBBS(SBsum))));
        vcfBuilder.genotypes(calledGenotypes);

        if (annotationDBBuilder != null) {
            annotationDBBuilder.attribute(GATKVCFConstants.SB_TABLE_KEY, SBsum);
            annotationDBBuilder.noGenotypes();
        }

        for (final Class c : allASAnnotations) {
            try {
                final InfoFieldAnnotation annotation = (InfoFieldAnnotation) c.getDeclaredConstructor().newInstance();
                if (annotation instanceof AS_StandardAnnotation && annotation instanceof ReducibleAnnotation) {
                    final ReducibleAnnotation ann = (ReducibleAnnotation) annotation;
                    //trim NON_REF out of AS values
                    //trim NON_REF out of AS values
                    if (variant.hasAttribute(ann.getRawKeyNames().get(0))) {
                        vcfBuilder.attribute(annotation.getKeyNames().get(0), trimASAnnotation(vcfBuilder.make(), targetAlleles, annotation.getKeyNames().get(0)));
                    }
                    if (!keepAllSites && variant.hasAttribute(ann.getRawKeyNames().get(0))) {
                        System.out.println("KCIBUL -- in the AS Check");

                        vcfBuilder.rmAttribute(ann.getRawKeyNames().get(0));
                    }
                }
            }
            catch (final Exception e) {
                throw new IllegalStateException("Something went wrong: ", e);
            }
        }
        //since AS_FS and AS_SOR share the same raw key, we have to wait to remove raw keys until all the finalized values are added
        if (!keepAllSites) {
            System.out.println("KCIBUL -- in the ??? Check");

            for (final Class c : allASAnnotations) {
                try {
                    final InfoFieldAnnotation annotation = (InfoFieldAnnotation) c.getDeclaredConstructor().newInstance();
                    if (annotation instanceof AS_StandardAnnotation && annotation instanceof ReducibleAnnotation) {
                        final ReducibleAnnotation ann = (ReducibleAnnotation) annotation;
                        for (final String rawKey : ann.getRawKeyNames()) {
                            if (variant.hasAttribute(rawKey)) {
                                vcfBuilder.rmAttribute(rawKey);
                            }
                        }
                    }
                } catch (final Exception e) {
                    throw new IllegalStateException("Something went wrong: ", e);
                }
            }
        }
        if (variant.hasAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY)) {
            vcfBuilder.attribute(GATKVCFConstants.AS_ALT_ALLELE_DEPTH_KEY, AS_QualByDepth.finalizeRawGVCFVarDPValues(
                    variant.getAttributeAsString(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, null), targetAlleles.size()));
            //Since this was added by ReblockGVCFs and not the AnnotationEngine, I feel justified in taking it out without using the AnnotationEngine
            vcfBuilder.rmAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY);
        }

        if (annotationDBBuilder != null) {
            annotationDBBuilder.alleles(targetAlleles);
        }

        vcfBuilder.alleles(targetAlleles);

        //instead of annotationDBBuilder.make(), we modify the builder passed in (if non-null)
        return vcfBuilder.make();
    }

    //assume input genotypes are diploid

    /**
     * Remove the NON_REF allele from the genotypes, updating PLs, ADs, and GT calls
     * @param vc the input variant with NON_REF
     * @return a GenotypesContext
     */
    @VisibleForTesting
    protected GenotypesContext iterateOnGenotypes(final VariantContext vc, final List<Allele> targetAlleles,
                                                final Map<Allele,Integer> targetAlleleCounts, final int[] SBsum,
                                                final boolean nonRefReturned, final boolean summarizePLs,
                                                final int[] rawGenotypeCounts) {
        final int maxAllelesToOutput = maxAltAllelesToOutput + 1;  //+1 for ref
        final List<Allele> inputAllelesWithNonRef = vc.getAlleles();
        if(nonRefReturned && !inputAllelesWithNonRef.get(inputAllelesWithNonRef.size()-1).equals(Allele.NON_REF_ALLELE)) {
            throw new IllegalStateException("This tool assumes that the NON_REF allele is listed last, as in HaplotypeCaller GVCF output,"
            + " but that was not the case at position " + vc.getContig() + ":" + vc.getStart() + ".");
        }
        final GenotypesContext mergedGenotypes = GenotypesContext.create();

        int newPLsize = -1;
        if (!summarizePLs) {
            final int maximumAlleleCount = inputAllelesWithNonRef.size();
            final int numConcreteAlts = maximumAlleleCount - 2; //-1 for NON_REF and -1 for ref
            if (maximumAlleleCount <= maxAllelesToOutput) {
                newPLsize = likelihoodSizeCache[numConcreteAlts + 1]; //cache is indexed by #alleles with ref; don't count NON_REF
            } else {
                newPLsize = GenotypeLikelihoods.numLikelihoods(numConcreteAlts + 1, ASSUMED_PLOIDY);
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
                    //If GenomicsDB returns no-call genotypes like CombineGVCFs (depending on the GenomicsDBExportConfiguration),
                    // then we need to actually find the GT from PLs
                    makeGenotypeCall(genotypeBuilder, GenotypeLikelihoods.fromPLs(PLs).getAsVector(), targetAlleles);
                }
            }
            final Map<String, Object> attrs = new HashMap<>(g.getExtendedAttributes());
            attrs.remove(GATKVCFConstants.MIN_DP_FORMAT_KEY);
            //attrs.remove(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
            calledGT = genotypeBuilder.attributes(attrs).make();
            mergedGenotypes.add(calledGT);

            if (g.hasAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                MathUtils.addToArrayInPlace(SBsum, getSBFieldAsIntArray(g));
            }

            //running total for AC values
            for (int i = 0; i < ASSUMED_PLOIDY; i++) {
                final Allele a = calledGT.getAllele(i);
                final int count = targetAlleleCounts.getOrDefault(a, 0);
                if (!a.equals(Allele.NO_CALL)) {
                    targetAlleleCounts.put(a,count+1);
                }
            }

            //re-tally genotype counts if they are missing from the original VC
            if (rawGenotypeCounts != null) {
                final int altCount = (int)g.getAlleles().stream().filter(a -> !a.isReference()).count();
                rawGenotypeCounts[altCount]++;
            }
        }
        return mergedGenotypes;
    }

    /**
     * For a genotype with likelihoods that has a no-call GT, determine the most likely genotype from PLs and set it
     * We use a GenotypeLikelihoodCalculator to convert from the best PL index to the indexes of the alleles for that
     * genotype so we can set the GenotypeBuilder with the alleles
     * @param gb    GenotypeBuilder to modify and pass back
     * @param genotypeLikelihoods   PLs to use to call genotype; count should agree with number of alleles in allelesToUse
     * @param allelesToUse  alleles in the parent VariantContext (with ref), because GenotypeBuilder needs the allele String rather than index
     */
    @VisibleForTesting
    protected void makeGenotypeCall(final GenotypeBuilder gb,
                                        final double[] genotypeLikelihoods,
                                        final List<Allele> allelesToUse) {
        final int maxAllelesToOutput = maxAltAllelesToOutput + 1; //+1 for ref

        if ( genotypeLikelihoods == null || !GATKVariantContextUtils.isInformative(genotypeLikelihoods) ) {
            gb.alleles(GATKVariantContextUtils.noCallAlleles(ASSUMED_PLOIDY)).noGQ();
        } else {
            final int maxLikelihoodIndex = MathUtils.maxElementIndex(genotypeLikelihoods);

            GenotypeLikelihoodCalculator glCalc;
            if ( allelesToUse.size() <= maxAllelesToOutput ) {
                glCalc = glcCache.get(allelesToUse.size());
            } else {
                final GenotypeLikelihoodCalculators GLCprovider = new GenotypeLikelihoodCalculators();
                glCalc = GLCprovider.getInstance(ASSUMED_PLOIDY, allelesToUse.size());
            }
            
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
            for (final int i : calledAllelePLPositions) {
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
                final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair PLalleleAltArrayIndexes = GenotypeLikelihoods.getAllelePair(i); //this call assumes ASSUMED_PLOIDY is 2 (diploid)
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
     * Some legacy versions of GenomicsDB report no-calls as `0` or `.`
     * @param g genotype to check
     * @return  true if this is a genotype that should be represented as a ploidy-aware, spec compliant no-call
     */
    private static boolean isGDBnoCall(final Genotype g) {
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
     *  Trim an annotation to the values representing the target alleles
     * @param variant   the VariantContext with annotations corresponding to the original alleles
     * @param targetAlleles the subset of alleles to retain
     * @param key   the key for the annotation of interest
     * @return  a String representing an array of allele-specific values matching targetAlleles
     */
    private static String trimASAnnotation(final VariantContext variant, final List<Allele> targetAlleles, final String key) {
        final int[] relevantIndices = targetAlleles.stream().filter(a -> !a.isReference()).mapToInt(a -> variant.getAlternateAlleles().indexOf(a)).toArray();
        if (!variant.hasAttribute(key)) {
            return null;
        }
        final List<String> annotationEntries = AnnotationUtils.decodeAnyASList(variant.getAttribute(key).toString());
        if (annotationEntries == null) {
            return null;
        }
        final List<String> returnString = new ArrayList<>();
        for (int i = 0; i < relevantIndices.length; i++) {
            if (relevantIndices[i] <= annotationEntries.size()-1) {
                returnString.add(annotationEntries.get(relevantIndices[i]));
            } else {
                returnString.add(VCFConstants.MISSING_VALUE_v4);
            }
        }
        return AnnotationUtils.encodeStringList(returnString);
    }


    /**
     *
     * @param vc input VariantContext
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
   * Parse SB field into an array of ints to pass to MathUtils. This method is required as we cannot always expect the
   * variant context to be fully decoded.
   * @param g genotype for a sample
   * @return int array for the given genotype
   */
  private static int[] getSBFieldAsIntArray(Genotype g) {
      Object sbbsObj = g.getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
      if (sbbsObj == null) {
          return new int[0];
      } else if (sbbsObj instanceof String) {
          final String[] sbbsStr = ((String)sbbsObj).split(",");
          try {
              return Arrays.stream(sbbsStr).map(String::trim).mapToInt(Integer::parseInt).toArray();
          } catch (final Exception ex) {
              throw new IllegalStateException("The GnarlyGenotyper tool assumes that input variants have SB FORMAT "
                      + " fields as a list of integers separated by commas.", ex);
          }
      } else {
          try {
              @SuppressWarnings("unchecked")
              final List<Integer> sbbsList = (ArrayList<Integer>) sbbsObj;
              return Ints.toArray(sbbsList);
          } catch (final ClassCastException e) {
              throw new IllegalStateException("The GnarlyGenotyper tool assumes that input variants have SB FORMAT " +
                      "fields parsed into ArrayLists.");
          }
      }
    }
}
