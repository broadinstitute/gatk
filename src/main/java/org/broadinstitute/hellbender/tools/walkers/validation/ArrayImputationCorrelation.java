package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.VariantIDsVariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import picard.sam.util.Pair;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.function.BiFunction;

@CommandLineProgramProperties(
        summary = ArrayImputationCorrelation.USAGE_SUMMARY,
        oneLineSummary = ArrayImputationCorrelation.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class ArrayImputationCorrelation extends AbstractConcordanceWalker {
    final static String USAGE_SUMMARY = "ArrayImputationCorrelation";
    final static String USAGE_ONE_LINE_SUMMARY = "ArrayImputationCorrelation";

    @Argument(fullName= StandardArgumentDefinitions.RESOURCE_LONG_NAME, doc="External resource VCF file", optional=true)
    public FeatureInput<VariantContext> af_resource;

    /**
     * List of IDs (or a .list file containing ids) to select. The tool will only select variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. If a file, the file
     * name must end in ".list", and the expected file format is simply plain text with one ID per line.
     */
    @Argument(fullName="keep-ids", shortName="ids", doc="List of variant rsIDs to select", optional=true)
    private Set<String> rsIDsToKeep = new HashSet<>();

    @Argument(fullName="dosage-field", optional=true)
    public String dosageField = null;

    @Argument(
            doc = "Output file for correlation.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile = null;

    @Argument(
            doc = "Output file for accuracy.",
            fullName = "output-accuracy",
            shortName = "OA"
    )
    private GATKPath outputAccuracyFile = null;

    @Argument(
        fullName = "af-annotations"
    )
    private List<String> afAnnotations;

    @Argument(
            fullName = "sample-map",
            optional = true
    )
    private List<String> sampleMappings;

    @Argument(
            optional = true
    )
    private int nBins = 14;

    @Argument(
            optional = true
    )
    private double firstBinRightEdge = 0.0005;

    @Argument(
            optional = true
    )
    private double minAfForAccuracyMetrics = 0;

    @Argument(optional = true)
    private boolean allowDifferingPloidy = false;

    final List<List<AFCorrelationAggregator>> aggregators = new ArrayList<>();
    final List<AccuracyMetrics> snpMetrics = new ArrayList<>();
    final List<AccuracyMetrics> indelMetrics = new ArrayList<>();
    final Map<String, String> afAnnotationsMap = new HashMap<>();
    final Set<String> afAnnotationSet = new HashSet<>();
    final Map<String, String> sampleMap = new HashMap<>();

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getEvalHeader();
        final List<String> samples = header.getSampleNamesInOrder();
        loadMapping(afAnnotations, afAnnotationsMap);
        afAnnotationSet.addAll(afAnnotationsMap.values());

        if (sampleMappings != null) {
            loadMapping(sampleMappings, sampleMap);
        }

        for (final String sample : samples) {
            String mappedSample = getMappedSample(sample);
            if (!afAnnotationsMap.containsKey(sample) && !afAnnotationsMap.containsKey(mappedSample)) {
                continue;
            }
            final List<AFCorrelationAggregator> theseAggregators = new ArrayList<>();
            theseAggregators.add(new AFCorrelationAggregator(firstBinRightEdge / 2, sample));
            for (int i = 1; i < nBins; i++) {
                final double log10_bin_width = (-Math.log10(firstBinRightEdge)) / (double) (nBins - 1);
                final double logBinCenter = Math.log10(firstBinRightEdge) + (i - 0.5) * log10_bin_width;
                final double binCenter = Math.pow(10, logBinCenter);
                theseAggregators.add(new AFCorrelationAggregator(binCenter, sample));
            }
            aggregators.add(theseAggregators);
            snpMetrics.add(new AccuracyMetrics(VariantContext.Type.SNP, sample));
            indelMetrics.add(new AccuracyMetrics(VariantContext.Type.INDEL, sample));
        }
    }

    private void loadMapping(final List<String> mappingStrings, final Map<String, String> map) {
        for (final String mappingString : mappingStrings) {
            String[] tokens = new String[2];
            StringUtil.split(mappingString, tokens,':');
            map.put(tokens[0], tokens[1]);
        }
    }

    private String getMappedSample(final String sample) {
        return sampleMap.getOrDefault(sample, sample);
    }

    private String getAfAnnotation(final String sample) {
        String afAnnotation = afAnnotationsMap.get(sample);
        if (afAnnotation == null) {
            afAnnotation = afAnnotationsMap.get(getMappedSample(sample));
        }

        return afAnnotation;
    }

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS::test;
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() {
        final VariantFilter filter;
        if (rsIDsToKeep != null && !rsIDsToKeep.isEmpty()) {
            filter = new VariantFilterLibrary.PassesFiltersVariantFilter().and(new VariantIDsVariantFilter(rsIDsToKeep));
        } else {
            filter = new VariantFilterLibrary.PassesFiltersVariantFilter();
        }
        return filter::test;
    }

    @Override
    protected boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval) {
        return truth.getAlleles().equals(eval.getAlleles());
    }

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {

        if (!truthVersusEval.hasEval() || !truthVersusEval.hasTruth()) {
            return;
        }

        final VariantContext evalVC = truthVersusEval.getEval();

        final VariantContext truthVC = truthVersusEval.getTruth();

        final Map<String, Double> afMap = buildAFMap(evalVC);

        for (int i=0; i<aggregators.size(); i++) {
            final String sample = snpMetrics.get(i).sampleName;
            final String mappedSample = getMappedSample(sample);
            final String afAnnotation = getAfAnnotation(sample);
            Double af = afMap.get(afAnnotation);
            if (af == null) {
                continue;
            }
            final int bin = getBin(af);
            final Genotype evalGenotype = evalVC.getGenotype(sample);
            final Genotype truthGenotype = truthVC.getGenotype(mappedSample);
            if (truthGenotype.isNoCall()) {
                continue;
            }
            final double evalRefFrac = getDosageFrac(evalGenotype, evalVC.getReference());
            final double truthRefFrac = getDosageFrac(truthGenotype, truthVC.getReference());
            final ConcordanceState concordanceState = getConcordanceState(truthGenotype, evalGenotype);
            if (evalVC.isSNP()) {
                aggregators.get(i).get(bin).snp_pearsonCorrelationAggregator.addEntry(evalRefFrac, truthRefFrac);
                if (af > minAfForAccuracyMetrics) {
                    snpMetrics.get(i).incrementMetrics(concordanceState);
                }
            } else if (evalVC.isIndel()) {
                aggregators.get(i).get(bin).indel_pearsonCorrelationAggregator.addEntry(evalRefFrac, truthRefFrac);
                if (af > minAfForAccuracyMetrics) {
                    indelMetrics.get(i).incrementMetrics(concordanceState);
                }
            }
        }
    }

    private Map<String, Double> buildAFMap(final VariantContext vc) {
        if (features != null) {
            final FeatureContext featureContext = new FeatureContext(features, new SimpleInterval(vc));
            final List<VariantContext> resourceFeatures = featureContext.getValues(af_resource, vc.getStart());
            for (final VariantContext feature : resourceFeatures) {
                if (feature.getReference().equals(vc.getReference())) {
                    final int altAlleleIndex = feature.getAlleleIndex(vc.getAlternateAllele(0));
                    if (altAlleleIndex >= 0) {
                        return buildAFMapForIndex(feature, altAlleleIndex);
                    }
                }
            }
            //haven't found a feature to extract afs from, so return empty map
            return Collections.emptyMap();
        } else {
            return buildAFMapForIndex(vc, 1);
        }
    }

    private Map<String, Double> buildAFMapForIndex(final VariantContext vc, final int index) {
        final Map<String, Double> afMap = new HashMap<>();
        for (final String afAnnotation : afAnnotationSet) {
            final List<Double> afList = vc.getAttributeAsDoubleList(afAnnotation, -99);
            if (afList.size() < index + 1) {
                continue;
            }
            final Double af = afList.get(index);
            if (af < 0 || af > 1) {
                throw new GATKException("Invalid AF value " + af + " at " + vc.getContig() + ":" + vc.getStart() + " for allele " + vc.getAlternateAllele(index-1));
            }
            afMap.put(afAnnotation, af);
        }

        return afMap;
    }


    private double getDosageFrac(final Genotype geno, final Allele refAllele) {
        if (dosageField != null) {
            final String dosageString = (String)geno.getExtendedAttribute(dosageField);
            return dosageString != null ? VCFUtils.parseVcfDouble(dosageString)/(double)geno.getPloidy() : getDosageFracFromGenotypeCall(geno, refAllele);
        }

        return getDosageFracFromGenotypeCall(geno, refAllele);
    }

    private double getDosageFracFromGenotypeCall(final Genotype geno, final Allele refAllele) {
        return 1d - (double)geno.countAllele(refAllele)/(double)geno.getPloidy();
    }

    private ConcordanceState getConcordanceState(final Genotype truth, final Genotype eval) {
        if (truth.getPloidy() != eval.getPloidy()) {
            //situation could arise from haploid vs diploid on X
            if (!allowDifferingPloidy) {
                throw new GATKException("sample " + eval.getSampleName() + " is ploidy " + eval.getPloidy() + " while truth sample " + truth.getSampleName() + " is ploidy " + truth.getPloidy() + "." +
                        "  This may be due to haploid vs diploid representation on X.  If you would like to allow for this type of data, use the allowDifferingPloidy argument.");
            }

            //build a set of all alleles in truth and eval genotypes.  If sets are the same, then true, if different, then false.  so 0/0 and 0 are concordant, 1/1 and 1 are concordant, but 0/1 and 1 are discordant.
            final Set<Allele> truthAlleles = new HashSet<>(truth.getAlleles());
            final Set<Allele> evalAlleles = new HashSet<>(eval.getAlleles());

            return evaluateConcordanceStateByAlleles(truthAlleles, evalAlleles);
        }

        //now normal situation, agree on ploidy.  we ignore phasing

        final List<Allele> sortedEvalAlleles = new ArrayList<>(eval.getAlleles());
        Collections.sort(sortedEvalAlleles);

        final List<Allele> sortedTruthAlleles = new ArrayList<>(eval.getAlleles());
        Collections.sort(sortedTruthAlleles);

        return evaluateConcordanceStateByAlleles(sortedTruthAlleles, sortedEvalAlleles);
    }

    private ConcordanceState evaluateConcordanceStateByAlleles(final Collection<Allele> truthAlleles, final Collection<Allele> evalAlleles) {
        if (truthAlleles.equals(evalAlleles)) {
            if (hasNonRefAlleles(evalAlleles)) {
                return ConcordanceState.TRUE_POSITIVE;
            } else {
                return ConcordanceState.TRUE_NEGATIVE;
            }
        } else {
            if (hasNonRefAlleles(evalAlleles)) {
                return ConcordanceState.FALSE_POSITIVE;
            } else {
                return ConcordanceState.FALSE_NEGATIVE;
            }
        }
    }

    private boolean hasNonRefAlleles(final Collection<Allele> alleles) {
        for (final Allele allele : alleles) {
            if (allele.isNonReference()) {
                return true;
            }
        }
         return false;
    }

    @Override
    public Object onTraversalSuccess() {
        try (final AFCorrelationWriter writer = new AFCorrelationWriter(outputFile.toPath())) {
            for (final List<AFCorrelationAggregator> theseAggregators : aggregators )
            writer.writeAllRecords(theseAggregators);
        } catch (final IOException ex) {
            throw new GATKException("error writing output", ex);
        }

        try (final AccuracyWriter writer = new AccuracyWriter(outputAccuracyFile.toPath())) {
            writer.writeAllRecords(snpMetrics);
            writer.writeAllRecords(indelMetrics);
        } catch (final IOException ex) {
            throw new GATKException("error writing output", ex);
        }
        return null;
    }


    private int getBin(double val) {
        if (val < firstBinRightEdge) {
            return 0;
        }

        return (int)Math.ceil((1-Math.log10(val)/Math.log10(firstBinRightEdge)) * (nBins - 1));
    }

    private static final class AFCorrelationAggregator {
        final double binCenter;
        final String sampleName;
        final PearsonCorrelationAggregator snp_pearsonCorrelationAggregator = new PearsonCorrelationAggregator();
        final PearsonCorrelationAggregator indel_pearsonCorrelationAggregator = new PearsonCorrelationAggregator();

        AFCorrelationAggregator(final double binCenter, final String sampleName) {
            this.binCenter = binCenter;
            this.sampleName = sampleName;
        }
    }

    private static class AFCorrelationWriter extends TableWriter<AFCorrelationAggregator> {
        AFCorrelationWriter(final Path file) throws IOException {
            super(file, new TableColumnCollection("bin_center", "sample_name", "snp_correlation", "snp_n_sites", "indel_correlation", "indel_n_sites"));
        }

        protected void composeLine(final AFCorrelationAggregator afCorrelationAggregator, final DataLine dataLine) {
            dataLine.set("bin_center", afCorrelationAggregator.binCenter)
                    .set("sample_name", afCorrelationAggregator.sampleName)
                    .set("snp_correlation", afCorrelationAggregator.snp_pearsonCorrelationAggregator.getCorrelation())
                    .set("snp_n_sites", afCorrelationAggregator.snp_pearsonCorrelationAggregator.n)
                    .set("indel_correlation", afCorrelationAggregator.indel_pearsonCorrelationAggregator.getCorrelation())
                    .set("indel_n_sites", afCorrelationAggregator.indel_pearsonCorrelationAggregator.n);
        }
    }

    static final class PearsonCorrelationAggregator {
        double sum_xy;
        double sum_x;
        double sum_y;
        double sum_x2;
        double sum_y2;
        int n;
        void addEntry(final double x, final double y) {
            sum_xy += x*y;
            sum_x += x;
            sum_y += y;
            sum_x2 += x*x;
            sum_y2 += y*y;
            n++;
        }

        double getCorrelation() {
            final double n_d = n;
            final double e_xy = sum_xy/n_d;
            final double e_x = sum_x/n_d;
            final double e_y = sum_y/n_d;
            final double e_x2 = sum_x2/n_d;
            final double e_y2 = sum_y2/n_d;

            final double r = (e_xy - e_x * e_y)/(Math.sqrt(e_x2 - e_x*e_x)* Math.sqrt(e_y2 - e_y*e_y));
            return r;
        }
    }

    private static class AccuracyWriter extends TableWriter<AccuracyMetrics> {
        AccuracyWriter(final Path file) throws IOException {
            super(file, new TableColumnCollection("type", "sample_name", "true_positives", "true_negatives", "false_positives", "false_negatives", "recall", "precision", "accuracy"));
        }

        protected void composeLine(final AccuracyMetrics metrics, final DataLine dataLine) {
            dataLine.set("type", metrics.type.toString())
                    .set("sample_name", metrics.sampleName)
                    .set("true_positives", metrics.true_positives)
                    .set("true_negatives", metrics.true_negatives)
                    .set("false_positives", metrics.false_positives)
                    .set("false_negatives", metrics.false_negatives)
                    .set("recall", metrics.getRecall())
                    .set("precision", metrics.getPrecision())
                    .set("accuracy", metrics.getAccuracy());
        }
    }

    static final class AccuracyMetrics extends MetricBase {
        final String sampleName;
        int true_positives;
        int true_negatives;
        int false_positives;
        int false_negatives;
        final VariantContext.Type type;

        AccuracyMetrics(final VariantContext.Type type, final String sampleName) {
            this.type = type;
            this.sampleName = sampleName;
        }
        double getRecall() {
            return (double)true_positives/((double)true_positives + (double) false_negatives);
        }

        double getPrecision() {
            return (double)true_positives/((double)true_positives + (double) false_positives);
        }

        double getAccuracy() {
            return (double)(true_positives + true_negatives)/(double)(true_positives + false_positives + true_negatives + false_negatives);
        }

        void incrementMetrics(final ConcordanceState concordanceState) {
            switch (concordanceState) {
                case TRUE_POSITIVE:
                    true_positives++;
                    break;
                case TRUE_NEGATIVE:
                case FILTERED_TRUE_NEGATIVE:
                    true_negatives++;
                    break;
                case FALSE_POSITIVE:
                    false_positives++;
                    break;
                case FALSE_NEGATIVE:
                case FILTERED_FALSE_NEGATIVE:
                    false_negatives++;
            }
        }
    }
}
