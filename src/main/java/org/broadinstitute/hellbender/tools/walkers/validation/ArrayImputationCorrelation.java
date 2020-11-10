package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
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
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.VariantIDsVariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = ArrayImputationCorrelation.USAGE_SUMMARY,
        oneLineSummary = ArrayImputationCorrelation.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class ArrayImputationCorrelation extends AbstractConcordanceWalker {
    final static String USAGE_SUMMARY = "ArrayImputationCorrelation";
    final static String USAGE_ONE_LINE_SUMMARY = "ArrayImputationCorrelation";

    @Argument(doc = "af annotation",
            shortName = "a")
    protected String afAnnotation;

    @Argument(fullName= StandardArgumentDefinitions.RESOURCE_LONG_NAME, doc="External resource VCF file", optional=true)
    public FeatureInput<VariantContext> af_resource;

    /**
     * List of IDs (or a .list file containing ids) to select. The tool will only select variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. If a file, the file
     * name must end in ".list", and the expected file format is simply plain text with one ID per line.
     */
    @Argument(fullName="keep-ids", shortName="ids", doc="List of variant rsIDs to select", optional=true)
    private Set<String> rsIDsToKeep = new HashSet<>();

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
            doc = "Output file for incorrectly genotyped sites",
            fullName = "output-incorrect-sites",
            shortName = "OI",
            optional = true
    )
    private GATKPath outputIncorrectSitesFile = null;

    final int nBins = 14;
    final double firstBinRightEdge = 0.0005;
    final List<AFCorrelationAggregator> aggregators = new ArrayList<>();
    final AccuracyMetrics snpMetrics = new AccuracyMetrics(VariantContext.Type.SNP);
    final AccuracyMetrics indelMetrics = new AccuracyMetrics(VariantContext.Type.INDEL);
    private VariantContext currentReferenceBlockVC;
    private IncorrectlyGenotypedSitesWriter incorrectlyGenotypedSitesWriter;

    @Override
    public void onTraversalStart() {
        aggregators.add(new AFCorrelationAggregator(firstBinRightEdge /2));
        for (int i = 1; i < nBins; i++) {
            final double log10_bin_width = (- Math.log10(firstBinRightEdge))/(double)(nBins - 1);
            final double logBinCenter = Math.log10(firstBinRightEdge) + (i - 0.5) * log10_bin_width;
            final double binCenter = Math.pow(10, logBinCenter);
            aggregators.add(new AFCorrelationAggregator(binCenter));
        }

        if (outputIncorrectSitesFile != null) {
            try {
                incorrectlyGenotypedSitesWriter = new IncorrectlyGenotypedSitesWriter(outputIncorrectSitesFile.toPath());
            }
            catch (final IOException ex) {
                throw new GATKException("Error writing out to file " + outputIncorrectSitesFile, ex);
            }
        }
    }

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        final VariantFilter filter = new VariantFilterLibrary.PassesFiltersVariantFilter();
        return filter::test;
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
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {

        if (!truthVersusEval.hasEval()) {
            return;
        }

        final VariantContext evalVC = truthVersusEval.getEval();

        final VariantContext truthVC;
        if (truthVersusEval.hasTruth() ) {
            truthVC = truthVersusEval.getTruth();
            if (truthVC.getGenotype(0).isHomRef()) {
                currentReferenceBlockVC = truthVC;
            } else {
                currentReferenceBlockVC = null;
            }
        } else if (currentReferenceBlockVC != null && currentReferenceBlockVC.overlaps(evalVC)) {
            truthVC = currentReferenceBlockVC;
        } else {
            return;
        }

        final FeatureContext featureContext = new FeatureContext(features, new SimpleInterval(evalVC));

        final List<VariantContext> resourceFeatures = featureContext.getValues(af_resource, evalVC.getStart());




        if ((truthVC.getReference().equals(evalVC.getReference()) && (truthVC.getAlternateAllele(0).equals(evalVC.getAlternateAllele(0)) || truthVC.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE))) ||
            truthVC == currentReferenceBlockVC) {

            Double af = null;
            for (final VariantContext vc : resourceFeatures) {
                if (vc.getReference().equals(evalVC.getReference())) {
                    for (int i = 0; i < vc.getNAlleles() - 1; i++) {
                        if (vc.getAlternateAllele(i).equals(evalVC.getAlternateAllele(0))) {
                            final List<Double> afList = vc.getAttributeAsDoubleList(afAnnotation, -99);
                            if (afList.isEmpty()) {
                                return;
                            }
                            af = afList.get(i);
                            break;
                        }
                    }
                }
            }
            if (af == null) {
                return;
            }
            if (af < 0 || af > 1) {
                throw new GATKException("invalid af value");
            }

            final int bin = getBin(af);
            final double evalRefFrac = (double)evalVC.getGenotype(0).countAllele(evalVC.getReference())/(double)evalVC.getGenotype(0).getPloidy();
            final double truthRefFrac = (double)truthVC.getGenotype(0).countAllele(truthVC.getReference())/(double)truthVC.getGenotype(0).getPloidy();
            final int truthAltCount = truthVC.getGenotype(0).getPloidy() - truthVC.getGenotype(0).countAllele(truthVC.getReference());
            final int evalAltCount = evalVC.getGenotype(0).getPloidy() - evalVC.getGenotype(0).countAllele(evalVC.getReference());
            if (evalVC.isSNP()) {
                aggregators.get(bin).snp_pearsonCorrelationAggregator.addEntry(evalRefFrac, truthRefFrac);
                snpMetrics.incrementMetrics(evalAltCount, truthAltCount);
            } else if (evalVC.isIndel()){
                aggregators.get(bin).indel_pearsonCorrelationAggregator.addEntry(evalRefFrac, truthRefFrac);
                indelMetrics.incrementMetrics(evalAltCount, truthAltCount);
            }

            if (incorrectlyGenotypedSitesWriter != null && truthAltCount != evalAltCount) {
                try {
                    incorrectlyGenotypedSitesWriter.writeRecord(evalVC);
                } catch (final IOException ex) {
                    throw new GATKException("error writing to incorrectly genotyped sites file", ex);
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try (final AFCorrelationWriter writer = new AFCorrelationWriter(outputFile.toPath())) {
            writer.writeAllRecords(aggregators);
        } catch (final IOException ex) {
            throw new GATKException("error writing output", ex);
        }

        try (final AccuracyWriter writer = new AccuracyWriter(outputAccuracyFile.toPath())) {
            writer.writeRecord(snpMetrics);
            writer.writeRecord(indelMetrics);
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

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        return true;
    }

    private static final class AFCorrelationAggregator {
        final double binCenter;
        final PearsonCorrelationAggregator snp_pearsonCorrelationAggregator = new PearsonCorrelationAggregator();
        final PearsonCorrelationAggregator indel_pearsonCorrelationAggregator = new PearsonCorrelationAggregator();

        AFCorrelationAggregator(final double binCenter) {
            this.binCenter = binCenter;
        }
    }

    private static class IncorrectlyGenotypedSitesWriter extends TableWriter<VariantContext> {
        IncorrectlyGenotypedSitesWriter(final Path file) throws IOException {
            super(file, new TableColumnCollection("chrom", "start", "end"));
        }

        protected void composeLine(final VariantContext vc, final DataLine dataLine) {
            dataLine.set("chrom", vc.getContig())
                    .set("start", vc.getStart())
                    .set("end", vc.getEnd());
        }
    }

    private static class AFCorrelationWriter extends TableWriter<AFCorrelationAggregator> {
        AFCorrelationWriter(final Path file) throws IOException {
            super(file, new TableColumnCollection("bin_center", "snp_correlation", "snp_n_sites", "indel_correlation", "indel_n_sites"));
        }

        protected void composeLine(final AFCorrelationAggregator afCorrelationAggregator, final DataLine dataLine) {
            dataLine.set("bin_center", afCorrelationAggregator.binCenter)
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
            super(file, new TableColumnCollection("type", "true_positives", "true_negatives", "false_positives", "false_negatives", "recall", "precision", "accuracy"));
        }

        protected void composeLine(final AccuracyMetrics metrics, final DataLine dataLine) {
            dataLine.set("type", metrics.type.toString())
                    .set("true_positives", metrics.true_positives)
                    .set("true_negatives", metrics.true_negatives)
                    .set("false_positives", metrics.false_positives)
                    .set("false_negatives", metrics.false_negatives)
                    .set("recall", metrics.getRecall())
                    .set("precision", metrics.getPrecision())
                    .set("accuracy", metrics.getAccuracy());
        }
    }

    static final class AccuracyMetrics {
        int true_positives;
        int true_negatives;
        int false_positives;
        int false_negatives;
        final VariantContext.Type type;

        AccuracyMetrics(final VariantContext.Type type) {
            this.type = type;
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

        void incrementMetrics(final int evalAltCount, final int truthAltCount) {
            if (evalAltCount == truthAltCount) {
                if (evalAltCount > 0) {
                    true_positives++;
                } else {
                    true_negatives++;
                }
            } else {
                if (evalAltCount > 0) {
                    false_positives++;
                } else {
                    false_negatives++;
                }
            }
        }
    }
}
