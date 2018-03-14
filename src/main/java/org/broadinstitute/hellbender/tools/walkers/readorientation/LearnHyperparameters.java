package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientationFilterConstants.*;

/**
 * Created by tsato on 10/16/17.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = CoverageAnalysisProgramGroup.class // TODO: check that this is correct
)
public class LearnHyperparameters extends CommandLineProgram {
    public static double DEFAULT_CONVERGENCE_THRESHOLD = 1e-4;
    public static int DEFAULT_MAX_ITERATIONS = 20;

    @Argument(fullName = CollectDataForReadOrientationFilter.REF_SITE_METRICS_LONG_NAME,
            shortName = CollectDataForReadOrientationFilter.REF_SITE_METRICS_SHORT_NAME,
            doc = "a tab-separated depth histogram over ref sites from CollectDataForReadOrientationFilter")
    private File refHistogramTable;

    @Argument(fullName = "alt-histogram",
            shortName = "alt-histogram",
            doc = "",
            optional = true)
    private File altHistogramTable = null;

    @Argument(fullName = CollectDataForReadOrientationFilter.ALT_DATA_TABLE_LONG_NAME,
            shortName = CollectDataForReadOrientationFilter.ALT_DATA_TABLE_SHORT_NAME,
            doc = "")
    private File altDataTable;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "table of hyperparameters")
    private File output;

    @Argument(fullName = "convergence-threshold",
            shortName = "c",
            doc = "convergence threshold for EM")
    private double converagenceThreshold = DEFAULT_CONVERGENCE_THRESHOLD;

    @Argument(fullName = "",
            shortName = "max",
            doc = "maximum iteration before we give up on EM")
    private int maxEMIterations = DEFAULT_MAX_ITERATIONS;


    List<Histogram<Integer>> refHistograms;
    List<Histogram<Integer>> altHistograms;

    @Override
    protected void onStartup(){
        final MetricsFile<?, Integer> referenceSiteMetrics = new MetricsFile<>();
        final Reader refHistogramReader = IOUtil.openFileForBufferedReading(refHistogramTable);
        referenceSiteMetrics.read(refHistogramReader);
        CloserUtil.close(refHistogramReader);

        final MetricsFile<?, Integer> altSiteMetrics = new MetricsFile<>();
        if (altHistogramTable != null){
            final Reader altHistogramReader = IOUtil.openFileForBufferedReading(altHistogramTable);
            altSiteMetrics.read(altHistogramReader);
            CloserUtil.close(altHistogramReader);
        }

        refHistograms = referenceSiteMetrics.getAllHistograms();
        altHistograms = altSiteMetrics.getAllHistograms();
    }

    @Override
    public Object doWork(){
        final int defaultInitialListSize = 1_000_000;

        final Map<String, List<AltSiteRecord>> altDesignMatrixByContext =
                AltSiteRecord.readAltSiteRecords(altDataTable, defaultInitialListSize)
                        .stream()
                        .collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));

        // Since AGT F1R2 is equivalent to ACT F1R2 (in the sense that the order of bases in the original molecule on which
        // the artifact befell and what the base changed to is the same)
        final List<Hyperparameters> hyperparametersAcrossContexts = new ArrayList<>((int) Math.pow(ReadOrientationFilterConstants.REGULAR_BASES.size(), ReadOrientationFilterConstants.REFERENCE_CONTEXT_SIZE)/2);

        // TODO: extract a method to account for reverse complement, and create a test
        for (final String refContext : ReadOrientationFilterConstants.CANONICAL_KMERS){
            final String reverseComplement = SequenceUtil.reverseComplement(refContext);

            final Histogram<Integer> refHistogram = refHistograms.stream()
                    .filter(h -> h.getValueLabel().equals(refContext))
                    .findFirst().orElseGet(() -> initializeRefHistogram(refContext));
            final Histogram<Integer> refHistogramRevComp = refHistograms.stream()
                    .filter(h -> h.getValueLabel().equals(reverseComplement))
                    .findFirst().orElseGet(() -> initializeRefHistogram(reverseComplement));

            final List<Histogram<Integer>> altHistogramsForContext = altHistograms.stream()
                    .filter(h -> h.getValueLabel().startsWith(refContext))
                    .collect(Collectors.toList());

            final List<Histogram<Integer>> altHistogramsRevComp = altHistograms.stream()
                    .filter(h -> h.getValueLabel().startsWith(reverseComplement))
                    .collect(Collectors.toList());

            final List<AltSiteRecord> altDesignMatrix = altDesignMatrixByContext
                    .getOrDefault(refContext, new ArrayList<>()); // Cannot use Collections.emptyList() here because we might add to it
            final List<AltSiteRecord> altDesignMatrixRevComp = altDesignMatrixByContext
                    .getOrDefault(reverseComplement, Collections.emptyList());
            // Warning: the below method will mutate the content of {@code altDesignMatrixRevComp} and append to {@code altDesignMatrix}
            mergeDesignMatrices(altDesignMatrix, altDesignMatrixRevComp);

            final Histogram<Integer> combinedRefHistograms = combineRefHistogramWithRC(refContext, refHistogram, refHistogramRevComp);
            final List<Histogram<Integer>> combinedAltHistograms = combineAltHistogramWithRC(altHistogramsForContext, altHistogramsRevComp);

            if (combinedRefHistograms.getSumOfValues() == 0 || altDesignMatrix.isEmpty()) {
                logger.info(String.format(
                        String.format("Skipping the reference context %s as we didn't find either the ref or alt table for the context", refContext)));
                continue;
            }

            final LearnHyperparametersEngine engine = new LearnHyperparametersEngine(
                    combinedRefHistograms,
                    combinedAltHistograms,
                    altDesignMatrix,
                    converagenceThreshold,
                    maxEMIterations,
                    logger);
            final Hyperparameters hyperparameters = engine.runEMAlgorithm();
            hyperparametersAcrossContexts.add(hyperparameters);
        }

        Hyperparameters.writeHyperparameters(hyperparametersAcrossContexts, output);
        return "SUCCESS";
    }

    @VisibleForTesting
    public static Histogram<Integer> combineRefHistogramWithRC(final String refContext,
                                                               final Histogram<Integer> refHistogram,
                                                               final Histogram<Integer> refHistogramRevComp){
        Utils.validateArg(refHistogram.getValueLabel()
                .equals(SequenceUtil.reverseComplement(refHistogramRevComp.getValueLabel())),
                "ref context = " + refHistogram.getValueLabel() + ", rev comp = " + refHistogramRevComp.getValueLabel());
        Utils.validateArg(refHistogram.getValueLabel().equals(refContext), "better match");

        final Histogram<Integer> combinedRefHistogram = new Histogram<>(binName, refContext);
        combinedRefHistogram.prefillBins(bins);

        for (final Integer depth : refHistogram.keySet()){
            final double newCount = refHistogram.get(depth).getValue() + refHistogramRevComp.get(depth).getValue();
            combinedRefHistogram.increment(depth, newCount);
        }

        return combinedRefHistogram;
    }

    @VisibleForTesting
    public static List<Histogram<Integer>> combineAltHistogramWithRC(final List<Histogram<Integer>> altHistograms,
                                                                     final List<Histogram<Integer>> altHistogramsRevComp){
        if (altHistograms.isEmpty() && altHistogramsRevComp.isEmpty()){
            return Collections.emptyList();
        }

        final String refContext = ! altHistograms.isEmpty() ?
                labelToContext(altHistograms.get(0).getValueLabel()).getLeft() :
                SequenceUtil.reverseComplement(labelToContext(altHistogramsRevComp.get(0).getValueLabel()).getLeft());

        // Contract: altHistogram must be of the canonical representaiton of the kmer
        Utils.validateArg(CANONICAL_KMERS.contains(refContext), "refContext must be the canonical representation but got " + refContext);

        final String reverseComplement = SequenceUtil.reverseComplement(refContext);
        final List<Histogram<Integer>> combinedHistograms = new ArrayList<>(numSubHistograms);


        for (Nucleotide altAllele : REGULAR_BASES){
            if (altAllele == Nucleotide.valueOf(refContext.substring(MIDDLE_INDEX, MIDDLE_INDEX+1))){
                continue;
            }

            final Nucleotide altAlleleRevComp = Nucleotide.valueOf(SequenceUtil.reverseComplement(altAllele.toString()));

            for (ArtifactType type : ArtifactType.values()) {
                final ArtifactType oppositeType = type == ArtifactType.F1R2 ? ArtifactType.F2R1 : ArtifactType.F1R2;
                final Histogram<Integer> altHistogram =
                        altHistograms.stream()
                                .filter(h -> h.getValueLabel().equals(contextToLabel(refContext, altAllele, type)))
                                .findFirst()
                                .orElseGet(() -> {
                                    final Histogram<Integer> h = new Histogram<>(binName, refContext);
                                    h.prefillBins(bins);
                                    return h;
                                });

                final Histogram<Integer> altHistogramRevComp =
                        altHistogramsRevComp.stream()
                                .filter(h -> h.getValueLabel().equals(contextToLabel(reverseComplement, altAlleleRevComp, oppositeType)))
                                .findFirst().orElseGet(() -> {
                                    final Histogram<Integer> h = new Histogram<>(binName, reverseComplement);
                                    h.prefillBins(bins);
                                    return h;
                                });

                final Histogram<Integer> combinedHistogram = new Histogram<>(binName, contextToLabel(refContext, altAllele, type));
                combinedHistogram.prefillBins(bins);

                // Add the histograms manually - I don't like the addHistogram() method because it does so with side-effect
                for (final Integer depth : altHistogram.keySet()){
                    final double newCount = altHistogram.get(depth).getValue() + altHistogramRevComp.get(depth).getValue();
                    combinedHistogram.increment(depth, newCount);
                }

                combinedHistograms.add(combinedHistogram);
            }
        }

        return combinedHistograms;
    }

    @VisibleForTesting
    public static void mergeDesignMatrices(final List<AltSiteRecord> altDesignMatrix, List<AltSiteRecord> altDesignMatrixRevComp){
        if (altDesignMatrix.isEmpty() && altDesignMatrixRevComp.isEmpty()){
            return;
        }

        // Order matters here. Assumes that all elements in the list have the same reference context
        Utils.validateArg(altDesignMatrix.isEmpty() ||
                CANONICAL_KMERS.contains(altDesignMatrix.get(0).getReferenceContext()),
                "altDesignMatrix must be in the parimary ref context...does that make sense?");

        final Optional<String> refContext = altDesignMatrix.isEmpty() ? Optional.empty() :
                Optional.of(altDesignMatrix.get(0).getReferenceContext());
        final Optional<String> revCompContext = altDesignMatrixRevComp.isEmpty() ? Optional.empty() :
                Optional.of(altDesignMatrixRevComp.get(0).getReferenceContext());
        if (refContext.isPresent() && revCompContext.isPresent()){
            // TOOD: write isReverseComplement(String context, String context2)
            Utils.validateArg(refContext.get().equals(SequenceUtil.reverseComplement(revCompContext.get())),
                    "context and rev comp don't match");
        }

        altDesignMatrix.addAll(
                altDesignMatrixRevComp.stream()
                        .map(AltSiteRecord::getReverseComplementOfRecord)
                        .collect(Collectors.toList())
        );
    }

}
