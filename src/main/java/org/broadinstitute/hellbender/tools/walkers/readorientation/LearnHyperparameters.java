package org.broadinstitute.hellbender.tools.walkers.readorientation;

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
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.CollectDataForReadOrientationFilter.*;

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


    List<Histogram<Integer>> histograms;
    List<AltSiteRecord> altDesignMatrix;

    @Override
    protected void onStartup(){
        altDesignMatrix = AltSiteRecord.readAltSiteRecords(altDataTable);

        final MetricsFile<?, Integer> referenceSiteMetrics = new MetricsFile<>();
        final Reader in = IOUtil.openFileForBufferedReading(refHistogramTable);
        referenceSiteMetrics.read(in);
        CloserUtil.close(in);

        histograms = referenceSiteMetrics.getAllHistograms();
    }

    @Override
    public Object doWork(){
        final int defaultInitialListSize = 1_000_000;

        final List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(altDataTable, defaultInitialListSize);
        final Map<String, List<AltSiteRecord>> altDesignMatrixByContext = altDesignMatrix.stream()
                .collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));

        // Since AGT F1R2 is equivalent to ACT F1R2 (in the sense that the order of bases in the original molecule on which
        // the artifact befell and what the base changed to is the same)
        final List<Hyperparameters> hyperparametersAcrossContexts = new ArrayList<>((int) Math.pow(REGULAR_BASES.size(), REFERENCE_CONTEXT_SIZE)/2);

        for (final String refContext : ALL_KMERS_MODULO_REVERSE_COMPLEMENT){
            final String reverseComplement = SequenceUtil.reverseComplement(refContext);

            // Contract: {@code CollectDataForReadOrientationFilter} outputs the ref histogram for all 4^K contexts,
            // even when some contexts did not appear in the bam
            final Histogram<Integer> histogram = histograms.stream()
                    .filter(h -> h.getValueLabel().equals(refContext))
                    .findFirst().get();
            final Histogram<Integer> histogramRevComp = histograms.stream()
                    .filter(h -> h.getValueLabel().equals(reverseComplement))
                    .findFirst().get();
            histogram.addHistogram(histogramRevComp);

            List<AltSiteRecord> altDesignMatrixForContext = altDesignMatrixByContext.getOrDefault(refContext, new ArrayList<>());
            List<AltSiteRecord> altDesignMatrixForRevComp = altDesignMatrixByContext.getOrDefault(reverseComplement, Collections.emptyList())
                    .stream()
                    .map(AltSiteRecord::getReverseComplementOfRecord)
                    .collect(Collectors.toList());
            altDesignMatrixForContext.addAll(altDesignMatrixForRevComp);

            if (histogram.getCount() == 0 || altDesignMatrixForContext.isEmpty()) {
                logger.info(String.format(
                        String.format("Skipping the reference context %s as we didn't find either the ref or alt table for the context", refContext)));
                continue;
            }

            final LearnHyperparametersEngine engine = new LearnHyperparametersEngine(
                    histogram,
                    altDesignMatrixForContext,
                    converagenceThreshold,
                    maxEMIterations,
                    logger);
            final Hyperparameters hyperparameters = engine.runEMAlgorithm();
            hyperparametersAcrossContexts.add(hyperparameters);
        }

        Hyperparameters.writeHyperparameters(hyperparametersAcrossContexts, output);
        return "SUCCESS";
    }

    @Override
    protected void onShutdown() {}
}
