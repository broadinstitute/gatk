package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Created by tsato on 10/16/17.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class // TODO: check that this is correct
)
public class LearnHyperparameters extends CommandLineProgram {
    public static double DEFAULT_CONVERGENCE_THRESHOLD = 1e-2;
    public static int DEFAULT_MAX_ITERATIONS = 20;

    @Argument(fullName = CollectDataForReadOrientationFilter.REF_HISTOGRAM_TABLE_LONG_NAME,
            shortName = CollectDataForReadOrientationFilter.REF_HISTOGRAM_TABLE_SHORT_NAME,
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


    List<RefSiteHistogram> histograms;
    List<AltSiteRecord> altDesignMatrix;

    @Override
    protected void onStartup(){
        altDesignMatrix = AltSiteRecord.readAltSiteRecords(altDataTable);
        histograms = RefSiteHistogram.readRefSiteHistograms(refHistogramTable);
    }

    @Override
    public Object doWork(){
        final List<String> all3mers = SequenceUtil.generateAllKmers(3).stream().map(String::new).collect(Collectors.toList());
        int numRows = 0;

        try (BufferedReader reader = new BufferedReader(new FileReader(altDataTable));){
            // Count the number of lines (i.e. number of rows in the design matrix) in the hope that initializing the
            // array list to the exact size will give us a performance boost
            while (reader.readLine() != null) {
                numRows++;
            }
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while counting the number of lines in %s", altDataTable.toString()), e);
        }

        List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(altDataTable, numRows);
        Map<String, List<AltSiteRecord>> altDesignMatrixByContext = altDesignMatrix.stream().collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));

        List<Hyperparameters> hyperparameterEstimates = new ArrayList<>(64);
        for (final String refContext : all3mers){
            Optional<RefSiteHistogram> histogram = histograms.stream().filter(h -> h.getReferenceContext().equals(refContext)).findFirst();
            List<AltSiteRecord> altDesignMatrixForContext = altDesignMatrixByContext.get(refContext);

            if (! histogram.isPresent() || altDesignMatrixForContext == null){
                logger.info(String.format(String.format("Skipping the reference context %s, as it was not found in either ref or alt table", refContext)));
                continue;
            }

            final LearnHyperparametersEngine engine = new LearnHyperparametersEngine(histogram.get(), altDesignMatrixForContext,
                    converagenceThreshold, maxEMIterations);
            final Hyperparameters hyperparameters = engine.runEMAlgorithm(logger);
            hyperparameterEstimates.add(hyperparameters);
        }

        Hyperparameters.writeHyperparameters(hyperparameterEstimates, output);
        return "SUCCESS";
    }

    @Override
    protected void onShutdown() {}
}
