package org.broadinstitute.hellbender.tools.exome.gcbias;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Correct coverage for sample-specific GC bias effects as described in {@link GCCorrector}.
 *
 * <p>
 * Inputs are
 * 1. read counts file (format described in {@link ReadCountCollectionUtils}).  Since this tool corrects
 *    for multiplicative GC bias effects, these counts should not be in log space i.e. they should represent
 *    raw coverage or proportional coverage upstream of tangent normalization.
 * 2. targets file containing GC content annotation as produced by {@link AnnotateTargets}.  Every target in the input
 *    read counts must be present in the targets file but the reverse need not be true.
 *
 * The resulting read counts file has the same targets (rows) and samples (columns) as the input, with entries denoting GC-bias-corrected coverage.
 * Coverage is represented by doubles in {@link ReadCountCollection}.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CorrectGCBias \
 *   --input coverage.tsv \
 *   --targets targets.annotated.tsv \
 *   --output coverage.gc_corrected.tsv
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        oneLineSummary = "Correct for sample-specific GC bias",
        summary = "Correct coverage in read counts by (i) estimating per-sample bias as a function of per-target GC " +
                "content and (ii) dividing input coverage by derived bias curves.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public class CorrectGCBias extends CommandLineProgram {
    public static final String INPUT_READ_COUNTS_FILE_LONG_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String INPUT_READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String OUTPUT_READ_COUNTS_FILE_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String OUTPUT_READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    @Argument(
            doc = "Read counts input file.",
            shortName = INPUT_READ_COUNTS_FILE_SHORT_NAME,
            fullName = INPUT_READ_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected File inputReadCountsFile;

    @Argument(
            doc = "Output file of GC-corrected read counts.",
            shortName = OUTPUT_READ_COUNTS_FILE_SHORT_NAME,
            fullName = OUTPUT_READ_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected File outputReadCountsFile;

    // arguments for GC-annotated targets
    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection();

    @Override
    protected Object doWork() {
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);
        final ReadCountCollection inputCounts = readInputCounts(targets);
        final double[] gcContentByTarget = gcContentsOfTargets(inputCounts, targets);
        final ReadCountCollection outputCounts = GCCorrector.correctCoverage(inputCounts, gcContentByTarget);
        writeOutputCounts(outputCounts);
        return "SUCCESS";
    }

    private ReadCountCollection readInputCounts(final TargetCollection<Target> targets) {
        final ReadCountCollection inputCounts;
        try {
            inputCounts = ReadCountCollectionUtils.parse(inputReadCountsFile, targets, false);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputReadCountsFile, ex.getMessage(), ex);
        }
        return inputCounts;
    }

    private void writeOutputCounts(final ReadCountCollection outputCounts) {
        try {
            ReadCountCollectionUtils.write(outputReadCountsFile, outputCounts);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputReadCountsFile, ex);
        }
    }

    private double[] gcContentsOfTargets(final ReadCountCollection inputCounts, final TargetCollection<Target> targets) {
        final List<Target> annotatedTargets = inputCounts.targets().stream().map(t -> targets.target(t.getName())).collect(Collectors.toList());
        if (!annotatedTargets.stream().allMatch(t -> t.getAnnotations().hasAnnotation(TargetAnnotation.GC_CONTENT))) {
            throw new UserException.BadInput("At least one target lacks a GC annotation.");
        }
        return annotatedTargets.stream().mapToDouble(t -> t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT)).toArray();
    }


}
