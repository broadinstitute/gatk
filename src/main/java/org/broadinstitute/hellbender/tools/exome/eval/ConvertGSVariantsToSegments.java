package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.math.IntRange;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecord;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecordWriter;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegment;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

/**
 * Converts Genome STRiP (GS) variant call format (VCF) variants into {@link HiddenStateSegmentRecord} CopyNumberTriState format calls.
 *
 * <p>
 *     CopyNumberTriState refers to the three copy number calls of plus (+), neutral (0) and minus (-).
 *     Use the resulting file with {@link EvaluateCopyNumberTriStateCalls}.
 *     Here are example CopyNumberTriState format calls:
 * </p>
 *
 * <pre>
 *     Sample  Chromosome      Start   End     Segment_Call    Num_Targets     Segment_Mean    Segment_Std     Exact_Quality   Some_Quality    Start_Quality   End_Quality     Event_Quality
 *     01C05110        1       21944797        21946084        0       1       2.3443  0.0000  44.4002 NaN     NaN     NaN     0.0000
 *     01C06980        1       21944797        21946084        0       1       2.2533  0.0000  59.8000 NaN     NaN     NaN     0.0000
 *     01C07132        1       21944797        21946084        0       1       2.1878  0.0000  68.7000 NaN     NaN     NaN     0.0000
 *     01C07215        1       21944797        21946084        0       1       1.6267  0.0000  96.5991 NaN     NaN     NaN     0.0000
 *     01C07248        1       21944797        21946084        0       1       1.7863  0.0000  103.2934        NaN     NaN     NaN     0.0000
 *     01C07286        1       21944797        21946084        0       1       2.1658  0.0000  68.7000 NaN     NaN     NaN     0.0000
 *     01C07431        1       21944797        21946084        0       1       2.1904  0.0000  66.1000 NaN     NaN     NaN     0.0000
 *     01C07450        1       21944797        21946084        0       1       2.1317  0.0000  75.0000 NaN     NaN     NaN     0.0000
 *     01C05110        1       46597740        46598439        0       0       2.2678  0.0000  34.2017 NaN     NaN     NaN     0.0000
 *     01C06980        1       46597740        46598439        0       0       2.0544  0.0000  49.5000 NaN     NaN     NaN     0.0000
 *     01C07132        1       46597740        46598439        +       0       2.7620  0.0000  5.3133  NaN     NaN     NaN     5.3000
 *     01C07215        1       46597740        46598439        0       0       2.5568  0.0000  13.6898 NaN     NaN     NaN     0.2000
 * </pre>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" ConvertGSVariantsToSegments \
 *   --variant genomestrip.vcf.gz \
 *   --output copynumbertristate.seg
 * </pre>
 *
 * <p>
 *     Optional argument --targets (see {@link org.broadinstitute.hellbender.tools.exome.convertbed.ConvertBedToTargetFile})
 *     limits data to the listed target intervals.
 * </p>
 * 
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Convert Genome STRiP VCF file contents into a CopyNumberTriState segment file. Genome STRiP can be found at http://www.broadinstitute.org/software/genomestrip",
        oneLineSummary = "(Internal) Convert Genome STRiP VCF file into a segment file",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class ConvertGSVariantsToSegments extends VariantWalker {

    public static final String OUTPUT_FILE_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
    public static final String OUTPUT_FILE_FULL_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String NEUTRAL_COPY_NUMBER_SHORT_NAME = "neutral";
    public static final String NEUTRAL_COPY_NUMBER_FULL_NAME = "neutralCopyNumber";

    public static final int NEUTRAL_COPY_NUMBER_DEFAULT = 2;

    public static final String GS_COPY_NUMBER_FORMAT = "CN";
    public static final String GS_COPY_NUMBER_FRACTION = "CNF";
    public static final String GS_COPY_NUMBER_POSTERIOR = "CNP";

    @ArgumentCollection
    protected static final TargetArgumentCollection targetArguments = new TargetArgumentCollection();

    @Argument(
            doc = "Output segment file",
            shortName = OUTPUT_FILE_SHORT_NAME,
            fullName = OUTPUT_FILE_FULL_NAME
    )
    protected File outputFile;

    @Argument(
            doc = "Reference copy number",
            shortName = NEUTRAL_COPY_NUMBER_SHORT_NAME,
            fullName = NEUTRAL_COPY_NUMBER_FULL_NAME,
            optional = true
    )
    protected int neutralCopyNumber = NEUTRAL_COPY_NUMBER_DEFAULT;

    protected HiddenStateSegmentRecordWriter<CopyNumberTriState, Target> outputWriter;

    protected TargetCollection<Target> targets;

    @Override
    public void onTraversalStart() {
        try {
            outputWriter = new HiddenStateSegmentRecordWriter<>(outputFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
        targets = targetArguments.readTargetCollection(false);
    }

    private void closeOutputWriter(final TableWriter<?> writer, final File file, final List<RuntimeException> exceptions) {
        if (writer == null) {
            return;
        }
        try {
            writer.close();
        } catch (final Exception ex) {
            exceptions.add(new UserException.CouldNotCreateOutputFile(file, ex));
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final List<RuntimeException> exceptionList = new ArrayList<>(2);
        closeOutputWriter(outputWriter, outputFile, exceptionList);
        if (!exceptionList.isEmpty()) {
            for (final Exception ex : exceptionList) {
                logger.error(ex.getMessage());
            }
            throw exceptionList.get(0);
        }
        return "SUCCESS";
    }

    @Override
    public boolean requiresReference() {
        return false;
    }

    @Override
    public boolean requiresReads() {
        return false;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final SimpleInterval interval = new SimpleInterval(variant);
        final int targetCount = targets.indexRange(interval).size();
        final int[] callCounts = new int[CopyNumberTriState.values().length];
        for (final Genotype genotype : variant.getGenotypes().iterateInSampleNameOrder()) {
            final String sample = genotype.getSampleName();
            final double mean = doubleFrom(genotype.getExtendedAttribute(GS_COPY_NUMBER_FRACTION));
            final int copyNumber = intFrom(genotype.getExtendedAttribute(GS_COPY_NUMBER_FORMAT));
            final CopyNumberTriState call = copyNumber == neutralCopyNumber ? CopyNumberTriState.NEUTRAL : (copyNumber < neutralCopyNumber) ? CopyNumberTriState.DELETION : CopyNumberTriState.DUPLICATION;
            callCounts[call.ordinal()]++;
            final double[] probs = doubleArrayFrom(genotype.getExtendedAttribute(GS_COPY_NUMBER_POSTERIOR));
            final double log10PostQualCall = calculateLog10CallQuality(probs, call);
            final double log10PostQualNonRef = calculateLog10CallQualityNonRef(probs);
            final double phredProbCall = -10.0 * log10PostQualCall;
            final double phredProbNonRef = -10.0 * log10PostQualNonRef;
            final HiddenStateSegment<CopyNumberTriState, Target> segment = new HiddenStateSegment<>(
                    interval,
                    targetCount,
                    mean,
                    0.0, // GS VCF does not contain any stddev or var estimate for coverage fraction.
                    call,
                    phredProbCall, // GS does not provide an EQ, we approximate it to be the 1 - sum of all call compatible CN corresponding posterior probs
                    Double.NaN, // GS does not provide a SQ, we leave is a NaN.
                    Double.NaN, // GS does not provide a START Q.
                    Double.NaN, // GS does not provide a END Q.
                    phredProbNonRef
            );

            final HiddenStateSegmentRecord<CopyNumberTriState, Target> record = new HiddenStateSegmentRecord<>(sample, segment);
            try {
                outputWriter.writeRecord(record);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
            }
        }
    }

    private double calculateLog10CallQualityNonRef(final double[] log10Probs) {
        return log10Probs.length > neutralCopyNumber ? log10Probs[neutralCopyNumber] : -1000.0;
    }

    private double calculateLog10CallQuality(final double[] log10Probs, final CopyNumberTriState call) {
        final IntRange callCopyNumberRange;
        switch (call) {
            case NEUTRAL:
                callCopyNumberRange = new IntRange(neutralCopyNumber, neutralCopyNumber);
                break;
            case DELETION:
                callCopyNumberRange = new IntRange(0, neutralCopyNumber - 1);
                break;
            case DUPLICATION:
                callCopyNumberRange = new IntRange(neutralCopyNumber + 1, log10Probs.length - 1);
                break;
            default:
                throw new GATKException("unexpected call");
        }

        // We aggregate the probs of any copy number that that does not correspond
        final double log10OneMinusProbCall = MathUtils.approximateLog10SumLog10(
                MathUtils.log10SumLog10(log10Probs, 0, Math.min(callCopyNumberRange.getMinimumInteger(), log10Probs.length)),
                MathUtils.log10SumLog10(log10Probs, callCopyNumberRange.getMaximumInteger() + 1, log10Probs.length)
        );

        final double log10ProbTotal = MathUtils.log10SumLog10(log10Probs);
        return log10OneMinusProbCall - log10ProbTotal;
    }

    private double[] doubleArrayFrom(final Object object) {
        final String stringValue = String.valueOf(object);
        final String[] stringValues = stringValue.split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        return Stream.of(stringValues).mapToDouble(Double::valueOf).toArray();
    }

    private double doubleFrom(final Object object) {
        return Double.valueOf(String.valueOf(object));
    }

    private int intFrom(final Object object) {
        return Integer.valueOf(String.valueOf(object));
    }

}
