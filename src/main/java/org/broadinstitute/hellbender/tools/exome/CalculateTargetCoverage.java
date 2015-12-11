package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TargetWalker;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Calculate the coverage (integer read counts) of targets.
 * <p>
 * The user must indicate the target table file using the argument -{@value TargetArgumentCollection#TARGET_FILE_SHORT_NAME}.
 * That file must follow the format described in {@link TargetTableReader}.
 * </p>
 * <p>
 * Alignment files must also be provided using -{@value StandardArgumentDefinitions#INPUT_SHORT_NAME}
 * </p>
 *
 * <p>
 *     The output follows the same format as the input target file including an additional column per each sample
 *     with the total number of reads that overlap each target by at least one base.
 * </p>
 *
 * <p>
 *     A read must match the following criteria to be taken into account:
 *     <ul>
 *         <li>must be <i>well-formed</i> as defined in {@link WellformedReadFilter},</li>
 *         <li>must be mapped (unmapped reads with mapped mates are not taken into account),</li>
 *         <li>must align with at least one base in the reference (all insertion reads are discarded) and</li>
 *         <li>must not be marked as duplicated.</li>
 *     </ul>
 * </p>
 *
 * <p>
 *   By default, all targets in the inputs will be present in the output. However the user can indicate
 *   the subset of targets to output using the {@value TargetArgumentCollection#TARGET_FILE_SHORT_NAME} argument.
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        oneLineSummary = "Calculate Target Coverage",
        summary = "Calculate Target Coverage",
        programGroup = ExomeAnalysisProgramGroup.class
)
public class CalculateTargetCoverage extends TargetWalker {

    /**
     * Short name for the {@link #maximumCoverage} argument.
     */
    public static final String MAXIMUM_COVERAGE_SHORT_NAME = "max";

    /**
     * Long name for the {@link #maximumCoverage} argument.
     */
    public static final String MAXIMUM_COVERAGE_FULL_NAME = "maximumCoverage";

    /**
     * Default value for the {@link #maximumCoverage} argument.
     */
    public static final long MAXIMUM_COVERAGE_DEFAULT = Long.MAX_VALUE;

    /**
     * Short name for the {@link #minimumMappingQuality} argument.
     */
    public static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";

    /**
     * Long name for the {@link #minimumMappingQuality} argument.
     */
    public static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";

    /**
     * Default value for the {@link #minimumMappingQuality} argument.
     */
    public static final int MINIMUM_MAPPING_QUALITY_DEFAULT = 0;

    @Argument(
            doc = "Output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional  = false
    )
    protected File outputFile;

    @Argument(
            doc = "Maximum coverage per target and coverage group",
            shortName = MAXIMUM_COVERAGE_SHORT_NAME,
            fullName = MAXIMUM_COVERAGE_FULL_NAME,
            optional = true
    )
    protected long maximumCoverage = MAXIMUM_COVERAGE_DEFAULT;

    @Argument(
            doc = "Minimum mapping quality",
            shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            fullName  = MINIMUM_MAPPING_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumMappingQuality = MINIMUM_MAPPING_QUALITY_DEFAULT;

    protected TableWriter<ReadCountRecord> outputTableWriter;

    protected ToIntFunction<GATKRead> readToColumn;

    protected CountingReadFilter readFilter;

    protected int countColumnCount = 0;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        final List<String> sampleList = getHeaderForReads().getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .filter(Objects::nonNull)
                .sorted()
                .distinct()
                .collect(Collectors.toList());
        final Map<String, Integer> readGroupToIndex = getHeaderForReads().getReadGroups().stream()
                .filter(rg -> rg.getSample() != null)
                .collect(Collectors.toMap(SAMReadGroupRecord::getId, rg -> sampleList.indexOf(rg.getSample())));

        readToColumn = (read) -> readGroupToIndex.getOrDefault(read.getReadGroup(), -1);
        countColumnCount = sampleList.size();
        readFilter = makeReadFilter();

        try  {
            final Writer outputWriter = createOutputWriter(outputFile);
            outputTableWriter = ReadCountCollectionUtils.writerWithIntervals(outputWriter, sampleList);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    @Override
    public Object onTraversalDone() {
        try {
            outputTableWriter.close();
            return null;
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problem closing the output:"  + ex.getMessage());
        } finally {
            outputTableWriter = null;
        }
    }

    private static Writer createOutputWriter(final File outputFile) {
        try {
            return new FileWriter(outputFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "could not open output file", ex);
        }
    }

    private CountingReadFilter makeReadFilter() {
        final CountingReadFilter baseFilter = new CountingReadFilter("Wellformed", new WellformedReadFilter(getHeaderForReads()))
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED))
                .and(new CountingReadFilter("Not_Duplicate", ReadFilterLibrary.NOT_DUPLICATE))
                .and(new CountingReadFilter("Non_Zero_Reference_Length", ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT));

        return minimumMappingQuality <= 0 ? baseFilter :
                baseFilter.and(new CountingReadFilter("MinMQ_" + minimumMappingQuality, read -> read.getMappingQuality() >= minimumMappingQuality));
    }

    @Override
    public void apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final long[] counts = new long[countColumnCount];
        StreamSupport.stream(readsContext.spliterator(), false)
                .filter(readFilter)
                .mapToInt(readToColumn)
                .filter(i -> i >= 0)
                .forEach(i -> counts[i]++);
        // cap by the maximum coverage allowed.
        for (int i = 0; i < counts.length; i++) {
            counts[i] = Math.min(counts[i], maximumCoverage);
        }
        try {
            outputTableWriter.writeRecord(new ReadCountRecord(target, counts));
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }
}
