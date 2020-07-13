package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Retrieves SV evidence records on a given set of intervals.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Evidence file URI
 *     </li>
 *     <li>
 *         Interval list
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Evidence file (local)
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk LocalizeSVEvidence \
 *       --evidence-file gs://bucket/batch.SR.txt.gz \
 *       -L intervals.bed \
 *       -O local.SR.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Retrieves SV evidence records",
        oneLineSummary = "Retrieves SV evidence records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class LocalizeSVEvidence extends IntervalWalker {

    public static final String EVIDENCE_FILE_NAME = "evidence-file";
    public static final String INCLUDE_HEADER_STRING = "include-header";
    public static final String MAX_QUERY_SIZE_NAME = "max-query";

    @Argument(
            doc = "Input file URI with extension '.SR.txt.gz', '.PE.txt.gz', '.BAF.txt.gz', or '.bincov.bed.gz'",
            fullName = EVIDENCE_FILE_NAME
    )
    private String inputFilePath;

    @Argument(
            doc = "Output file. Note that files ending in '.gz' will NOT be block compressed.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFilePath;

    @Argument(
            doc = "Include header if it exists",
            fullName = INCLUDE_HEADER_STRING
    )
    private boolean includeHeader = false;

    @Advanced
    @Argument(
            doc = "Maximum query size, in bases. Lowering this can reduce memory usage but also increase run time.",
            fullName = MAX_QUERY_SIZE_NAME
    )
    private int maxQuerySize = 10000;

    private File outputFile;
    private PrintStream printStream;
    private FeatureDataSource<SVEvidence> source;

    private static final int QUERY_LOOKAHEAD = 0;

    @Override
    public void onTraversalStart() {
        outputFile = new File(outputFilePath);
        try {
            printStream = IOUtils.makePrintStreamMaybeGzipped(outputFile);
        }
        catch(IOException e) {
            throw new UserException.CouldNotCreateOutputFile(e.getMessage(), e);
        }
        initializeEvidenceDataSource();
    }

    private void initializeEvidenceDataSource() {
        source = new FeatureDataSource<>(
                inputFilePath,
                "inputFile",
                QUERY_LOOKAHEAD,
                SVEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
        doHeader(source);
    }

    private void doHeader(final FeatureDataSource<SVEvidence> source) {
        if (includeHeader) {
            final Object header = source.getHeader();
            if (header != null) {
                if (header instanceof String) {
                    printStream.println((String) header);
                } else {
                    throw new GATKException.ShouldNeverReachHereException("Expected header object of type " + String.class.getSimpleName());
                }
            } else {
                logger.warn("Header not found");
            }
        }
    }

    @Override
    public void apply(final SimpleInterval interval,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final List<SimpleInterval> intervals = partitionInterval(interval);
        for (final SimpleInterval i : intervals) {
            write(source.queryAndPrefetch(i));
        }
    }

    // Partitions interval into one or more intervals with given max length
    private List<SimpleInterval> partitionInterval(final SimpleInterval interval) {
        final int numIntervals = (int) Math.ceil(interval.getLengthOnReference() / (double) maxQuerySize);
        final List<SimpleInterval> intervals = new ArrayList<>(numIntervals);
        final String contig = interval.getContig();
        int start = interval.getStart();
        while (start < interval.getEnd()) {
            int end = Math.min(start + maxQuerySize, interval.getEnd());
            intervals.add(new SimpleInterval(contig, start, end));
            start = end;
        }
        return intervals;
    }

    private void write(final List<SVEvidence> data) {
        for (final SVEvidence d : data) {
            printStream.println(d.encode());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        printStream.close();
        return null;
    }
}
