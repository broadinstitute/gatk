package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
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
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.function.Function;

/**
 * Retrieves SV evidence from a remote source on a given set of intervals.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Remote evidence file
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
 *         Evidence file
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk LocalizeSVEvidence --split-reads-file gs://bucket/batch.split.txt.gz -I intervals.bed -O local.split.txt.gz
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
    public static final String QUERY_LOOKAHEAD_NAME = "query-lookahead";

    @Argument(
            doc = "Input file URI.",
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
            doc = "Remote query lookahead, in bases.",
            fullName = QUERY_LOOKAHEAD_NAME
    )
    private int queryLookahead = 0;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private FeatureDataSource<DepthEvidence> depthSource;
    private FeatureDataSource<BafEvidence> bafSource;

    private File outputFile;
    private PrintStream printStream;

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
        if (inputFilePath.endsWith(SplitReadEvidenceCodec.FORMAT_SUFFIX)) {
            splitReadSource = new FeatureDataSource<>(
                    inputFilePath,
                    "inputFilePath",
                    queryLookahead,
                    SplitReadEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            doHeader(splitReadSource);
        } else if (inputFilePath.endsWith(DiscordantPairEvidenceCodec.FORMAT_SUFFIX)) {
            discordantPairSource = new FeatureDataSource<>(
                    inputFilePath,
                    "inputFilePath",
                    queryLookahead,
                    DiscordantPairEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            doHeader(discordantPairSource);
        } else if (inputFilePath.endsWith(DepthEvidenceCodec.FORMAT_SUFFIX)) {
            depthSource = new FeatureDataSource<>(
                    inputFilePath,
                    "inputFilePath",
                    queryLookahead,
                    DepthEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            doHeader(depthSource);
        } else if (inputFilePath.endsWith(BafEvidenceCodec.FORMAT_SUFFIX)) {
            bafSource = new FeatureDataSource<>(
                    inputFilePath,
                    "inputFilePath",
                    queryLookahead,
                    BafEvidence.class,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer);
            doHeader(bafSource);
        } else {
            throw new UserException.CouldNotReadInputFile("Input file suffix must be one of: " +
                String.join(
                        " ",
                        SplitReadEvidenceCodec.FORMAT_SUFFIX,
                        DiscordantPairEvidenceCodec.FORMAT_SUFFIX,
                        DepthEvidenceCodec.FORMAT_SUFFIX,
                        BafEvidenceCodec.FORMAT_SUFFIX
                )
            );
        }
    }

    private void doHeader(final FeatureDataSource<? extends Feature> source) {
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
        if (splitReadSource != null) {
            final List<SplitReadEvidence> data = splitReadSource.queryAndPrefetch(interval);
            write(data, SplitReadEvidenceCodec::encode);
        } else if (discordantPairSource != null) {
            final List<DiscordantPairEvidence> data = discordantPairSource.queryAndPrefetch(interval);
            write(data, DiscordantPairEvidenceCodec::encode);
        } else if (depthSource != null) {
            final List<DepthEvidence> data = depthSource.queryAndPrefetch(interval);
            write(data, DepthEvidenceCodec::encode);
        } else if (bafSource != null) {
            final List<BafEvidence> data = bafSource.queryAndPrefetch(interval);
            write(data, BafEvidenceCodec::encode);
        } else {
            throw new GATKException.ShouldNeverReachHereException("No data sources were initialized");
        }
    }

    final <T extends Feature> void write(final List<T> data, final Function<T, String> encoder) {
        for (final T d : data) {
            printStream.println(encoder.apply(d));
        }
    }

    @Override
    public Object onTraversalSuccess() {
        printStream.close();
        return null;
    }
}
