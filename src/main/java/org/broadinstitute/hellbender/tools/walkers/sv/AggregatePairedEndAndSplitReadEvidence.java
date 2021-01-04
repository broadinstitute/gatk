package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Retrieves PE/SR evidence and performs breakpoint refinement
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 *     <li>
 *         Mean depth table
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk AggregatePairedEndAndSplitReadEvidence
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Adds PE/SR evidence to structural variant records",
        oneLineSummary = "Adds PE/SR evidence to structural variant records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class AggregatePairedEndAndSplitReadEvidence extends VariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
    public static final String PE_INNER_WINDOW_LONG_NAME = "pe-inner-window";
    public static final String PE_OUTER_WINDOW_LONG_NAME = "pe-outer-window";
    public static final String SR_WINDOW_LONG_NAME = "sr-window";

    @Argument(
            doc = "Split reads evidence file",
            fullName = SPLIT_READ_LONG_NAME,
            optional = true
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs evidence file",
            fullName = DISCORDANT_PAIRS_LONG_NAME,
            optional = true
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Tab-delimited table with sample IDs in the first column and expected per-base coverage per sample " +
                    "in the second column. Required if a split read evidence file is provided.",
            fullName = SAMPLE_COVERAGE_LONG_NAME,
            optional = true
    )
    private GATKPath sampleCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Inner discordant pair window size (bp)",
            fullName = PE_INNER_WINDOW_LONG_NAME,
            optional = true
    )
    private int innerWindow = 50;

    @Argument(
            doc = "Outer discordant pair window size (bp)",
            fullName = PE_OUTER_WINDOW_LONG_NAME,
            optional = true
    )
    private int outerWindow = 500;

    @Argument(
            doc = "Split read window size (bp)",
            fullName = SR_WINDOW_LONG_NAME,
            optional = true
    )
    private int splitReadWindow = 200;

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private BreakpointRefiner breakpointRefiner;
    private CachingSVEvidenceAggregator<DiscordantPairEvidence> discordantPairCollector;
    private CachingSVEvidenceAggregator<SplitReadEvidence> startSplitCollector;
    private CachingSVEvidenceAggregator<SplitReadEvidence> endSplitCollector;
    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;
    private VCFHeader header;

    private List<SVCallRecordWithEvidence> records;

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());

        if (!splitReadCollectionEnabled() && !discordantPairCollectionEnabled()) {
            throw new UserException.BadInput("At least one evidence file must be provided");
        }
        if (discordantPairCollectionEnabled()) {
            initializeDiscordantPairCollection();
        }
        if (splitReadCollectionEnabled()) {
            initializeSplitReadCollection();
        }
        records = new ArrayList<>();
        writer = createVCFWriter(Paths.get(outputFile));
        header = getVCFHeader();
        writer.writeHeader(header);
    }

    private void initializeDiscordantPairCollection() {
        initializeDiscordantPairDataSource();
        discordantPairCollector = new DiscordantPairEvidenceAggregator(discordantPairSource, dictionary, innerWindow, outerWindow);
    }

    private void initializeSplitReadCollection() {
        if (sampleCoverageFile == null) {
            throw new UserException.BadInput("Sample coverage table is required when a split read evidence file is present");
        }
        initializeSplitReadEvidenceDataSource();
        loadSampleCoverage();
        startSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, true);
        endSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, false);
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap, dictionary);
    }

    private void initializeDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile.toString(),
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile.toString(),
                "splitReadsFile",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadSampleCoverage() {
        final String fileString = sampleCoverageFile.toString();
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(fileString, e);
        }
    }

    @Override
    public void closeTool() {
        if (discordantPairSource != null) {
            discordantPairSource.close();
        }
        if (splitReadSource != null) {
            splitReadSource.close();
        }
        if (writer != null) {
            writer.close();
        }
        super.closeTool();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        records.add(new SVCallRecordWithEvidence(SVCallRecordUtils.create(variant, true)));
    }

    @Override
    public Object onTraversalSuccess() {
        List<SVCallRecordWithEvidence> currentRecords = records;

        // Get PE evidence
        if (discordantPairCollectionEnabled()) {
            discordantPairCollector.startProgressMeter();
            currentRecords = discordantPairCollector
                    .collectEvidence(currentRecords)
                    .collect(Collectors.toList());
            discordantPairCollector.stopProgressMeter();
        }

        // Get SR evidence
        if (splitReadCollectionEnabled()) {
            startSplitCollector.startProgressMeter();
            // Sort by end position to minimize cache misses for end site queries
            currentRecords = startSplitCollector.collectEvidence(currentRecords)
                    .sorted(SVCallRecordUtils.getSVLocatableComparatorByEnds(dictionary))
                    .collect(Collectors.toList());
            startSplitCollector.stopProgressMeter();
            endSplitCollector.startProgressMeter();
            // Refine breakpoints and clean up
            currentRecords = endSplitCollector.collectEvidence(currentRecords).map(r -> breakpointRefiner.refineCall(r))
                    .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                    .collect(Collectors.toList());
            endSplitCollector.stopProgressMeter();
        }

        // Write
        write(currentRecords);
        return super.onTraversalSuccess();
    }

    private void write(final List<SVCallRecordWithEvidence> records) {
        records.stream().map(this::buildVariantContext).forEachOrdered(writer::add);
    }

    private VCFHeader getVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        if (splitReadCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
        }
        if (discordantPairCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
        }
        return header;
    }

    private boolean splitReadCollectionEnabled() {
        return splitReadsFile != null;
    }

    private boolean discordantPairCollectionEnabled() {
        return discordantPairsFile != null;
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        return SVCallRecordUtils.createBuilderWithEvidence(call, splitReadCollectionEnabled(), discordantPairCollectionEnabled()).make();
    }
}
