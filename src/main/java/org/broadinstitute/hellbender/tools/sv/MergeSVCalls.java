package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Creates sparsely formatted file of structural variants.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Standardized SV VCFs
 *     </li>
 *     <li>
 *         gCNV segments VCFs
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         TSV variants file
 *     </li>
 *     <li>
 *         Small CNV interval list (optional)
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk MergeSVCalls
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Creates sparse structural variants file",
        oneLineSummary = "Creates sparse structural variants file",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class MergeSVCalls extends GATKTool {
    public static final String MIN_GCNV_QUALITY_LONG_NAME = "min-gcnv-quality";
    public static final String SMALL_CNV_SIZE_LONG_NAME = "small-cnv-size";
    public static final String SMALL_CNV_PADDING_LONG_NAME = "small-cnv-padding";
    public static final String SMALL_CNV_OUTPUT_LONG_NAME = "small-cnv-output";
    public static final String IGNORE_DICTIONARY_LONG_NAME = "ignore-dict";

    public static final int DEFAULT_SMALL_CNV_SIZE = 5000;
    public static final int DEFAULT_SMALL_CNV_PADDING = 1000;

    @Argument(
            doc = "Input standardized SV and gCNV segments VCFs",
            fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME
    )
    private List<String> inputFiles;

    @Argument(
            doc = "Output TSV file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Argument(
            doc = "Small CNV intervals output",
            fullName = SMALL_CNV_OUTPUT_LONG_NAME,
            optional = true
    )
    private File smallCNVOutput;

    @Argument(
            doc = "Small CNV size",
            fullName = SMALL_CNV_SIZE_LONG_NAME,
            optional = true
    )
    private int smallCNVSize = DEFAULT_SMALL_CNV_SIZE;

    @Argument(
            doc = "Small CNV interval padding",
            fullName = SMALL_CNV_PADDING_LONG_NAME,
            optional = true
    )
    private int smallCNVPadding = DEFAULT_SMALL_CNV_PADDING;

    @Argument(
            doc = "Skip VCF sequence dictionary check",
            fullName = IGNORE_DICTIONARY_LONG_NAME,
            optional = true
    )
    private Boolean skipVcfDictionaryCheck = false;

    @Argument(
            doc = "Min gCNV quality (QS)",
            fullName = MIN_GCNV_QUALITY_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minGCNVQuality = 60;

    private List<SVCallRecord> records;
    private SAMSequenceDictionary dictionary;


    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        records = new ArrayList<>();
    }

    @Override
    public Object onTraversalSuccess() {
        records.sort(IntervalUtils.getDictionaryOrderComparator(dictionary));
        if (smallCNVOutput != null) {
            generateSmallEvents();
        }
        writeVariants();
        return null;
    }

    @Override
    public void traverse() {
        for (int i = 0; i < inputFiles.size(); i++) {
            processFile(inputFiles.get(i), i);
        }
    }

    private final void processFile(final String file, final int index) {
        final FeatureDataSource<VariantContext> source = getFeatureDataSource(file, "file_" + index);
        final VCFHeader header = getHeaderFromFeatureSource(source);
        final VCFHeaderLine headerSource = header.getMetaDataLine("source");
        final Stream<VariantContext> inputVariants = StreamSupport.stream(source.spliterator(), false);
        final Stream<SVCallRecord> inputRecords;
        if (headerSource != null && headerSource.getValue().equals(PostprocessGermlineCNVCalls.class.getSimpleName())) {
            inputRecords = inputVariants
                    .map(v -> SVCallRecord.createDepthOnlyFromGCNV(v, minGCNVQuality))
                    .filter(r -> r != null);
        } else {
            inputRecords = inputVariants.map(SVCallRecord::create);
        }
        final Collection<SVCallRecord> inputRecordCollection = inputRecords.collect(Collectors.toList());
        inputRecordCollection.stream().forEachOrdered(r -> progressMeter.update(r.getStartAsInterval()));
        records.addAll(inputRecordCollection);
    }

    private VCFHeader getHeaderFromFeatureSource(final FeatureDataSource<VariantContext> source) {
        final Object header = source.getHeader();
        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        return (VCFHeader)header;
    }

    private FeatureDataSource<VariantContext> getFeatureDataSource(final String vcf, final String name) {
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(
                vcf, name, 100000, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
        if (!skipVcfDictionaryCheck) {
            featureDataSource.getSequenceDictionary().assertSameDictionary(dictionary);
        }
        featureDataSource.setIntervalsForTraversal(getTraversalIntervals());
        return featureDataSource;
    }

    private void generateSmallEvents() {
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final List<SimpleInterval> smallEventIntervals = records.stream().filter(this::isSmallCNV)
                .map(r -> new SimpleInterval(r.getContig(), r.getStart(), r.getEnd())).collect(Collectors.toList());
        final List<GenomeLoc> smallEventLocs = IntervalUtils.genomeLocsFromLocatables(parser, smallEventIntervals);
        final GenomeLocSortedSet intervalSet = IntervalUtils.sortAndMergeIntervals(parser, smallEventLocs, IntervalMergingRule.ALL);
        final IntervalList intervalList = new IntervalList(dictionary);
        intervalList.addall(intervalSet.stream().map(Interval::new).collect(Collectors.toList()));
        intervalList.write(smallCNVOutput);
    }

    private static boolean isCNV(final SVCallRecord call) {
        return call.getType().equals(StructuralVariantType.DEL) || call.getType().equals(StructuralVariantType.DUP)
                || (call.getType().equals(StructuralVariantType.BND) && call.getContig().equals(call.getEndContig())
                    && call.getStartStrand() != call.getEndStrand());
    }

    private boolean isSmallCNV(final SVCallRecord call) {
        if (!isCNV(call)) return false;
        return call.getEnd() - call.getStart() <= smallCNVSize;
    }

    private void writeVariants() {
        final SVCallRecordCodec callRecordCodec = new SVCallRecordCodec();
        try (final PrintStream writer = new PrintStream(outputFile)) {
            records.stream().forEachOrdered(r -> writer.println(callRecordCodec.encode(r)));
        } catch(final IOException e) {
            throw new GATKException("Error writing output file", e);
        }
    }
}
