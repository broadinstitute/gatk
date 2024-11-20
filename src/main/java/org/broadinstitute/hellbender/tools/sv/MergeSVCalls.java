package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
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
 *     <li>
 *         cnMOPs call tables
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         VCF
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
@ExperimentalFeature
@DocumentedFeature
public final class MergeSVCalls extends GATKTool {
    public static final String MIN_GCNV_QUALITY_LONG_NAME = "min-gcnv-quality";
    public static final String SMALL_CNV_SIZE_LONG_NAME = "small-cnv-size";
    public static final String SMALL_CNV_PADDING_LONG_NAME = "small-cnv-padding";
    public static final String SMALL_CNV_OUTPUT_LONG_NAME = "small-cnv-output";
    public static final String IGNORE_DICTIONARY_LONG_NAME = "ignore-dict";
    public static final String CNMOPS_INPUT_LONG_NAME = "cnmops";
    public static final String COMPRESSION_LEVEL_LONG_NAME = "compression-level";
    public static final String CREATE_INDEX_LONG_NAME = "create-index";

    public static final int DEFAULT_SMALL_CNV_SIZE = 5000;
    public static final int DEFAULT_SMALL_CNV_PADDING = 1000;

    @Argument(
            doc = "Input standardized SV and gCNV segments VCFs",
            fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME
    )
    private List<String> inputFiles;

    @Argument(
            doc = "cnMOPS input files in tabular format",
            fullName = CNMOPS_INPUT_LONG_NAME
    )
    private List<String> cnmopsFiles;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

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
    private Set<String> samples;
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
        samples = records.stream().flatMap(r -> r.getCalledSamples().stream())
                .sorted().collect(Collectors.toCollection(LinkedHashSet::new));

        // Sort variants
        records = records.stream().sorted(SVCallRecordUtils.getCallComparator(dictionary)).collect(Collectors.toList());

        // Deduplicate
        final Function<Collection<SVCallRecord>, SVCallRecord> collapser = SVCallRecordUtils::deduplicateWithRawCallAttribute;
        final SVCallRecordDeduplicator<SVCallRecord> deduplicator = new SVCallRecordDeduplicator<>(collapser, dictionary);
        records = deduplicator.deduplicateItems(records);

        // Write
        writeVariants();
        return null;
    }

    @Override
    public void traverse() {
        for (final String path : cnmopsFiles) {
            processCNMOPSFile(path);
        }
        for (int i = 0; i < inputFiles.size(); i++) {
            processVcf(inputFiles.get(i), i);
        }
    }

    private void processCNMOPSFile(final String path) {
        try {
            final BufferedReader reader = new BufferedReader(IOUtils.makeReaderMaybeGzipped(Paths.get(path)));
            final String headerString = reader.readLine();
            if (headerString == null) {
                logger.warn("Skipping empty cnMOPS file: " + path);
                return;
            }
            logger.info("Parsing cnMOPS file: " + path);
            if (!headerString.startsWith("#")) {
                throw new UserException.BadInput("Expected first line to be a header starting with \"#\" but found: \"" + headerString + "\"");
            }
            reader.lines().map(this::cnmopsRecordParser).forEach(records::add);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Encountered exception while reading cnMOPS file: " + path, e);
        }
    }

    private SVCallRecord cnmopsRecordParser(final String line) {
        final String[] tokens = line.split("\t");
        final String contig = tokens[0];
        final int start = Integer.valueOf(tokens[1]);
        final int end = Integer.valueOf(tokens[2]);
        final String sample = tokens[4];
        final StructuralVariantType type = StructuralVariantType.valueOf(tokens[5]);
        if (!type.equals(StructuralVariantType.DEL) && !type.equals(StructuralVariantType.DUP)) {
            throw new UserException.BadInput("Unexpected cnMOPS record type: " + type.name());
        }
        final int length = end - start;
        final boolean isDel = type.equals(StructuralVariantType.DEL);
        final boolean startStrand = isDel;
        final boolean endStrand = !isDel;
        final List<String> algorithms = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);
        // Make genotypes het to play nicely with the defragmenter
        final List<Genotype> genotypes = Collections.singletonList(buildRawCallGenotype(sample));
        final String id = String.format("cnmops_%s_%s_%d_%d", type.toString(), contig, start, end);
        return new SVCallRecord(id, contig, start, startStrand, contig, end, endStrand, type, length, algorithms, genotypes);
    }

    private static Genotype buildRawCallGenotype(final String sample) {
        final GenotypeBuilder builder = new GenotypeBuilder(sample);
        builder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
        return builder.make();
    }

    private void processVcf(final String file, final int index) {
        final FeatureDataSource<VariantContext> source = getFeatureDataSource(file, "file_" + index);
        final VCFHeader header = getHeaderFromFeatureSource(source);
        final VCFHeaderLine headerSource = header.getMetaDataLine("source");
        final Stream<VariantContext> inputVariants = StreamSupport.stream(source.spliterator(), false);
        final Stream<SVCallRecord> inputRecords;
        if (headerSource != null && headerSource.getValue().equals(PostprocessGermlineCNVCalls.class.getSimpleName())) {
            inputRecords = inputVariants
                    .map(v -> SVCallRecordUtils.createDepthOnlyFromGCNVWithOriginalGenotypes(v, minGCNVQuality))
                    .filter(r -> r != null);
        } else {
            inputRecords = inputVariants.map(SVCallRecordUtils::create);
        }
        inputRecords.map(this::sanitizeVcfGenotypes).forEachOrdered(records::add);
    }

    private SVCallRecord sanitizeVcfGenotypes(final SVCallRecord record) {
        final Predicate<Genotype> filter = g -> SVCallRecord.isCarrier(g);
        final Function<Genotype, Map<String,Object>> attributeGenerator = g -> Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
        final GenotypesContext genotypes = SVCallRecordUtils.filterAndAddGenotypeAttributes(record.getGenotypes(), filter, attributeGenerator, true);
        return SVCallRecordUtils.copyCallWithNewGenotypes(record, genotypes);
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

    private void writeVariants() {
        final VariantContextWriter writer = createVCFWriter(outputFile);
        writer.writeHeader(getVcfHeader());
        records.stream().map(this::createVariant).forEachOrdered(writer::add);
        writer.close();
    }

    private VariantContext createVariant(final SVCallRecord call) {
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(call);
        builder.genotypes(SVCallRecordUtils.fillMissingSamplesWithEmptyGenotypes(builder.getGenotypes(), samples));
        builder.genotypes(sanitizeCreateVariantGenotypes(builder.getGenotypes()));
        return builder.make();
    }

    private GenotypesContext sanitizeCreateVariantGenotypes(final GenotypesContext genotypes) {
        final Predicate<Genotype> filter = g -> true;
        final Function<Genotype, Map<String,Object>> attributeGenerator =
                g -> g.hasAnyAttribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE) ? Collections.emptyMap() : Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
        return SVCallRecordUtils.filterAndAddGenotypeAttributes(genotypes, filter, attributeGenerator, false);
    }

    private VCFHeader getVcfHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRANDS_ATTRIBUTE, 1, VCFHeaderLineType.String, "First and second strands"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, 1, VCFHeaderLineType.String, "List of calling algorithms"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Sample non-reference in raw calls"));
        return header;
    }
}
