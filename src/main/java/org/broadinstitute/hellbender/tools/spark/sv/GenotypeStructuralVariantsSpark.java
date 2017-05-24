package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.avro.mapred.Pair;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ArrayUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * Created by valentin on 4/20/17.
 */
@CommandLineProgramProperties(summary = "genotype SV variant call files",
        oneLineSummary = "genotype SV variant call files",
        programGroup = StructuralVariationSparkProgramGroup.class)
public class GenotypeStructuralVariantsSpark extends GATKSparkTool {


    public static final String FASTQ_FILE_DIR_SHORT_NAME = "fastqDir";
    public static final String FASTQ_FILE_DIR_FULL_NAME = "fastqAssemblyDirectory";
    public static final String FASTQ_FILE_NAME_PATTERN_SHORT_NAME = "fastqName";
    public static final String FASTQ_FILE_NAME_PATTERN_FULL_NAME = "fastqNameFormat";
    public static final String PADDING_SHORT_NAME = "padding";
    public static final String PADDING_FULL_NAME = "padding";
    public static final String INSERT_SIZE_DISTR_SHORT_NAME = "insSize";
    public static final String INSERT_SIZE_DISTR_FULL_NAME = "insertSizeDistribution";

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private RequiredVariantInputArgumentCollection variantArguments = new RequiredVariantInputArgumentCollection();

    @Argument(doc = "fastq files location",
            shortName = FASTQ_FILE_DIR_SHORT_NAME,
            fullName  = FASTQ_FILE_DIR_FULL_NAME)
    private File fastqDir = null;

    @Argument(doc = "fastq files name pattern",
            shortName = FASTQ_FILE_NAME_PATTERN_SHORT_NAME,
            fullName  = FASTQ_FILE_NAME_PATTERN_FULL_NAME)
    private String fastqFilePattern = "assembly%2d.fastq";

    @Argument(doc = "output VCF file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    private File outputFile = null;

    @Argument(doc = "padding",
            shortName = PADDING_SHORT_NAME,
            fullName = PADDING_FULL_NAME,
            optional = true)
    private int padding = 300;

    @Argument(doc = "insert size distribution",
              shortName = INSERT_SIZE_DISTR_SHORT_NAME,
              fullName = INSERT_SIZE_DISTR_FULL_NAME,
              optional = true)
    private InsertSizeDistribution dist = new InsertSizeDistribution("N(300,100)");

    @Argument(doc = "pair-hmm implementation")
    private StructuralVariantPairHMMImplementation pairHmm = new StructuralVariantPairHMMImplementation("Affine(45,10)");

    private VariantsSparkSource variantsSource;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private void setUp(final JavaSparkContext ctx) {
        if (!outputFile.getParentFile().isDirectory()) {
            throw new CommandLineException.BadArgumentValue(StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                    "the output file location is not a directory:");
        } else if (outputFile.exists() && !outputFile.isFile()) {
            throw new CommandLineException.BadArgumentValue(StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                    "the output file makes reference to something that is not a file");
        }
        variantsSource = new VariantsSparkSource(ctx);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        setUp(ctx);
        final JavaRDD<StructuralVariantContext> variants = variantsSource.getParallelVariantContexts(
                variantArguments.variantFiles.get(0).getFeaturePath(), getIntervals())
                .map(StructuralVariantContext::new);
        final JavaRDD<StructuralVariantContext> outputVariants = processVariants(variants, ctx);
        final VCFHeader header = composeOutputHeader();
        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputFile.getParent(), outputFile.getName(),
                referenceArguments.getReferenceFile().getAbsolutePath(), outputVariants, header, logger);
        tearDown(ctx);
    }

    private VCFHeader composeOutputHeader() {
        final SAMFileHeader readHeader = getHeaderForReads();
        final List<String> samples = readHeader.getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .distinct()
                .sorted()
                .collect(Collectors.toList());
        final VCFHeader result = new VCFHeader(Collections.emptySet(), samples);
        result.setSequenceDictionary(getReferenceSequenceDictionary());
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "last base position of the variant"));
        return result;
    }

    private JavaRDD<StructuralVariantContext> processVariants(final JavaRDD<StructuralVariantContext> variants, final JavaSparkContext ctx) {
        final JavaPairRDD<StructuralVariantContext, List<Integer>> variantAndAssemblyIds = variants.mapToPair(v -> new Tuple2<>(v, v.assemblyIDs()));
        final File fastqDir = this.fastqDir;
        final String fastqFilePattern = this.fastqFilePattern;

        final JavaPairRDD<StructuralVariantContext, List<File>> variantAndAssemblyFiles = variantAndAssemblyIds
                .mapValues(v -> v.stream().map(id -> new File(fastqDir, String.format(fastqFilePattern, id))).collect(Collectors.toList()));
        final JavaPairRDD<StructuralVariantContext, List<SVFastqUtils.FastqRead>> variantReads = variantAndAssemblyFiles
                .mapValues(v -> v.stream().flatMap(file -> SVFastqUtils.readFastqFile(file.getAbsolutePath(), null).stream()).collect(Collectors.toList()));
        final JavaPairRDD<StructuralVariantContext, Map<String, List<SVFastqUtils.FastqRead>>> variantReadsByName = variantReads
                .mapValues(v -> v.stream().collect(Collectors.groupingBy(r -> templateName(r)._1())));
        final JavaPairRDD<StructuralVariantContext, List<Template>> variantTemplates = variantReadsByName.mapValues(map ->
            map.entrySet().stream().map(entry ->
                Template.create(entry.getKey(), entry.getValue(), read -> new Template.Fragment(read.getBases(), ArrayUtils.toInts(read.getQuals(), false)))
            ).collect(Collectors.toList()));
        variantTemplates.foreach(vt -> {
            if (structuralVariantAlleleIsSupported(vt._1().getStructuralVariantAllele())) {
                System.err.println("" + vt._1().getContig() + ":" + vt._1().getStart() + " " + vt._2().size() + " templates ");
            } else {
                System.err.println("" + vt._1().getContig() + ":" + vt._1().getStart() + " not supported ");
            }
        });
        return variants;
    }

    private static boolean structuralVariantAlleleIsSupported(final StructuralVariantAllele structuralVariantAllele) {
        switch (structuralVariantAllele) {
            case INS:
            case DEL:
                return true;
            default:
                return false;
        }
    }

    private static Pattern FASTQ_READ_NAME_FORMAT = Pattern.compile("^(.*)\\/(\\d+)$");

    private static Tuple2<String, Integer> templateName(final SVFastqUtils.FastqRead read) {
        final String readName = read.getName();
        final Matcher matcher = FASTQ_READ_NAME_FORMAT.matcher(readName);
        final String templateName;
        final int fragmentNumber;
        if (matcher.find()) {
            templateName = matcher.group(1);
            fragmentNumber = Integer.parseInt(matcher.group(2));
        } else {
            templateName = readName;
            fragmentNumber = -1;
        }
        return new Tuple2<>(templateName, fragmentNumber);
    }

    private void tearDown(final JavaSparkContext ctx) {
        variantsSource = null;
    }
}
