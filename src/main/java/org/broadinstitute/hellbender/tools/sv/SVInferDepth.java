package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Annotate a VCF with scores from a Convolutional Neural Network (CNN).
 *
 * This tool streams variants and their reference context to a python program,
 * which evaluates a pre-trained neural network on each variant.
 * The default models were trained on single-sample VCFs.
 * The default model should not be used on VCFs with annotations from joint call-sets.
 *
 * <h3>1D Model with pre-trained architecture</h3>
 *
 * <pre>
 * gatk CNNScoreVariants \
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf
 * </pre>
 *
 */
@ExperimentalFeature
@DocumentedFeature
@CommandLineProgramProperties(
        summary = SVTrainGenotyping.USAGE_SUMMARY,
        oneLineSummary = SVTrainGenotyping.USAGE_ONE_LINE_SUMMARY,
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)

public class SVInferDepth extends GATKTool {

    private static final String COLUMN_SEPARATOR = "\t";
    private static final String FIRST_DIM_SEPARATOR = ";";
    private static final String SECOND_DIM_SEPARATOR = ",";

    public static final String SV_GENOTYPE_PYTHON_SCRIPT = "svgenotyper_cnv_infer.py";

    static final String USAGE_ONE_LINE_SUMMARY = "Run model to infer CNVs";
    static final String USAGE_SUMMARY = "Run model to infer CNVs";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF")
    private GATKPath outputVcf;

    @Argument(fullName = "model-name", doc = "Model name")
    private String modelName;

    @Argument(fullName = "model-dir", doc = "Model directory")
    private GATKPath modelDir;

    @Argument(fullName = "device", doc = "Device for Torch backend (e.g. \"cpu\", \"cuda\")", optional = true)
    private String device = "cpu";

    @Argument(fullName = "random-seed", doc = "PRNG seed", optional = true)
    private int randomSeed = 92837488;

    @Argument(fullName = "predictive-samples", doc = "Number of samples per iteration for predictive distribution", optional = true)
    private int predictiveSamples = 100;

    @Argument(fullName = "predictive-iter", doc = "Number of iterations for predictive distribution", optional = true)
    private int predictiveIter = 10;

    @Argument(fullName = "discrete-samples", doc = "Number of samples per iteration for discrete distribution", optional = true)
    private int discreteSamples = 1000;

    @Argument(fullName = "discrete-log-freq", doc = "Number of iterations between log messages for discrete sampling", optional = true)
    private int discreteLogFreq = 100;

    @Argument(fullName = "jit", doc = "Enable JIT compilation", optional = true)
    private boolean enableJit = false;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);

    private VariantContextWriter vcfWriter;
    private BufferedReader modelOutput;
    private List<String> sampleList;

    @Override
    public void onStartup() {
        super.onStartup();

        logger.info("Reading genotype model sample list...");
        final Path sampleListPath = Paths.get(modelDir.toString(), modelName + ".sample_ids.list");
        try (final BufferedReader file = new BufferedReader(IOUtils.makeReaderMaybeGzipped(sampleListPath))) {
            sampleList = file.lines().collect(Collectors.toList());
        } catch (final IOException e) {
            throw new RuntimeException("Error reading from genotype model samples list: "
                    + sampleListPath.toAbsolutePath().toString());
        }

        // Execute Python code to initialize training
        final File tempDir = IOUtils.createTempDir(modelName + ".");
        final File tempFile = new File(Paths.get(tempDir.getAbsolutePath(), modelName + ".genotypes.tsv").toString());
        logger.info("Executing inference script...");
        final boolean result = pythonExecutor.executeScript(
                new Resource(SV_GENOTYPE_PYTHON_SCRIPT, getClass()),
                null,
                generatePythonArguments(tempFile));
        if (!result) {
            throw new GATKException("Python process returned non-zero exit code");
        }
        logger.info("Script completed with normal exit code.");

        logger.info("Reading output file...");
        try {
            modelOutput = new BufferedReader(IOUtils.makeReaderMaybeGzipped(tempFile.toPath()));
            final String header = modelOutput.readLine();
            if (!header.startsWith("#")) {
                throw new RuntimeException("Expected Python output file header starting with #");
            }
        } catch (final IOException e) {
            throw new RuntimeException("Error reading from Python output file in: " + tempFile.getAbsolutePath());
        }
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(composeVariantContextHeader());
    }

    public VCFHeader composeVariantContextHeader() {
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException.BadInput("Sequence dictionary file required.");
        }
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        final Set<VCFHeaderLine> vcfDefaultToolHeaderLines = getDefaultToolVCFHeaderLines();
        final VCFHeader result = new VCFHeader(Collections.emptySet(), sampleList);

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        result.setSequenceDictionary(sequenceDictionary);

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(result::addMetaDataLine);

        /* header lines related to genotype formatting */
        result.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        result.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FIELD, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Integer, "Max likelihood copy number"));
        result.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Integer, "Copy number log posterior (in Phred-scale) rounded down"));
        result.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DEPTH_P_HARDY_WEINBERG_LOSS_FIELD, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Float, "Hardy-Weinberg probability for CNV loss"));
        result.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DEPTH_P_HARDY_WEINBERG_GAIN_FIELD, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Float, "Hardy-Weinberg probability for CNV gain"));
        result.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DEPTH_BACKGROUND_FIELD, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Float, "Depth background signal fraction"));
        result.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DEPTH_MEAN_BIAS_FIELD, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Float, "Depth mean bias"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));
        return result;
    }

    private VariantContext createVariantContextFromOutputLine(final String line) {
        final String[] tokens = line.trim().split(COLUMN_SEPARATOR);
        Utils.validate(tokens.length == 8, "Expected 8 columns but found " + tokens.length + " in line: " + line);
        final String contig = tokens[0];
        final int start = Integer.parseInt(tokens[1]);
        final int end = start + Integer.parseInt(tokens[2]);
        final String[] stateProbStringArray = tokens[3].split(FIRST_DIM_SEPARATOR);
        final double eps = Double.parseDouble(tokens[4]);
        final double phiBin = Double.parseDouble(tokens[5]);
        final double pHWLoss = Double.parseDouble(tokens[6]);
        final double pHWGain = Double.parseDouble(tokens[7]);
        final String id = String.join("_", Arrays.asList(contig, String.valueOf(start), String.valueOf(end)));

        final VariantContextBuilder builder = new VariantContextBuilder("", contig, start, end, Arrays.asList(Allele.REF_N, Allele.SV_SIMPLE_CNV));
        builder.id(id);
        builder.attribute(GATKSVVCFConstants.DEPTH_P_HARDY_WEINBERG_LOSS_FIELD, pHWLoss);
        builder.attribute(GATKSVVCFConstants.DEPTH_P_HARDY_WEINBERG_GAIN_FIELD, pHWGain);
        builder.attribute(GATKSVVCFConstants.DEPTH_BACKGROUND_FIELD, eps);
        builder.attribute(GATKSVVCFConstants.DEPTH_MEAN_BIAS_FIELD, phiBin);

        if (stateProbStringArray.length != sampleList.size()) {
            throw new UserException.BadInput("Encountered line with " + stateProbStringArray.length + " sample state posteriors but the sample list is of length " + sampleList.size());
        }
        final List<Genotype> genotypes = new ArrayList<>(stateProbStringArray.length);
        for (int i = 0; i < stateProbStringArray.length; i++) {
            final double[] sampleStateProbs = Arrays.stream(stateProbStringArray[i].split(SECOND_DIM_SEPARATOR)).mapToDouble(Double::valueOf).toArray();
            final int[] sampleStatePhreds = DoubleStream.of(sampleStateProbs).mapToInt(p -> Math.round(QualityUtil.getPhredScoreFromErrorProbability(Math.max(Double.MIN_VALUE, p)))).toArray();
            final int copyState = MathUtils.maxElementIndex(sampleStateProbs);
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleList.get(i));
            genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, sampleStatePhreds);
            genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_FIELD, copyState);
            genotypeBuilder.alleles(Collections.singletonList(Allele.NO_CALL));
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);

        return builder.make();
    }

    @Override
    public void traverse() {
        modelOutput.lines().forEachOrdered(line -> vcfWriter.add(createVariantContextFromOutputLine(line)));
    }

    @Override
    public Object onTraversalSuccess() {
        vcfWriter.close();
        try {
            modelOutput.close();
        } catch (final IOException e) {
            throw new GATKException("Error closing model output file", e);
        }
        return null;
    }

    private List<String> generatePythonArguments(final File output) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("--output=" + output.getAbsolutePath());
        arguments.add("--model_name=" + modelName);
        arguments.add("--model_dir=" + modelDir);
        arguments.add("--device=" + device);
        arguments.add("--random_seed=" + randomSeed);
        arguments.add("--infer_predictive_samples=" + predictiveSamples);
        arguments.add("--infer_predictive_iter=" + predictiveIter);
        arguments.add("--infer_discrete_samples=" + discreteSamples);
        arguments.add("--infer_discrete_log_freq=" + discreteLogFreq);
        if (enableJit) {
            arguments.add("--jit");
        }
        return arguments;
    }

}
