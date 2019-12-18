package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

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

public class SVGenotype extends TwoPassVariantWalker {

    private final static String NL = String.format("%n");

    public static final String SV_GENOTYPE_PYTHON_SCRIPT = "svgenotyper_genotype.py";

    static final String USAGE_ONE_LINE_SUMMARY = "Run model to genotype structural variants";
    static final String USAGE_SUMMARY = "Run model to genotype structural variants";

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

    @Argument(fullName = "discrete-samples", doc = "Number of samples per iteration for discrete distribution", optional = true)
    private int discreteSamples = 1000;

    @Argument(fullName = "discrete-log-freq", doc = "Number of iterations between log messages for discrete sampling", optional = true)
    private int discreteLogFreq = 100;

    @Argument(fullName = "jit", doc = "Enable JIT compilation", optional = true)
    private boolean enableJit = false;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);

    private StructuralVariantType svType;
    private VariantContextWriter vcfWriter;
    private SVGenotypeEngineFromModel genotypeEngine;
    private BufferedReader modelOutput;
    private List<String> sampleList;

    public static List<String> FORMAT_FIELDS = Lists.newArrayList(
            GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE,
            GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE,
            GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE,
            GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY,
            GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY
    );

    @Override
    public void onTraversalStart() {}

    @Override
    public void firstPassApply(final VariantContext variant,
                               final ReadsContext readsContext,
                               final ReferenceContext referenceContext,
                               final FeatureContext featureContext) {
        validateRecord(variant);
    }

    private void validateRecord(final VariantContext variant) {
        if (svType == null) {
            svType = variant.getStructuralVariantType();
        } else {
            if (!variant.getStructuralVariantType().equals(svType)) {
                throw new UserException.BadInput("Variants must all have the same SVTYPE. First variant was "
                        + svType.name() + " but found " + variant.getStructuralVariantType().name() + " for record " + variant.getID());
            }
        }
    }

    @Override
    public void afterFirstPass() {

        logger.info("Reading genotype model sample list...");
        final Path sampleListPath = Paths.get(modelDir.toString(), modelName + ".sample_ids.list");
        try (final BufferedReader file = new BufferedReader(IOUtils.makeReaderMaybeGzipped(sampleListPath))) {
            sampleList = file.lines().collect(Collectors.toList());
        } catch (final IOException e) {
            throw new RuntimeException("Error reading from genotype model samples list: "
                    + sampleListPath.toAbsolutePath().toString());
        }

        logger.info("Reading VCF sample list...");
        final List<String> vcfSampleList = getHeaderForVariants().getGenotypeSamples();
        if (!(sampleList.size() == vcfSampleList.size() && new HashSet<>(sampleList).containsAll(vcfSampleList))) {
            throw new UserException.BadInput("VCF and genotype model sample sets are not identical");
        }

        // Execute Python code to initialize training
        final File tempDir = IOUtils.createTempDir(modelName + ".");
        final File tempFile = new File(Paths.get(tempDir.getAbsolutePath(), modelName + ".genotypes.tsv").toString());
        logger.info("Executing training script...");
        final boolean result = pythonExecutor.executeScript(
                new Resource(SV_GENOTYPE_PYTHON_SCRIPT, getClass()),
                null,
                generatePythonArguments(tempFile));
        if (!result) {
            throw new GATKException("Python process returned non-zero exit code");
        }
        logger.info("Script completed with normal exit code.");

        logger.info("Reading output file...");
        genotypeEngine = new SVGenotypeEngineFromModel();
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
        final VCFHeader header = getHeaderForVariants();
        SVGenotypeEngineFromModel.getVcfHeaderMetadata().stream().forEach(line -> header.addMetaDataLine(line));
        vcfWriter.writeHeader(header);
    }

    @Override
    public void secondPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext) {
        try {
            vcfWriter.add(genotypeEngine.genotypeFromModel(variant, modelOutput.readLine(), sampleList));
        } catch (final IOException e) {
            throw new GATKException("Error reading model output file", e);
        }
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
        arguments.add("--svtype=" + svType.name());
        arguments.add("--device=" + device);
        arguments.add("--random_seed=" + randomSeed);
        arguments.add("--genotype_discrete_samples=" + discreteSamples);
        arguments.add("--genotype_discrete_log_freq=" + discreteLogFreq);
        if (enableJit) {
            arguments.add("--jit");
        }
        return arguments;
    }

}
