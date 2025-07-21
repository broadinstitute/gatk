package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonExecutorBase;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

/**
 * Annotate a VCF with scores from a PyTorch-based Convolutional Neural Network (CNN).
 *
 * It contains both a 1D model that uses only the reference sequence and variant annotations,
 * and a 2D model that uses reads in addition to the reference sequence and variant annotations.
 *
 * The scores for each variant record will be placed in an INFO field annotation named CNN_1D
 * (if using the 1D model) or CNN_2D (if using the 2D model). These scores represent the
 * log odds of being a true variant versus being false under the trained convolutional neural
 * network.
 *
 * The provided models were trained on short-read human sequencing data, and will likely not perform
 * well for other kinds of sequencing data, or for non-human data. A companion training tool for
 * NVScoreVariants will be released in the future to support users who need to train their own models.
 *
 * Example command for running with the 1D model:
 *
 * <pre>
 * gatk NVScoreVariants \
 *     -V src/test/resources/large/VQSR/recalibrated_chr20_start.vcf \
 *     -R src/test/resources/large/human_g1k_v37.20.21.fasta \
 *     -O output.vcf
 * </pre>
 *
 * Example command for running with the 2D model:
 *
 * <pre>
 * gatk NVScoreVariants \
 *     -V src/test/resources/large/VQSR/recalibrated_chr20_start.vcf \
 *     -R src/test/resources/large/human_g1k_v37.20.21.fasta \
 *     --tensor-type read_tensor \
 *     -I src/test/resources/large/VQSR/g94982_contig_20_start_bamout.bam \
 *     -O output.vcf
 * </pre>
 *
 * <b><i>The PyTorch Python code that this tool relies upon was contributed by engineers at
 * <a href="https://github.com/NVIDIA-Genomics-Research">NVIDIA Genomics Research</a>.
 * We would like to give particular thanks to Babak Zamirai of NVIDIA, who authored
 * the tool, as well as to Ankit Sethia, Mehrzad Samadi, and George Vacek (also of NVIDIA),
 * without whom this project would not have been possible.</i></b>
 */
@CommandLineProgramProperties(
    summary = "Annotate a VCF with scores from a PyTorch-based Convolutional Neural Network (CNN)",
    oneLineSummary = "Annotate a VCF with scores from a PyTorch-based Convolutional Neural Network (CNN)",
    programGroup = VariantFilteringProgramGroup.class
)
@BetaFeature
public class NVScoreVariants extends CommandLineProgram {

    public static final String NV_SCORE_VARIANTS_PACKAGE = "scorevariants";
    public static final String NV_SCORE_VARIANTS_SCRIPT = "nvscorevariants.py";
    public static final String NV_SCORE_VARIANTS_1D_MODEL_FILENAME = "1d_cnn_mix_train_full_bn.pt";
    public static final String NV_SCORE_VARIANTS_2D_MODEL_FILENAME = "small_2d.pt";
    public static final String NV_SCORE_VARIANTS_1D_MODEL = Resource.LARGE_RUNTIME_RESOURCES_PATH + "/nvscorevariants/" + NV_SCORE_VARIANTS_1D_MODEL_FILENAME;
    public static final String NV_SCORE_VARIANTS_2D_MODEL = Resource.LARGE_RUNTIME_RESOURCES_PATH + "/nvscorevariants/" + NV_SCORE_VARIANTS_2D_MODEL_FILENAME;

    public enum TensorType {
        reference,
        read_tensor
    }

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output VCF file")
    private File outputVCF;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "Input VCF file containing variants to score")
    private File inputVCF;

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file")
    private File reference;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "BAM file containing reads, if using the 2D model", optional = true)
    private File bam;

    @Argument(fullName = "tensor-type", doc = "Name of the tensors to generate: reference for 1D reference tensors and read_tensor for 2D tensors.", optional = true)
    private TensorType tensorType = TensorType.reference;

    @Argument(fullName = "batch-size", doc = "Batch size", optional = true)
    private int batchSize = 64;

    @Argument(fullName = "random-seed", doc = "Seed to initialize the random number generator", optional = true)
    private int randomSeed = 724;

    @Argument(fullName = "tmp-file", doc = "The temporary VCF-like file where variants scores will be written", optional = true)
    private File tmpFile;

    @Argument(fullName = "accelerator", doc = "Type of hardware accelerator to use (auto, cpu, cuda, mps, tpu, etc)", optional = true)
    private String accelerator = "auto";

    @Override
    protected void onStartup() {
        PythonScriptExecutor.checkIfRunningInGatkLiteDocker();
        PythonScriptExecutor.checkPythonEnvironmentForPackage(NV_SCORE_VARIANTS_PACKAGE);
    }

    @Override
    protected Object doWork() {
        final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(PythonExecutorBase.PythonExecutableName.PYTHON3, true);
        final Resource pythonScriptResource = new Resource(NV_SCORE_VARIANTS_SCRIPT, NVScoreVariants.class);
        final File extractedModelDirectory = extractModelFilesToTempDirectory();

        if ( tmpFile == null ) {
            tmpFile = IOUtils.createTempFile("NVScoreVariants_tmp", ".txt");
        }

        final List<String> arguments = new ArrayList<>(Arrays.asList(
            "--output-file", outputVCF.getAbsolutePath(),
            "--vcf-file", inputVCF.getAbsolutePath(),
            "--ref-file", reference.getAbsolutePath(),
            "--tensor-type", tensorType.name(),
            "--batch-size", Integer.toString(batchSize),
            "--seed", Integer.toString(randomSeed),
            "--tmp-file", tmpFile.getAbsolutePath(),
            "--model-directory", extractedModelDirectory.getAbsolutePath()
        ));

        if (accelerator != null) {
            arguments.addAll(List.of("--accelerator",accelerator));
        }

        if ( tensorType == TensorType.reference && bam != null ) {
            throw new UserException.BadInput("--" + StandardArgumentDefinitions.INPUT_LONG_NAME +
                    " should only be specified when running with --tensor-type " + TensorType.read_tensor.name());
        }
        else if ( tensorType == TensorType.read_tensor && bam == null ) {
            throw new UserException.BadInput("Need to specify a BAM file via --" + StandardArgumentDefinitions.INPUT_LONG_NAME +
                    " when running with --tensor-type " + TensorType.read_tensor.name());
        }

        if ( bam != null ) {
            arguments.addAll(Arrays.asList("--input-file", bam.getAbsolutePath()));
        }

        logger.info("Running Python NVScoreVariants module with arguments: " + arguments);
        final ProcessOutput pythonOutput = pythonExecutor.executeScriptAndGetOutput(
                pythonScriptResource,
                null,
                arguments
        );

        if ( pythonOutput.getExitValue() != 0 ) {
            logger.error("Error running NVScoreVariants Python command:\n" + pythonOutput.getStatusSummary(true));
        }
        
        return pythonOutput.getExitValue();
    }

    private File extractModelFilesToTempDirectory() {
        final File extracted1DModel = IOUtils.writeTempResourceFromPath(NV_SCORE_VARIANTS_1D_MODEL, null);
        final File extracted2DModel = IOUtils.writeTempResourceFromPath(NV_SCORE_VARIANTS_2D_MODEL, null);
        final File modelDirectory = IOUtils.createTempDir("NVScoreVariants_models");

        if ( ! extracted1DModel.renameTo(new File(modelDirectory, NV_SCORE_VARIANTS_1D_MODEL_FILENAME)) ) {
            throw new UserException("Error moving " + extracted1DModel.getAbsolutePath() + " to " + modelDirectory.getAbsolutePath());
        }
        if ( ! extracted2DModel.renameTo(new File(modelDirectory, NV_SCORE_VARIANTS_2D_MODEL_FILENAME)) ) {
            throw new UserException("Error moving " + extracted2DModel.getAbsolutePath() + " to " + modelDirectory.getAbsolutePath());
        }

        logger.info("Extracted models to: " + modelDirectory.getAbsolutePath());
        return modelDirectory;
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();
    }
}
