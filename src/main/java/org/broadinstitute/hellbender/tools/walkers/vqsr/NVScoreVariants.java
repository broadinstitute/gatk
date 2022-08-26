package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
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

/**
 * Annotate a VCF with scores from a PyTorch-based Convolutional Neural Network (CNN).
 *
 * It contains both a 1D model that uses only the reference sequence and variant annotations,
 * and a 2D model that uses reads in addition to the reference sequence and variant annotations.
 *
 * Running this tool currently requires activating a custom conda environment, stored in
 * scripts/nvscorevariants_environment.yml. To do this, run these commands:
 *
 * <pre>
 *     conda env create -f scripts/nvscorevariants_environment.yml
 *     conda activate scorevariants
 *     pip install build/gatkPythonPackageArchive.zip  (created via ./gradlew pythonPackageArchive)
 * </pre>
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
 */
@CommandLineProgramProperties(
    summary = "Annotate a VCF with scores from a PyTorch-based Convolutional Neural Network (CNN)",
    oneLineSummary = "Annotate a VCF with scores from a PyTorch-based Convolutional Neural Network (CNN)",
    programGroup = VariantFilteringProgramGroup.class
)
@ExperimentalFeature
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
    private int batchSize = 32;

    @Argument(fullName = "random-seed", doc = "Seed to initialize the random number generator")
    private int randomSeed = 724;

    @Argument(fullName = "tmp-file", doc = "The temporary VCF-like file where variants scores will be written", optional = true)
    private File tmpFile;

    // TODO: this argument does not appear to be hooked up in the underlying Python script,
    // TODO: so we won't expose it in this tool until it is:
    //
    // @Argument(fullName = "num-gpus", doc = "Number of GPUs", optional = true)
    // private int numGPUs;
    
    @Override
    protected void onStartup() {
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
