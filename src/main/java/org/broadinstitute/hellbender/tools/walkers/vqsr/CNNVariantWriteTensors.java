package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Write variant tensors for training a Convolutional Neural Network (CNN) for filtering variants.
 * After running this tool, a model can be trained with the {@link CNNVariantTrain} tool.
 *
 *
 * <h3>Inputs</h3>
 * <ul>
 *      <li>The input variants to make into tensors.
 *      These variant calls must be annotated with the standard best practices annotations.</li>
 *      <li>The truth VCF has validated variant calls, like those in the genomes in a bottle,
 *      platinum genomes, or CHM VCFs.  Variants in both the input VCF and the truth VCF
 *      will be used as positive training data.</li>
 *      <li>The truth BED is a bed file define the confident region for the validated calls.
 *      Variants from the input VCF inside this region, but not included in the truth VCF
 *      will be used as negative training data.</li>
 *      <li>The --tensor-type argument determines what types of tensors will be written.
 *      Set it to "reference" to write 1D tensors or "read_tensor" to write 2D tensors.</li>
 *      <li>The bam-file argument is necessary to write 2D tensors which incorporate read data.</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 * <li>data-dir This directory is created and populated with variant tensors.
 *  it will be divided into training, validation and test sets and each set will be further divided into
 *  positive and negative SNPs and INDELs.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Write Reference Tensors</h4>
 * <pre>
 * gatk CNNVariantWriteTensors \
 *   -R reference.fasta \
 *   -V input.vcf.gz \
 *   --truth-vcf platinum-genomes.vcf \
 *   --truth-bed platinum-confident-region.bed \
 *   --tensor-type reference \
 *   --output-tensor-dir my-tensor-folder
 * </pre>
 *
 * <h4>Write Read Tensors</h4>
 * <pre>
 * gatk CNNVariantWriteTensors \
 *   -R reference.fasta \
 *   -V input.vcf.gz \
 *   --truth-vcf platinum-genomes.vcf \
 *   --truth-bed platinum-confident-region.bed \
 *   --tensor-type read_tensor \
 *   --bam-file input.bam \
 *   --output-tensor-dir my-tensor-folder
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Write variant tensors for training a CNN to filter variants",
        oneLineSummary = "Write variant tensors for training a CNN to filter variants",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class CNNVariantWriteTensors extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
            doc = "Reference fasta file.")
    private String reference;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            doc = "Input VCF file")
    private String inputVcf;

    @Argument(fullName = "output-tensor-dir", shortName = "output-tensor-dir", doc = "Directory of training tensors. Subdivided into train, valid and test sets.")
    private String outputTensorsDir;

    @Argument(fullName = "truth-vcf", shortName = "truth-vcf", doc = "Validated VCF file.")
    private String truthVcf;

    @Argument(fullName = "truth-bed", shortName = "truth-bed", doc = "Confident region of the validated VCF file.")
    private String truthBed;

    @Argument(fullName = "bam-file", shortName = "bam-file", doc = "BAM or BAMout file to use for read data when generating 2D tensors.", optional = true)
    private String bamFile = "";

    @Argument(fullName = "tensor-type", shortName = "tensor-type", doc = "Name of the tensors to generate.")
    private TensorType tensorType = TensorType.reference;

    @Advanced
    @Argument(fullName = "channels-last", shortName = "channels-last", doc = "Store the channels in the last axis of tensors, tensorflow->true, theano->false", optional = true)
    private boolean channelsLast = true;

    @Advanced
    @Argument(fullName = "annotation-set", shortName = "annotation-set", doc = "Which set of annotations to use.", optional = true)
    private String annotationSet = "best_practices";

    @Argument(fullName = "max-tensors", shortName = "max-tensors", doc = "Maximum number of tensors to write.", optional = true, minValue = 0)
    private int maxTensors = 1000000;

    // Start the Python executor. This does not actually start the Python process, but fails if python can't be located
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);

    @Override
    protected void onStartup() {
        PythonScriptExecutor.checkPythonEnvironmentForPackage("vqsr_cnn");
    }

    @Override
    protected Object doWork() {
        final Resource pythonScriptResource = new Resource("training.py", CNNVariantWriteTensors.class);
        List<String> arguments = new ArrayList<>(Arrays.asList(
                "--reference_fasta", reference,
                "--input_vcf", inputVcf,
                "--bam_file", bamFile,
                "--train_vcf", truthVcf,
                "--bed_file", truthBed,
                "--tensor_name", tensorType.name(),
                "--annotation_set", annotationSet,
                "--samples", Integer.toString(maxTensors),
                "--data_dir", outputTensorsDir));

        if(channelsLast){
            arguments.add("--channels_last");
        } else{
            arguments.add("--channels_first");
        }

        if (tensorType == TensorType.reference) {
            arguments.addAll(Arrays.asList("--mode", "write_reference_and_annotation_tensors"));
        } else if (tensorType == TensorType.read_tensor) {
            arguments.addAll(Arrays.asList("--mode", "write_read_and_annotation_tensors"));
        } else {
            throw new GATKException("Unknown tensor mapping mode:"+ tensorType.name());
        }

        logger.info("Args are:"+ Arrays.toString(arguments.toArray()));
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                arguments
        );
        return pythonReturnCode;
    }

}
