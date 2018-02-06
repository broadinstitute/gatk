package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by sam on 2/5/18.
 */
@CommandLineProgramProperties(
        summary = "Write tensors",
        oneLineSummary = "Write tensors",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class NeuralNetTrain extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
            doc = "Reference fasta file.", optional = true)
    private String reference = "";


    @Argument(fullName = "write-tensors", shortName = "wt", doc = "Maximum number of truth VCF sites to check.", optional = true)
    private boolean writeTensors = true;

    @Argument(fullName ="input-vcf", shortName = "iv", doc = "Input VCF file", optional = true)
    private String inputVcf = "";

    @Argument(fullName = "data-dir", shortName = "dd", doc = "Directory of training tensors, created if write-tensors is true, otherwise read.", optional = true)
    private String dataDir = "";

    @Argument(fullName = "output-dir", shortName = "od", doc = "Directory where models and plots will be saved.", optional = true)
    private String outputDir = "./";

    @Argument(fullName = "truth-vcf", shortName = "tv", doc = "Validated VCF file.", optional = true)
    private String truthVcf = "";

    @Argument(fullName = "truth-bed", shortName = "tb", doc = "Confident region of the validated VCF file.", optional = true)
    private String truthBed = "";

    @Argument(fullName = "tensor-name", shortName = "tn", doc = "Name of the tensors to generate.")
    private String tensorName = "reference";

    @Argument(fullName = "annotation-set", shortName = "as", doc = "Which set of annotations to use.", optional = true)
    private String annotationSet = "gatk";

    @Argument(fullName = "bam-file", shortName = "bf", doc = "BAM or BAMout file to use for read data when generating 2D tensors.", optional = true)
    private String bamFile = "";

    @Argument(fullName = "identifier", shortName = "id", doc = "Name of the model to be trained.", optional = true)
    private String id = "no_name";

    @Argument(fullName = "samples", shortName = "s", doc = "Maximum number of truth VCF sites to check.", optional = true)
    private int samples = 1200;

    @Argument(fullName = "epochs", shortName = "e", doc = "Maximum number of training epochs.", optional = true)
    private int epochs = 10;

    // Start the Python executor. This does not actually start the Python process, but fails if python can't be located
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);


    @Override
    protected void onStartup() {
        PythonScriptExecutor.checkPythonEnvironmentForPackage("vqsr_cnn");

    }

    @Override
    protected Object doWork() {
        final Resource pythonScriptResource = new Resource("training.py", org.broadinstitute.hellbender.tools.walkers.vqsr.NeuralNetTranches.class);
        List<String> arguments = new ArrayList<>(Arrays.asList(
                "--reference_fasta", reference,
                "--input_vcf", inputVcf,
                "--bam_file", bamFile,
                "--train_vcf", truthVcf,
                "--bed_file", truthBed,
                "--tensor_name", tensorName,
                "--annotation_set", annotationSet,
                "--samples", Integer.toString(samples),
                "--epochs", Integer.toString(epochs),
                "--data_dir", dataDir,
                "--output_dir", outputDir,
                "--id", id));

        if (writeTensors && tensorName.equals("reference")) {
            arguments.addAll(Arrays.asList("--mode", "write_reference_and_annotation_tensors"));
        } else if (writeTensors && tensorName.equals("read_tensor")) {
            arguments.addAll(Arrays.asList("--mode", "write_read_and_annotation_tensors"));
        } else if (tensorName.equals("reference")) {
            arguments.addAll(Arrays.asList("--mode", "train_on_reference_tensors_and_annotations"));
        } else if (tensorName.equals("read_tensor")) {
            arguments.addAll(Arrays.asList("--mode", "train_on_read_tensors_and_annotations"));
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
