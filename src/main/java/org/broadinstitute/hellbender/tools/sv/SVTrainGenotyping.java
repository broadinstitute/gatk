package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
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

public class SVTrainGenotyping extends VariantWalker {

    private final static String NL = String.format("%n");
    private static final String DATA_VALUE_SEPARATOR = ";";
    private static final String DATA_TYPE_SEPARATOR = "\t";

    public static final String SV_TRAIN_GENOTYPING_PYTHON_SCRIPT = "svgenotyper_train.py";

    static final String USAGE_ONE_LINE_SUMMARY = "Train model to genotype structural variants";
    static final String USAGE_SUMMARY = "Runs training on a set of variants and generates a genotyping model.";

    @Argument(fullName = "coverage-file", doc = "Tab-delimited table of sample mean coverage")
    private File coverageFile;

    @Argument(fullName = "output-name", doc = "Output name")
    private String outputName;

    @Argument(fullName = "output-dir", doc = "Output directory")
    private String outputDir;

    @Argument(fullName = "device", doc = "Device for Torch backend (e.g. \"cpu\", \"cuda\")", optional = true)
    private String device = "cpu";

    @Argument(fullName = "random-seed", doc = "PRNG seed", optional = true)
    private int randomSeed = 92837488;

    @Argument(fullName = "num-states", doc = "Max number of genotype states (0 for auto)", optional = true)
    private int numStates = 0;

    @Argument(fullName = "depth-dilution-factor", doc = "Dilution factor for depth posteriors", optional = true)
    private double depthDilutionFactor = 1e-6;

    @Argument(fullName = "eps-pe", doc = "Mean of PE noise prior", optional = true)
    private double epsilonPE  = 0.01;

    @Argument(fullName = "eps-sr1", doc = "Mean of SR1 noise prior", optional = true)
    private double epsilonSR1  = 0.01;

    @Argument(fullName = "eps-sr2", doc = "Mean of SR2 noise prior", optional = true)
    private double epsilonSR2  = 0.01;

    @Argument(fullName = "phi-pe", doc = "Variance of PE bias prior", optional = true)
    private double phiPE  = 0.1;

    @Argument(fullName = "phi-sr1", doc = "Variance of SR1 bias prior", optional = true)
    private double phiSR1  = 0.1;

    @Argument(fullName = "phi-sr2", doc = "Variance of SR2 bias prior", optional = true)
    private double phiSR2  = 0.1;

    @Argument(fullName = "read-length", doc = "Library read length", optional = true)
    private int readLength  = 150;

    @Argument(fullName = "lr-decay", doc = "Learning rate decay constant (lower is faster)", optional = true)
    private double lrDecay = 1000;

    @Argument(fullName = "lr-min", doc = "Minimum learning rate", optional = true)
    private double lrMin = 1e-3;

    @Argument(fullName = "lr-init", doc = "Initial learning rate", optional = true)
    private double lrInit = 0.01;

    @Argument(fullName = "adam-beta1", doc = "ADAM beta1 constant", optional = true)
    private double adamBeta1 = 0.9;

    @Argument(fullName = "adam-beta2", doc = "ADAM beta2 constant", optional = true)
    private double adamBeta2 = 0.999;

    @Argument(fullName = "max-iter", doc = "Max number of training iterations", optional = true)
    private int maxIter = 2000;

    @Argument(fullName = "iter-log-freq", doc ="Number of iterations between log messages", optional = true)
    private int iterLogFreq = 50;

    @Argument(fullName = "jit", doc = "Enable JIT compilation", optional = true)
    private boolean enableJit = false;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);

    private File samplesFile;
    private StructuralVariantType svType;
    private File variantsFile;
    private PrintStream variantsFileStream;

    private static final int INTRACHROMOSOMAL_LENGTH = Integer.MAX_VALUE;

    // Assumes INV have been converted to BND with SVCluster
    public static List<StructuralVariantType> SV_TYPES = Lists.newArrayList(
            StructuralVariantType.DEL,
            StructuralVariantType.DUP,
            StructuralVariantType.INS,
            StructuralVariantType.BND
    );

    public static List<String> FORMAT_FIELDS = Lists.newArrayList(
            GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE,
            GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE,
            GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE,
            GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY,
            GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY
    );

    @Override
    public void onTraversalStart() {
        samplesFile = createSampleList();
        variantsFile = IOUtils.createTempFile(outputName + ".variants", ".tsv");
        try {
            variantsFileStream = new PrintStream(variantsFile);
        } catch (final IOException e) {
            throw new GATKException("Could not create temporary file: " + variantsFile.getAbsolutePath());
        }
        PythonScriptExecutor.checkPythonEnvironmentForPackage("svgenotyper");
    }

    @Override
    public void apply(final VariantContext variant,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        validateRecord(variant);
        variantsFileStream.println(encodeVariant(variant));
    }

    private void validateRecord(final VariantContext variant) {
        if (isDepthOnly(variant)) {
            throw new UserException.BadInput("Depth-only variant not supported:" + variant.getID());
        }

        //Expect records to have uniform SVTYPE
        final StructuralVariantType variantType = variant.getStructuralVariantType();
        if (svType == null) {
            if (!SV_TYPES.contains(variantType)) {
                final String svTypesString = String.join(",", SV_TYPES.stream().map(StructuralVariantType::name).collect(Collectors.toList()));
                throw new UserException.BadInput("Unsupported variant type in first record:" + variantType.name() + ". Must be one of: " + svTypesString);
            } else {
                svType = variantType;
            }
        } else if (!svType.equals(variantType)) {
            throw new UserException.BadInput("Variants must all have the same SVTYPE. First variant was "
                    + svType.name() + " but found " + variantType.name() + " for record " + variant.getID());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        variantsFileStream.close();
        logger.info("Executing training script...");
        final boolean result = pythonExecutor.executeScript(
                new Resource(SV_TRAIN_GENOTYPING_PYTHON_SCRIPT, getClass()),
                null,
                generatePythonArguments());
        if (!result) {
            throw new GATKException("Python process returned non-zero exit code");
        }
        logger.info("Script completed with normal exit code.");
        return null;
    }

    private String encodeVariant(final VariantContext variant) {
        final Map<String, StringBuilder> stringBuilderMap = FORMAT_FIELDS.stream()
                .collect(Collectors.toMap(s -> s, s -> new StringBuilder()));
        final int numGenotypes = variant.getNSamples();
        int genotypeIndex = 0;
        for (final Genotype genotype : variant.getGenotypes()) {
            final String separator = genotypeIndex == numGenotypes - 1 ? "" : DATA_VALUE_SEPARATOR;
            for (final Map.Entry<String, StringBuilder> entry : stringBuilderMap.entrySet()) {
                entry.getValue().append(genotype.getExtendedAttribute(entry.getKey()) + separator);
            }
            genotypeIndex++;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(variant.getID() + DATA_TYPE_SEPARATOR);
        if (!variant.hasAttribute(GATKSVVCFConstants.SVLEN)) {
            throw new UserException.BadInput("Expected variant attribute " + GATKSVVCFConstants.SVLEN + " not found");
        }

        final int length;
        if (svType.equals(StructuralVariantType.DEL) || svType.equals(StructuralVariantType.DUP)
            || svType.equals(StructuralVariantType.INV)) {
            length = Integer.valueOf(((String)variant.getAttribute(GATKSVVCFConstants.SVLEN)));
        } else if (svType.equals(StructuralVariantType.INS)) {
            length = 0;
        } else if (svType.equals(StructuralVariantType.BND)) {
            final String chr2 = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            if (chr2 == null) {
                throw new UserException.BadInput("Missing BND variant attribute: " + GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
            }
            if (chr2.equals(variant.getContig())) {
                length = variant.getEnd() - variant.getStart();
            } else {
                length = INTRACHROMOSOMAL_LENGTH;
            }
        } else {
            throw new UserException.BadInput("Length calculation not supported for type: " + svType.name());
        }
        stringBuilder.append(length + DATA_TYPE_SEPARATOR);

        final int numFields = FORMAT_FIELDS.size();
        int attributeIndex = 0;
        for (final String attribute : FORMAT_FIELDS) {
            final String separator = attributeIndex == numFields - 1 ? "" : DATA_TYPE_SEPARATOR;
            stringBuilder.append(stringBuilderMap.get(attribute).toString() + separator);
            attributeIndex++;
        }
        return stringBuilder.toString();
    }

    private boolean isDepthOnly(final VariantContext variant) {
        return variant.getAttributeAsString(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, "").equals(GATKSVVCFConstants.DEPTH_ALGORITHM);
    }

    private File createSampleList() {
        final List<String> samples = getHeaderForVariants().getSampleNamesInOrder();
        return IOUtils.writeTempFile(samples, outputName + ".samples", ".tmp");
    }

    private List<String> generatePythonArguments() {
        final List<String> arguments = new ArrayList<>();
        arguments.add("--variants_file=" + variantsFile.getAbsolutePath());
        arguments.add("--coverage_file=" + coverageFile.getAbsolutePath());
        arguments.add("--samples_file=" + samplesFile.getAbsolutePath());
        arguments.add("--output_name=" + outputName);
        arguments.add("--output_dir=" + outputDir);
        arguments.add("--svtype=" + svType.name());
        arguments.add("--device=" + device);

        final int numStatesArg;
        if (svType.equals(StructuralVariantType.DEL) || svType.equals(StructuralVariantType.INS) || svType.equals(StructuralVariantType.BND)) {
            numStatesArg = 3;
        } else if (svType.equals(StructuralVariantType.DUP)) {
            numStatesArg = 5;
        } else {
            throw new UserException.BadInput("Unsupported SV type: " + svType.name());
        }
        arguments.add("--num_states=" + numStatesArg);

        arguments.add("--random_seed=" + randomSeed);
        arguments.add("--depth_dilution_factor=" + depthDilutionFactor);
        arguments.add("--eps_pe=" + epsilonPE);
        arguments.add("--eps_sr1=" + epsilonSR1);
        arguments.add("--eps_sr2=" + epsilonSR2);
        arguments.add("--phi_pe=" + phiPE);
        arguments.add("--phi_sr1=" + phiSR1);
        arguments.add("--phi_sr2=" + phiSR2);
        arguments.add("--read_length=" + readLength);
        arguments.add("--lr_decay=" + lrDecay);
        arguments.add("--lr_min=" + lrMin);
        arguments.add("--lr_init=" + lrInit);
        arguments.add("--adam_beta1=" + adamBeta1);
        arguments.add("--adam_beta2=" + adamBeta2);
        arguments.add("--max_iter=" + maxIter);
        arguments.add("--iter_log_freq=" + iterLogFreq);
        if (enableJit) {
            arguments.add("--jit");
        }
        return arguments;
    }
}
