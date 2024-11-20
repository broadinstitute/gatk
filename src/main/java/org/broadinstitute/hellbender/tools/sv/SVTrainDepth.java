package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVIntervalVariantDecoder;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
        summary = SVTrainDepth.USAGE_SUMMARY,
        oneLineSummary = SVTrainDepth.USAGE_ONE_LINE_SUMMARY,
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)

public class SVTrainDepth extends FeatureWalker<DepthEvidence> {

    private static final String DATA_VALUE_SEPARATOR = ";";
    private static final String DATA_TYPE_SEPARATOR = "\t";

    public static final String CNV_TRAIN_GENOTYPING_PYTHON_SCRIPT = "svgenotyper_cnv_train.py";

    static final String USAGE_ONE_LINE_SUMMARY = "Train model to genotype structural variants";
    static final String USAGE_SUMMARY = "Runs training on a set of variants and generates a genotyping model.";

    @Argument(fullName = "coverage-file", doc = "Tab-delimited table of per-base sample depth of the reference state")
    private GATKPath sampleDepthFilePath;

    @Argument(fullName = "depth-file", doc = "Read depth evidence file")
    private GATKPath depthEvidenceFilePath;

    @Argument(fullName = SVCopyNumberPosteriors.CONTIG_PLOIDY_CALLS_LONG_NAME, doc = "Contig ploidy calls file. Can be specified for multiple samples.")
    private List<GATKPath> contigPloidyCallFilePaths;

    @Argument(fullName = "output-name", doc = "Output name")
    private String outputName;

    @Argument(fullName = "output-dir", doc = "Output directory")
    private GATKPath outputPath;

    @Argument(fullName = "device", doc = "Device for Torch backend (e.g. \"cpu\", \"cuda\")", optional = true)
    private String device = "cpu";

    @Argument(fullName = "random-seed", doc = "PRNG seed", optional = true)
    private int randomSeed = 92837488;

    @Argument(fullName = "num-states", doc = "Max number of genotype states (must be 5)", optional = true)
    private int numStates = 5;

    @Argument(fullName = "mu-eps", doc = "Mean of PE noise prior", optional = true)
    private double muEps  = 0.01;

    @Argument(fullName = "var-phi-sample", doc = "Variance of per-sample count bias prior", optional = true)
    private double varPhiSample  = 0.1;

    @Argument(fullName = "var-phi-bin", doc = "Variance of per-bin count bias prior", optional = true)
    private double varPhiBin  = 0.1;

    @Argument(fullName = "alpha-ref", doc = "Hardy-Weinberg Dirichlet concentration for reference state", optional = true)
    private double alphaRef  = 18.;

    @Argument(fullName = "alpha-non-ref", doc = "Hardy-Weinberg Dirichlet concentration for non-reference states", optional = true)
    private double alphaNonRef  = 1.;

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
    private File modelDataFile;
    private PrintStream modelDataFileStream;
    private List<Integer> sampleIndexes;
    private List<String> samplesList;
    private Map<String, Map<String,Integer>> sampleContigPloidyMap;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(BafEvidence.class) || featureType.equals(DepthEvidence.class)
                || featureType.equals(DiscordantPairEvidence.class) || featureType.equals(SplitReadEvidence.class);
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return depthEvidenceFilePath;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        try {
            samplesList = new BufferedReader(new FileReader(sampleDepthFilePath.toPath().toFile())).lines()
                    .map(s -> s.split("\t")[0])
                    .collect(Collectors.toList());
        } catch (final IOException e) {
            throw new GATKException("Error reading sample depth file", e);
        }
        samplesFile = createSampleList();
        final List<String> headerSamples = getSamplesFromHeader();
        final Map<String, Integer> headerSampleIndexes = IntStream.range(0, headerSamples.size())
                .mapToObj(i -> new AbstractMap.SimpleImmutableEntry<>(headerSamples.get(i), i))
                .collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue()));
        sampleIndexes = samplesList.stream().map(headerSampleIndexes::get).collect(Collectors.toList());

        modelDataFile = IOUtils.createTempFile(outputName + ".model_data", ".tsv");
        try {
            modelDataFileStream = new PrintStream(modelDataFile);
        } catch (final IOException e) {
            throw new GATKException("Could not create temporary file: " + modelDataFile.getAbsolutePath());
        }
        PythonScriptExecutor.checkPythonEnvironmentForPackage("svgenotyper");

        final Collection<CalledContigPloidyCollection> contigPloidyCollections = contigPloidyCallFilePaths.stream()
                .map(p -> new CalledContigPloidyCollection(p.toPath().toFile()))
                .collect(Collectors.toList());
        final GermlineCNVIntervalVariantDecoder cnvDecoder = new GermlineCNVIntervalVariantDecoder(contigPloidyCollections);
        sampleContigPloidyMap = cnvDecoder.getSampleContigPloidyMap();
    }

    private List<String> getSamplesFromHeader() {
        final Object header = getDrivingFeaturesHeader();
        if (!(header instanceof DepthEvidenceCodec.DepthEvidenceMetadata)) {
            throw new GATKException("Unexpected header type");
        }
        return ((DepthEvidenceCodec.DepthEvidenceMetadata)header).getSamples();
    }

    @Override
    public void apply(final DepthEvidence feature,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext ) {
        modelDataFileStream.println(encode(feature));
    }

    protected final String encode(final DepthEvidence feature) {
        final StringBuilder sb = new StringBuilder(6);
        final String contig = feature.getContig();
        sb.append(contig);
        sb.append(DATA_TYPE_SEPARATOR);
        sb.append(feature.getStart());
        sb.append(DATA_TYPE_SEPARATOR);
        sb.append(feature.getEnd() - feature.getStart() + 1);
        sb.append(DATA_TYPE_SEPARATOR);
        final int[] allCounts = feature.getCounts();
        final List<String> counts = sampleIndexes.stream().map(i -> String.valueOf(allCounts[i])).collect(Collectors.toList());
        sb.append(String.join(DATA_VALUE_SEPARATOR, counts));
        sb.append(DATA_TYPE_SEPARATOR);
        final List<String> ploidies = samplesList.stream().map(s -> sampleContigPloidyMap.get(s).get(contig).toString()).collect(Collectors.toList());
        sb.append(String.join(DATA_VALUE_SEPARATOR, ploidies));
        return sb.toString();
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        modelDataFileStream.close();

        logger.info("Executing training script...");
        final boolean result = pythonExecutor.executeScript(
                new Resource(CNV_TRAIN_GENOTYPING_PYTHON_SCRIPT, getClass()),
                null,
                generatePythonArguments());
        if (!result) {
            throw new GATKException("Python process returned non-zero exit code");
        }
        logger.info("Script completed with normal exit code.");

        return null;
    }

    private File createSampleList() {
        return IOUtils.writeTempFile(samplesList, outputName + ".samples", ".tmp");
    }

    private List<String> generatePythonArguments() {
        final List<String> arguments = new ArrayList<>();
        arguments.add("--data_file=" + modelDataFile.getAbsolutePath());
        arguments.add("--sample_depth_file=" + sampleDepthFilePath.toString());
        arguments.add("--samples_file=" + samplesFile.getAbsolutePath());
        arguments.add("--output_name=" + outputName);
        arguments.add("--output_dir=" + outputPath);
        arguments.add("--device=" + device);
        arguments.add("--num_states=" + numStates);

        arguments.add("--random_seed=" + randomSeed);
        arguments.add("--mu_eps=" + muEps);
        arguments.add("--var_phi_sample=" + varPhiSample);
        arguments.add("--var_phi_bin=" + varPhiBin);
        arguments.add("--alpha_ref=" + alphaRef);
        arguments.add("--alpha_non_ref=" + alphaNonRef);
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
