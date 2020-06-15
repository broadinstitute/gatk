package org.broadinstitute.hellbender.tools.walkers.vqsr;

import java.util.*;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.io.IOException;
import java.util.stream.Collectors;
import java.io.UnsupportedEncodingException;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsamplingIterator;
import org.broadinstitute.hellbender.utils.downsampling.ReservoirDownsampler;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriter;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;

import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import com.intel.gkl.IntelGKLUtils;

/**
 * Annotate a VCF with scores from a Convolutional Neural Network (CNN).
 *
 * This tool streams variants and their reference context to a python program,
 * which evaluates a pre-trained neural network on each variant.
 * The default models were trained on single-sample VCFs.
 * The default model should not be used on VCFs with annotations from joint call-sets.
 *
 * The neural network performs convolutions over the reference sequence surrounding the variant
 * and combines those features with a multilayer perceptron on the variant annotations.
 *
 * 2D models convolve over aligned reads as well as the reference sequence, and variant annotations.
 * 2D models require a SAM/BAM file as input and for the --tensor-type argument to be set
 * to a tensor type which requires reads, as in the example below.
 *
 * Pre-trained 1D and 2D models are included in the distribution.
 * It is possible to train your own models with the tools:
 * {@link CNNVariantWriteTensors} and {@link CNNVariantTrain}.
 * CNNVariantTrain will create a json architecture file and an hd5 weights file, which you can use with this tool.
 *
 * The advanced argument `info-annotation-keys` is available for models trained with different sets info field annotations.
 * In order to do this you must first train your own model with the tools {@link CNNVariantWriteTensors} and {@link CNNVariantTrain}.
 * Otherwise, providing this argument with anything but the standard set of annotations will result in an error.
 *
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
 * <h3>2D Model with pre-trained architecture</h3>
 *
 * <pre>
 * gatk CNNScoreVariants \
 *   -I aligned_reads.bam \
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   -tensor-type read-tensor
 * </pre>
 *
 * <h3>1D Model with user-supplied architecture and weights:</h3>
 *
 * <pre>
 * gatk CNNScoreVariants \
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   -architecture path/to/my_model_folder/1dmodel.json
 *   -weights path/to/my_model_folder/1dmodel.hd5
 * </pre>
 *
 * <h3>2D Model with user-supplied model architecture and weights:</h3>
 *
 * <pre>
 * gatk CNNScoreVariants \
 *   -I aligned_reads.bam \
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   -tensor-type read-tensor \
 *   -architecture path/to/my_model_folder/2dmodel.json
 *   -weights path/to/my_model_folder/2dmodel.hd5
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = CNNScoreVariants.USAGE_SUMMARY,
        oneLineSummary = CNNScoreVariants.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantFilteringProgramGroup.class
)

public class CNNScoreVariants extends TwoPassVariantWalker {
    private final static String NL = String.format("%n");
    static final String USAGE_ONE_LINE_SUMMARY = "Apply a Convolutional Neural Net to filter annotated variants";
    static final String USAGE_SUMMARY = "Annotate a VCF with scores from a Convolutional Neural Network (CNN)." +
            "The CNN determines a Log Odds Score for each variant." +
            "Pre-trained models (1D or 2D) are specified via the architecture argument." +
            "1D models will look at the reference sequence and variant annotations." +
            "2D models look at aligned reads, reference sequence, and variant annotations." +
            "2D models require a BAM file as input as well as the tensor-type argument to be set.";
    static final String DISABLE_AVX_CHECK_NAME = "disable-avx-check";
    static final String AVXREQUIRED_ERROR = "This tool requires AVX instruction set support by default due to its dependency on recent versions of the TensorFlow library.\n" +
            " If you have an older (pre-1.6) version of TensorFlow installed that does not require AVX you may attempt to re-run the tool with the %s argument to bypass this check.\n" +
            " Note that such configurations are not officially supported.";

    private static final int CONTIG_INDEX = 0;
    private static final int POS_INDEX = 1;
    private static final int REF_INDEX = 2;
    private static final int ALT_INDEX = 3;
    private static final int KEY_INDEX = 4;
    private static final int FIFO_STRING_INITIAL_CAPACITY = 1024;
    private static final int MAX_BATCH_SIZE_1D = 1024;
    private static final int MAX_BATCH_SIZE_2D = 64;

    // These constants correspond to constants in the python code set in defines.py. They must be kept in sync.
    private static final String DATA_VALUE_SEPARATOR = ","; // If changed make change in defines.py
    private static final String DATA_TYPE_SEPARATOR = "\t"; // If changed make change in defines.py
    private static final String ANNOTATION_SEPARATOR = ";"; // If changed make change in defines.py
    private static final String ANNOTATION_SET_STRING = "=";// If changed make change in defines.py

    private List<String> defaultAnnotationKeys = new ArrayList<>(Arrays.asList("MQ", "DP", "SOR", "FS", "QD", "MQRankSum", "ReadPosRankSum"));

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private GATKPath outputFile;

    @Argument(fullName = "architecture", shortName = "architecture", doc = "Neural Net architecture configuration json file", optional = true)
    private String architecture;

    @Argument(fullName = "weights", shortName = "weights", doc = "Keras model HD5 file with neural net weights.", optional = true)
    private String weights;

    @Argument(fullName = "tensor-type", shortName = "tensor-type", doc = "Name of the tensors to generate, reference for 1D reference tensors and read_tensor for 2D tensors.", optional = true)
    private TensorType tensorType = TensorType.reference;

    @Argument(fullName = "window-size", shortName = "window-size", doc = "Neural Net input window size", minValue = 0, optional = true)
    private int windowSize = 128;

    @Argument(fullName = "read-limit", shortName = "read-limit", doc = "Maximum number of reads to encode in a tensor, for 2D models only.", minValue = 0, optional = true)
    private int readLimit = 128;

    @Argument(fullName = "filter-symbolic-and-sv", shortName = "filter-symbolic-and-sv", doc = "If set will filter symbolic and and structural variants from the input VCF", optional = true)
    private boolean filterSymbolicAndSV = false;

    @Advanced
    @Argument(fullName="info-annotation-keys", shortName="info-annotation-keys", doc="The VCF info fields to send to python.  This should only be changed if a new model has been trained which expects the annotations provided here.", optional=true)
    private List<String> annotationKeys = defaultAnnotationKeys;

    @Advanced
    @Argument(fullName = "inference-batch-size", shortName = "inference-batch-size", doc = "Size of batches for python to do inference on.", minValue = 1, maxValue = 4096, optional = true)
    private int inferenceBatchSize = 256;

    @Advanced
    @Argument(fullName = "transfer-batch-size", shortName = "transfer-batch-size", doc = "Size of data to queue for python streaming.", minValue = 1, maxValue = 8192, optional = true)
    private int transferBatchSize = 512;

    @Advanced
    @Argument(fullName = "inter-op-threads", shortName = "inter-op-threads", doc = "Number of inter-op parallelism threads to use for Tensorflow", minValue = 0, maxValue = 4096, optional = true)
    private int interOpThreads = 0;

    @Advanced
    @Argument(fullName = "intra-op-threads", shortName = "intra-op-threads", doc = "Number of intra-op parallelism threads to use for Tensorflow", minValue = 0, maxValue = 4096, optional = true)
    private int intraOpThreads = 0;

    @Advanced
    @Argument(fullName = "output-tensor-dir", shortName = "output-tensor-dir", doc = "Optional directory where tensors can be saved for debugging or visualization.", optional = true)
    private String outputTensorsDir = "";

    @Advanced
    @Argument(fullName = DISABLE_AVX_CHECK_NAME, shortName = DISABLE_AVX_CHECK_NAME, doc = "If set, no check will be made for AVX support.  " +
            "Use only if you have installed a pre-1.6 TensorFlow build. ", optional = true)
    private boolean disableAVXCheck = false;

    @Hidden
    @Argument(fullName = "enable-journal", shortName = "enable-journal", doc = "Enable streaming process journal.", optional = true)
    private boolean enableJournal = false;

    @Hidden
    @Argument(fullName = "keep-temp-file", shortName = "keep-temp-file", doc = "Keep the temporary file that python writes scores to.", optional = true)
    private boolean keepTempFile = false;

    @Hidden
    @Argument(fullName = "python-profile", shortName = "python-profile", doc = "Run the tool with the Python CProfiler on and write results to this file.", optional = true)
    private File pythonProfileResults;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final StreamingPythonScriptExecutor<String> pythonExecutor = new StreamingPythonScriptExecutor<>(true);

    private List<String> batchList = new ArrayList<>(inferenceBatchSize);

    private int curBatchSize = 0;
    private int windowEnd = windowSize / 2;
    private int windowStart = windowSize / 2;
    private boolean waitforBatchCompletion = false;

    private File scoreFile; // use java.nio.File here because python code needs to write to this
    private String scoreKey;
    private Scanner scoreScan;
    private VariantContextWriter vcfWriter;
    private String annotationSetString;

    private static String resourcePathReadTensor = Resource.LARGE_RUNTIME_RESOURCES_PATH + "/cnn_score_variants/small_2d.json";
    private static String resourcePathReferenceTensor = Resource.LARGE_RUNTIME_RESOURCES_PATH + "/cnn_score_variants/1d_cnn_mix_train_full_bn.json";

    @Override
    protected String[] customCommandLineValidation() {
        if (tensorType.equals(TensorType.read_tensor)){
            transferBatchSize = Math.max(transferBatchSize, MAX_BATCH_SIZE_2D);
            inferenceBatchSize = Math.max(inferenceBatchSize, MAX_BATCH_SIZE_2D);
        } else if (tensorType.equals(TensorType.reference)){
            transferBatchSize = Math.max(transferBatchSize, MAX_BATCH_SIZE_1D);
            inferenceBatchSize = Math.max(inferenceBatchSize, MAX_BATCH_SIZE_1D);
        }

        if (inferenceBatchSize > transferBatchSize) {
            return new String[]{"Inference batch size must be less than or equal to transfer batch size."};
        }

        if (architecture == null || weights == null){
            if (!tensorType.equals(TensorType.read_tensor) && !tensorType.equals(TensorType.reference)){
                return new String[]{"No default architecture for tensor type:" + tensorType.name()};
            }
        }
        return null;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected CountingVariantFilter makeVariantFilter() {
        return new CountingVariantFilter(
                filterSymbolicAndSV ?
                        VariantFilterLibrary.NOT_SV_OR_SYMBOLIC:
                        VariantFilterLibrary.ALLOW_ALL_VARIANTS
        );
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> readFilters = new ArrayList<>();
        readFilters.addAll(super.getDefaultReadFilters());
        List<String> filterList = new ArrayList<>();
        filterList.add("ID:" + HaplotypeBAMWriter.DEFAULT_HAPLOTYPE_READ_GROUP_ID);
        filterList.add("ID:" + HaplotypeBAMWriter.DEFAULT_GATK3_HAPLOTYPE_READ_GROUP_ID);
        readFilters.add(new ReadGroupBlackListReadFilter(filterList, null));
        return readFilters;
    }

    @Override
    public void onTraversalStart() {
        // Users can disable the AVX check to allow an older version of TF that doesn't require AVX to be used.
        if(this.disableAVXCheck == false) {
            IntelGKLUtils utils = new IntelGKLUtils();
            utils.load(null);
            if (utils.isAvxSupported() == false) {
                // Give user the bad news, suggest remedies.
                throw new UserException.HardwareFeatureException(String.format(CNNScoreVariants.AVXREQUIRED_ERROR, DISABLE_AVX_CHECK_NAME));
            }
        }

        final VCFHeader inputHeader = getHeaderForVariants();
        if (inputHeader.getGenotypeSamples().size() > 1) {
            logger.warn("CNNScoreVariants is a single sample tool but the input VCF has more than 1 sample.");
        }

        if (!annotationKeys.equals(defaultAnnotationKeys)){
            logger.warn("Annotation keys are not the default you must also provide a trained model that expects these annotations.");
        }

        // Start the Python process and initialize a stream writer for streaming data to the Python code
        pythonExecutor.start(Collections.emptyList(), enableJournal, pythonProfileResults);
        pythonExecutor.initStreamWriter(AsynchronousStreamWriter.stringSerializer);

        batchList = new ArrayList<>(transferBatchSize);

        // Execute Python code to open our output file, where it will write the contents of everything it reads
        // from the stream.
        try {
            // create a local temp that python code can write to
            scoreFile = File.createTempFile(outputFile.getBaseName().get(), ".temp");
            if (!keepTempFile) {
                scoreFile.deleteOnExit();
            } else {
                logger.info("Saving temp file from python:" + scoreFile.getAbsolutePath());
            }
            pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'w+')" + NL, scoreFile.getAbsolutePath()));
            pythonExecutor.sendSynchronousCommand("import vqsr_cnn" + NL);

            scoreKey = getScoreKeyAndCheckModelAndReadsHarmony();
            annotationSetString = annotationKeys.stream().collect(Collectors.joining(DATA_VALUE_SEPARATOR));
            initializePythonArgsAndModel();
        } catch (IOException e) {
            throw new GATKException("Error when creating temp file and initializing python executor.", e);
        }

    }

    @Override
    public void firstPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        referenceContext.setWindow(windowStart, windowEnd);
        if (tensorType.isReadsRequired()) {
            transferReadsToPythonViaFifo(variant, readsContext, referenceContext);
        } else {
            transferToPythonViaFifo(variant, referenceContext);
        }
        sendBatchIfReady();
    }

    @Override
    public void afterFirstPass() {
        if (waitforBatchCompletion) {
            pythonExecutor.waitForPreviousBatchCompletion();
        }
        if (curBatchSize > 0) {
            executePythonCommand();
            pythonExecutor.waitForPreviousBatchCompletion();
        }

        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);
        pythonExecutor.terminate();

        try {
            scoreScan = new Scanner(scoreFile);
            vcfWriter = createVCFWriter(outputFile);
            scoreScan.useDelimiter("\\n");
            writeVCFHeader(vcfWriter);
        } catch (IOException e) {
            throw new GATKException("Error when trying to temporary score file scanner.", e);
        }

    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        String sv = scoreScan.nextLine();
        String[] scoredVariant = sv.split("\\t");

        if (variant.getContig().equals(scoredVariant[CONTIG_INDEX])
                && Integer.toString(variant.getStart()).equals(scoredVariant[POS_INDEX])
                && variant.getReference().getBaseString().equals(scoredVariant[REF_INDEX])
                && variant.getAlternateAlleles().toString().equals(scoredVariant[ALT_INDEX])) {

            final VariantContextBuilder builder = new VariantContextBuilder(variant);
            if (scoredVariant.length > KEY_INDEX) {
                builder.attribute(scoreKey, scoredVariant[KEY_INDEX]);
            }
            vcfWriter.add(builder.make());

        } else {
            String errorMsg = "Score file out of sync with original VCF. Score file has:" + sv;
            errorMsg += "\n But VCF has:" + variant.toStringWithoutGenotypes();
            throw new GATKException(errorMsg);
        }
    }

    @Override
    public void closeTool() {
        logger.info("Done scoring variants with CNN.");
        if (vcfWriter != null) {
            vcfWriter.close();
        }
        if (scoreScan != null){
            scoreScan.close();
        }
    }

    private void transferToPythonViaFifo(final VariantContext variant, final ReferenceContext referenceContext) {
        try {
            final String outDat = String.format("%s%s%s%s%s%s%s\n",
                    getVariantDataString(variant), DATA_TYPE_SEPARATOR,
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"), DATA_TYPE_SEPARATOR,
                    getVariantInfoString(variant), DATA_TYPE_SEPARATOR,
                    variant.isSNP() ? "SNP" : variant.isIndel() ? "INDEL" : "OTHER");
            batchList.add(outDat);
            curBatchSize++;
        } catch (UnsupportedEncodingException e) {
            throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
        }

    }

    private void sendBatchIfReady() {
        if (curBatchSize == transferBatchSize) {
            if (waitforBatchCompletion == true) {
                // wait for the last batch to complete before we start a new one
                pythonExecutor.waitForPreviousBatchCompletion();
                waitforBatchCompletion = false;
            }
            executePythonCommand();
            waitforBatchCompletion = true;
            curBatchSize = 0;
            batchList = new ArrayList<>(transferBatchSize);
        }
    }

    private void transferReadsToPythonViaFifo(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext) {
        StringBuilder sb = new StringBuilder(FIFO_STRING_INITIAL_CAPACITY);
        try {
            sb.append(String.format("%s%s%s%s%s%s%s%s",
                    getVariantDataString(variant), DATA_TYPE_SEPARATOR,
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"), DATA_TYPE_SEPARATOR,
                    getVariantInfoString(variant), DATA_TYPE_SEPARATOR,
                    variant.isSNP() ? "SNP" : variant.isIndel() ? "INDEL" : "OTHER", DATA_TYPE_SEPARATOR));
        } catch (UnsupportedEncodingException e) {
            throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
        }
        Iterator<GATKRead> readIt = new ReadsDownsamplingIterator(readsContext.iterator(), new ReservoirDownsampler(readLimit));
        if (!readIt.hasNext()) {
            logger.warn("No reads at contig:" + variant.getContig() + " site:" + String.valueOf(variant.getStart()));
        }

        while (readIt.hasNext()) {
            sb.append(GATKReadToString(readIt.next()));
        }
        sb.append(NL);
        batchList.add(sb.toString());
        curBatchSize++;
    }

    private String GATKReadToString(GATKRead read) {
        StringBuilder sb = new StringBuilder(FIFO_STRING_INITIAL_CAPACITY);
        sb.append(read.getBasesString() + DATA_TYPE_SEPARATOR);

        appendQualityBytes(sb, read.getBaseQualities());
        sb.append(read.getCigar().toString() + DATA_TYPE_SEPARATOR);
        sb.append(read.isReverseStrand() + DATA_TYPE_SEPARATOR);
        sb.append((read.isPaired() ? read.mateIsReverseStrand() : "false") + DATA_TYPE_SEPARATOR);
        sb.append(read.isFirstOfPair() + DATA_TYPE_SEPARATOR);
        sb.append(read.getMappingQuality() + DATA_TYPE_SEPARATOR);
        sb.append(Integer.toString(read.getUnclippedStart()) + DATA_TYPE_SEPARATOR);
        return sb.toString();
    }

    private void appendQualityBytes(StringBuilder sb, byte[] qualities) {
        if(qualities.length == 0) {
            sb.append(DATA_TYPE_SEPARATOR);
            return;
        }

        for (int i = 0; i < qualities.length - 1; i++) {
            sb.append(Integer.toString(qualities[i]) + DATA_VALUE_SEPARATOR);
        }
        sb.append(Integer.toString(qualities[qualities.length - 1]) + DATA_TYPE_SEPARATOR);
    }

    private String getVariantDataString(final VariantContext variant) {
        return String.format("%s%s%d%s%s%s%s",
                variant.getContig(), DATA_TYPE_SEPARATOR,
                variant.getStart(), DATA_TYPE_SEPARATOR,
                variant.getReference().getBaseString(), DATA_TYPE_SEPARATOR,
                variant.getAlternateAlleles().toString()
        );
    }

    private String getVariantInfoString(final VariantContext variant) {
        // Create a string that will easily be parsed as a python dictionary
        StringBuilder sb = new StringBuilder(FIFO_STRING_INITIAL_CAPACITY);
        for (final String attributeKey : annotationKeys) {
            if (variant.hasAttribute(attributeKey)) {
                sb.append(attributeKey);
                sb.append(ANNOTATION_SET_STRING);
                sb.append(variant.getAttributeAsString(attributeKey, "0"));
                sb.append(ANNOTATION_SEPARATOR);
            }
        }
        return sb.toString();
    }

    private void executePythonCommand() {
        final String pythonCommand = String.format(
                "vqsr_cnn.score_and_write_batch(model, tempFile, %d, %d, '%s', '%s', %d, %d, '%s')",
                curBatchSize,
                inferenceBatchSize,
                tensorType,
                annotationSetString,
                windowSize,
                readLimit,
                outputTensorsDir) + NL;
        pythonExecutor.startBatchWrite(pythonCommand, batchList);
    }

    private void writeVCFHeader(VariantContextWriter vcfWriter) {
        // setup the header fields
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(scoreKey));
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    private String getScoreKeyAndCheckModelAndReadsHarmony() {
        if (tensorType.isReadsRequired() && this.hasReads()) {
            return GATKVCFConstants.CNN_2D_KEY;
        } else if (!tensorType.isReadsRequired() && this.hasReads()) {
            logger.warn(String.format("Reads are available, but tensor type %s does not use them.", tensorType.name()));
            return GATKVCFConstants.CNN_1D_KEY;
        } else if (!tensorType.isReadsRequired()) {
            return GATKVCFConstants.CNN_1D_KEY;
        } else {
            throw new GATKException("2D Models require a SAM/BAM file specified via -I (-input) argument.");
        }
    }

    private void initializePythonArgsAndModel() {
        if (architecture == null && weights == null) {
            if (tensorType.equals(TensorType.read_tensor)) {
                architecture = IOUtils.writeTempResourceFromPath(resourcePathReadTensor, null).getAbsolutePath();
                weights = IOUtils.writeTempResourceFromPath(
                        resourcePathReadTensor.replace(".json", ".hd5"),
                        null).getAbsolutePath();
            } else if (tensorType.equals(TensorType.reference)) {
                architecture = IOUtils.writeTempResourceFromPath(resourcePathReferenceTensor, null).getAbsolutePath();
                weights = IOUtils.writeTempResourceFromPath(
                        resourcePathReferenceTensor.replace(".json", ".hd5"), null).getAbsolutePath();
            } else {
                throw new GATKException("No default architecture for tensor type:" + tensorType.name());
            }
        } else if (weights == null) {
            weights = architecture.replace(".json", ".hd5");
        } else if (architecture == null) {
            architecture = weights.replace(".hd5", ".json");
        }

        String getArgsAndModel = String.format("args, model = vqsr_cnn.start_session_get_args_and_model(%d, %d, '%s', weights_hd5='%s')",
                intraOpThreads, interOpThreads, architecture, weights) + NL;
        logger.info("Using key:" + scoreKey + " for CNN architecture:" + architecture + " and weights:" + weights);
        pythonExecutor.sendSynchronousCommand(getArgsAndModel);
    }
}
