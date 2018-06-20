package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import  org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.*;
import java.util.*;
import java.util.stream.StreamSupport;


/**
 * Annotate a VCF with scores from a Convolutional Neural Network (CNN).
 *
 * This tool streams variants and their reference context to a python program
 * which evaluates a pre-trained neural network on each variant.
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
 *   -inference-batch-size 2 \
 *   -transfer-batch-size 2 \
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
 *   -architecture path/to/my_model.json \
 *   -weights path/to/my_weights.hd5
 * </pre>
 *
 * <h3>2D Model with user-supplied architecture and weights:</h3>
 *
 * <pre>
 * gatk CNNScoreVariants \
 *   -I aligned_reads.bam \
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   -inference-batch-size 2 \
 *   -transfer-batch-size 2 \
 *   -tensor-type read-tensor \
 *   -architecture path/to/my_model.json \
 *   -weights path/to/my_weights.hd5
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = CNNScoreVariants.USAGE_SUMMARY,
        oneLineSummary = CNNScoreVariants.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantFilteringProgramGroup.class
)

public class CNNScoreVariants extends VariantWalker {
    private final static String NL = String.format("%n");
    static final String USAGE_ONE_LINE_SUMMARY = "Apply a Convolutional Neural Net to filter annotated variants";
    static final String USAGE_SUMMARY = "Annotate a VCF with scores from a Convolutional Neural Network (CNN)." +
            "The CNN determines a Log Odds Score for each variant." +
            "Pre-trained models (1D or 2D) are specified via the architecture argument." +
            "1D models will look at the reference sequence and variant annotations." +
            "2D models look at aligned reads, reference sequence, and variant annotations." +
            "2D models require a BAM file as input as well as the tensor-type argument to be set.";

    private static final int CONTIG_INDEX = 0;
    private static final int POS_INDEX = 1;
    private static final int REF_INDEX = 2;
    private static final int ALT_INDEX = 3;
    private static final int KEY_INDEX = 4;
    private static final int FIFO_STRING_INITIAL_CAPACITY = 1024;
    private static final int MAX_READ_BATCH = 4098;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private String outputFile;

    @Argument(fullName = "architecture", shortName = "architecture", doc = "Neural Net architecture configuration json file", optional = true)
    private String architecture;

    @Argument(fullName = "weights", shortName = "weights", doc = "Keras model HD5 file with neural net weights.", optional = true)
    private String weights;

    @Argument(fullName = "tensor-type", shortName = "tensor-type", doc = "Name of the tensors to generate, reference for 1D reference tensors and read_tensor for 2D tensors.", optional = true)
    private TensorType tensorType = TensorType.reference;

    @Argument(fullName = "window-size", shortName = "window-size", doc = "Neural Net input window size", minValue = 0, optional = true)
    private int windowSize = 128;

    @Argument(fullName = "filter-symbolic-and-sv", shortName = "filter-symbolic-and-sv", doc = "If set will filter symbolic and and structural variants from the input VCF", optional = true)
    private boolean filterSymbolicAndSV = false;

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
    private File scoreFile;

    private String scoreKey;

    private static String resourcePathReadTensor = Resource.LARGE_RUNTIME_RESOURCES_PATH + "/cnn_score_variants/small_2d.json";
    private static String resourcePathReferenceTensor = Resource.LARGE_RUNTIME_RESOURCES_PATH + "/cnn_score_variants/1d_cnn_mix_train_full_bn.json";

    @Override
    protected String[] customCommandLineValidation() {
        if (inferenceBatchSize > transferBatchSize) {
            return new String[]{"Inference batch size must be less than or equal to transfer batch size."};
        }

        if (weights == null && architecture == null){
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
    protected VariantFilter makeVariantFilter(){
        if (filterSymbolicAndSV) {
            return VariantFilterLibrary.NOT_SV_OR_SYMBOLIC;
        } else {
            return VariantFilterLibrary.ALLOW_ALL_VARIANTS;
        }
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
        scoreKey = getScoreKeyAndCheckModelAndReadsHarmony();
        if (architecture == null && weights == null) {
            setArchitectureAndWeightsFromResources();
        }

        // Start the Python process and initialize a stream writer for streaming data to the Python code
        pythonExecutor.start(Collections.emptyList(), enableJournal, pythonProfileResults);

        pythonExecutor.initStreamWriter(AsynchronousStreamWriter.stringSerializer);
        batchList = new ArrayList<>(transferBatchSize);

        // Execute Python code to open our output file, where it will write the contents of everything it reads
        // from the stream.
        try {
            scoreFile = File.createTempFile(outputFile, ".temp");
            if (!keepTempFile) {
                scoreFile.deleteOnExit();
            } else {
                logger.info("Saving temp file from python:" + scoreFile.getAbsolutePath());
            }

            pythonExecutor.sendSynchronousCommand("from keras import backend" + NL);
            pythonExecutor.sendSynchronousCommand(String.format("backend.set_session(backend.tf.Session(config=backend.tf.ConfigProto(intra_op_parallelism_threads=%d, inter_op_parallelism_threads=%d)))" + NL, intraOpThreads, interOpThreads));

            pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'w+')" + NL, scoreFile.getAbsolutePath()));
            pythonExecutor.sendSynchronousCommand("import vqsr_cnn" + NL);

            String getArgsAndModel;
            if (weights != null && architecture != null) {
                getArgsAndModel = String.format("args, model = vqsr_cnn.args_and_model_from_semantics('%s', weights_hd5='%s')", architecture, weights) + NL;
                logger.info("Using key:" + scoreKey + " for CNN architecture:" + architecture + " and weights:" + weights);
            } else if (architecture == null) {
                getArgsAndModel = String.format("args, model = vqsr_cnn.args_and_model_from_semantics(None, weights_hd5='%s', tensor_type='%s')", weights, tensorType.name()) + NL;
                logger.info("Using key:" + scoreKey + " for CNN weights:" + weights);
            } else {
                getArgsAndModel = String.format("args, model = vqsr_cnn.args_and_model_from_semantics('%s')", architecture) + NL;
                logger.info("Using key:" + scoreKey + " for CNN architecture:" + architecture);
            }
            pythonExecutor.sendSynchronousCommand(getArgsAndModel);

        } catch (IOException e) {
            throw new GATKException("Error when creating temp file and initializing python executor.", e);
        }

    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        referenceContext.setWindow(windowStart, windowEnd);
        if (tensorType.isReadsRequired()) {
            transferReadsToPythonViaFifo(variant, readsContext, referenceContext);
        } else {
            transferToPythonViaFifo(variant, referenceContext);
        }
        sendBatchIfReady();
    }

    private void transferToPythonViaFifo(final VariantContext variant, final ReferenceContext referenceContext) {
        try {
            final String outDat = String.format("%s\t%s\t%s\t%s\n",
                    getVariantDataString(variant),
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"),
                    getVariantInfoString(variant),
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
            sb.append(String.format("%s\t%s\t%s\t%s\t",
                    getVariantDataString(variant),
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"),
                    getVariantInfoString(variant),
                    variant.isSNP() ? "SNP" : variant.isIndel() ? "INDEL" : "OTHER"));
        } catch (UnsupportedEncodingException e) {
            throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
        }
        Iterator<GATKRead> readIt = readsContext.iterator();
        if (!readIt.hasNext()) {
            logger.warn("No reads at contig:" + variant.getContig() + " site:" + String.valueOf(variant.getStart()));
        }
        while (readIt.hasNext()) {
            sb.append(GATKReadToString(readIt.next()));
        }
        sb.append("\n");
        batchList.add(sb.toString());
        curBatchSize++;
    }

    private String GATKReadToString(GATKRead read) {
        StringBuilder sb = new StringBuilder(FIFO_STRING_INITIAL_CAPACITY);
        sb.append(read.getBasesString() + "\t");

        appendQualityBytes(sb, read.getBaseQualities());
        sb.append(read.getCigar().toString() + "\t");
        sb.append(read.isReverseStrand() + "\t");
        sb.append((read.isPaired() ? read.mateIsReverseStrand() : "false") + "\t");
        sb.append(read.isFirstOfPair() + "\t");
        sb.append(read.getMappingQuality() + "\t");
        sb.append(Integer.toString(read.getUnclippedStart()) + "\t");
        return sb.toString();
    }

    private void appendQualityBytes(StringBuilder sb, byte[] qualities) {
        if(qualities.length == 0) {
            sb.append("\t");
            return;
        }

        for (int i = 0; i < qualities.length - 1; i++) {
            sb.append(Integer.toString(qualities[i]) + ",");
        }
        sb.append(Integer.toString(qualities[qualities.length - 1]) + "\t");
    }

    private String getVariantDataString(final VariantContext variant) {
        return String.format("%s\t%d\t%s\t%s",
                variant.getContig(),
                variant.getStart(),
                variant.getReference().getBaseString(),
                variant.getAlternateAlleles().toString()
        );
    }

    private String getVariantInfoString(final VariantContext variant) {
        // Create a string that will easily be parsed as a python dictionary
        String varInfo = "";
        for (final String attributeKey : variant.getAttributes().keySet()) {
            varInfo += attributeKey + "=" + variant.getAttribute(attributeKey).toString().replace(" ", "").replace("[", "").replace("]", "") + ";";
        }
        return varInfo;
    }

    @Override
    public Object onTraversalSuccess() {
        if (waitforBatchCompletion) {
            pythonExecutor.waitForPreviousBatchCompletion();
        }
        if (curBatchSize > 0) {
            executePythonCommand();
            pythonExecutor.waitForPreviousBatchCompletion();
        }

        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);
        pythonExecutor.terminate();

        writeOutputVCFWithScores();

        return true;
    }

    private void executePythonCommand() {
        final String pythonCommand = String.format(
                "vqsr_cnn.score_and_write_batch(args, model, tempFile, %d, %d, '%s')",
                curBatchSize,
                inferenceBatchSize,
                outputTensorsDir) + NL;
        pythonExecutor.startBatchWrite(pythonCommand, batchList);
    }


    private void writeOutputVCFWithScores() {
        try (final Scanner scoreScan = new Scanner(scoreFile);
             final VariantContextWriter vcfWriter = createVCFWriter(new File(outputFile))) {
            scoreScan.useDelimiter("\\n");
            writeVCFHeader(vcfWriter);
            final VariantFilter variantfilter = makeVariantFilter();

            // Annotate each variant in the input stream, as in variantWalkerBase.traverse()
            StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                    .filter(variantfilter)
                    .forEach(variant -> {
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
                    });

        } catch (IOException e) {
            throw new GATKException("Error when trying to write annotated VCF.", e);
        }

    }

    private void writeVCFHeader(VariantContextWriter vcfWriter) {
        // setup the header fields
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(scoreKey));
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
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

    private void setArchitectureAndWeightsFromResources() {
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
    }

}

