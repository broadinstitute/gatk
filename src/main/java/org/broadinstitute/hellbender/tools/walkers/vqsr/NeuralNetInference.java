package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriterService;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.StreamSupport;

/**
 * Annotate a VCF with scores from 1D Convolutional Neural Network.
 *
 * This tool streams variants and their reference context to a python program
 * which evaluates a pretrained neural network on each variant.
 * The neural network performs 1D convolutions over the reference sequence surrounding the variant
 * and combines those features with a multilayer perceptron on the variant annotations.
 *
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk NeuralNetInference\
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   --architecture src/main/resources/org/broadinstitute/hellbender/tools/walkers/vqsr/cnn_1d_annotations.hd5
 * </pre>
 *
 */
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = NeuralNetInference.USAGE_SUMMARY,
        oneLineSummary = NeuralNetInference.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)

public class NeuralNetInference extends VariantWalker {
    private final static String NL = String.format("%n");
    static final String USAGE_ONE_LINE_SUMMARY = "Apply 1d Convolutional Neural Net to filter annotated variants";
    static final String USAGE_SUMMARY = "Annotate a VCF with scores from 1d Convolutional Neural Network (CNN)." +
            "The CNN will look at the reference sequence and variant annotations to determine a Log Odds Score for each variant.";

    private static final int CONTIG_INDEX = 0;
    private static final int POS_INDEX = 1;
    private static final int REF_INDEX = 2;
    private static final int ALT_INDEX = 3;
    private static final int KEY_INDEX = 4;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private String outputFile = null; // output file produced by Python code

    @Argument(fullName = "architecture", shortName = "a", doc = "Neural Net architecture and weights hd5 file", optional = false)
    private String architecture = null;

    @Argument(fullName = "window-size", shortName = "ws", doc = "Neural Net input window size", minValue = 0, optional = true)
    private int windowSize = 128;

    @Advanced
    @Argument(fullName = "inference-batch-size", shortName = "ibs", doc = "Size of batches for python to do inference on.", minValue = 1, maxValue=4096, optional = true)
    private int inferenceBatchSize = 256;

    @Advanced
    @Argument(fullName = "transfer-batch-size", shortName = "tbs", doc = "Size of data to queue for python streaming.", minValue = 1, maxValue=8192, optional = true)
    private int transferBatchSize = 512;

    @Hidden
    @Argument(fullName = "enable-journal", shortName = "journal", doc = "Enable streaming process journal.", optional = true)
    private boolean enableJournal = false;

    @Hidden
    @Argument(fullName = "keep-temp-file", shortName = "ktf", doc = "Keep the temporary file that python writes scores to.", optional = true)
    private boolean keepTempFile = false;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final StreamingPythonScriptExecutor pythonExecutor = new StreamingPythonScriptExecutor(true);

    private FileOutputStream fifoWriter;
    private AsynchronousStreamWriterService<String> asyncWriter = null;
    private List<String> batchList = new ArrayList<>(inferenceBatchSize);

    private int curBatchSize = 0;
    private int windowEnd = windowSize/2;
    private int windowStart = (windowSize/2)-1;
    private boolean waitforBatchCompletion = false;
    private File scoreFile;

    @Override
    protected String[] customCommandLineValidation(){
        if(inferenceBatchSize > transferBatchSize){
            String[] errors = {"Inference batch size must be less than or equal to transfer batch size."};
            return errors;
        }
        return null;
    }

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public void onTraversalStart() {
        // Start the Python process, and get a FIFO from the executor to use to send data to Python. The lifetime
        // of the FIFO is managed by the executor; the FIFO will be destroyed when the executor is destroyed.
        pythonExecutor.start(Collections.emptyList(), enableJournal);
        final File fifoFile = pythonExecutor.getFIFOForWrite();

        // Open the FIFO for writing. Opening a FIFO for read or write will block until there is reader/writer
        // on the other end, so before we open it, send an ASYNCHRONOUS command (one that doesn't wait for a
        // response) to the Python process to open the FIFO for reading. The Python process will then block until
        // we open the FIFO. We can then call getAccumulatedOutput.
        pythonExecutor.sendAsynchronousCommand(String.format("fifoFile = open('%s', 'r')" + NL, fifoFile.getAbsolutePath()));
        try {
            fifoWriter = new FileOutputStream(fifoFile);
        } catch ( IOException e ) {
            throw new GATKException("Failure opening FIFO for writing", e);
        }

        pythonExecutor.getAccumulatedOutput();
        asyncWriter = pythonExecutor.getAsynchronousStreamWriterService(fifoWriter, AsynchronousStreamWriterService.stringSerializer);
        batchList = new ArrayList<>(transferBatchSize);

        // Also, ask Python to open our output file, where it will write the contents of everything it reads
        // from the FIFO. <code sendSynchronousCommand/>
        try {
            scoreFile = File.createTempFile(outputFile, ".temp");
            if (!keepTempFile){
                scoreFile.deleteOnExit();
            }
            pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'w+')" + NL, scoreFile.getAbsolutePath()));
            pythonExecutor.sendSynchronousCommand("from keras.models import load_model" + NL);
            pythonExecutor.sendSynchronousCommand("import vqsr_cnn" + NL);
            pythonExecutor.sendSynchronousCommand(String.format("model = load_model('%s', custom_objects=vqsr_cnn.get_metric_dict())", architecture) + NL);
            logger.info("Loaded CNN architecture:"+architecture);
        } catch (IOException e) {
            throw new GATKException("Error when creating temp file and initializing python executor.", e);
        }

    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        referenceContext.setWindow(windowStart, windowEnd);
        transferToPythonViaFifo(variant, referenceContext);
    }

    private void transferToPythonViaFifo(final VariantContext variant, final ReferenceContext referenceContext) {
        try {
            final String outDat = String.format("%s\t%s\t%s\t%s\n",
                    getVariantDataString(variant),
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"),
                    getVariantInfoString(variant),
                    variant.isSNP() ? "SNP" : variant.isIndel() ? "INDEL" : "OTHER");

            if (curBatchSize == transferBatchSize) {
                if (waitforBatchCompletion == true) {
                    // wait for the last batch to complete before we start a new one
                    asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
                    waitforBatchCompletion = false;
                    pythonExecutor.getAccumulatedOutput();
                }
                executePythonCommand();
                waitforBatchCompletion = true;
                curBatchSize = 0;
                batchList = new ArrayList<>(transferBatchSize);
            }

            batchList.add(outDat);
            curBatchSize++;
        } catch (UnsupportedEncodingException e) {
            throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
        }

    }

    private String getVariantDataString(final VariantContext variant){
        return String.format("%s\t%d\t%s\t%s",
                variant.getContig(),
                variant.getStart(),
                variant.getReference().getBaseString(),
                variant.getAlternateAlleles().toString()
        );

    }

    private String getVariantInfoString(final VariantContext variant){
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
            asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
            pythonExecutor.getAccumulatedOutput();
        }
        if (curBatchSize > 0){
            executePythonCommand();
            asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
            pythonExecutor.getAccumulatedOutput();
        }

        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);
        pythonExecutor.sendSynchronousCommand("fifoFile.close()" + NL);
        pythonExecutor.terminate();

        writeOutputVCFWithScores();

        return true;
    }

    private void executePythonCommand(){
        final String pythonCommand = String.format(
                "vqsr_cnn.score_and_write_batch(model, tempFile, fifoFile, %d, %d)", curBatchSize, inferenceBatchSize) + NL;
        pythonExecutor.sendAsynchronousCommand(pythonCommand);
        asyncWriter.startAsynchronousBatchWrite(batchList);
    }


    private void writeOutputVCFWithScores(){


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
                        if(variant.getContig().equals(scoredVariant[CONTIG_INDEX])
                                && Integer.toString(variant.getStart()).equals(scoredVariant[POS_INDEX])
                                && variant.getReference().getBaseString().equals(scoredVariant[REF_INDEX])
                                && variant.getAlternateAlleles().toString().equals(scoredVariant[ALT_INDEX])) {
                            final VariantContextBuilder builder = new VariantContextBuilder(variant);
                            builder.attribute(GATKVCFConstants.CNN_1D_KEY, scoredVariant[KEY_INDEX]);
                            vcfWriter.add(builder.make());
                        } else {
                            String errorMsg = "Score file out of sync with original VCF. Score file has:"+sv;
                            errorMsg += "\n But VCF has:" + variant.toStringWithoutGenotypes();
                            throw new GATKException(errorMsg);
                        }
                    });

        } catch (IOException e) {
            throw new GATKException("Error when trying to write annotated VCF.", e);
        }

    }

    private void writeVCFHeader(VariantContextWriter vcfWriter){
        // setup the header fields
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CNN_1D_KEY));
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

}
