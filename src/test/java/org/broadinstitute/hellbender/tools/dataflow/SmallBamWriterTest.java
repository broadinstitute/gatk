package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.SmallBamWriter;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.channels.Channels;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests for SmallBamWriter.
 */
public class SmallBamWriterTest extends BaseTest {
    private static final String LOCAL_INPUT = "src/test/resources/org/broadinstitute/hellbender/engine/reads_data_source_test2.bam";
    private static final String THIS_TEST_FOLDER = "org/broadinstitute/hellbender/engine/reads_data_source_test2.bam";

    private SAMFileHeader header;

    @Test
    public void checkLocal() throws Exception {
        File out = createTempFile("temp",".bam");
        String outputFile = out.getPath();
        testReadAndWrite(LOCAL_INPUT, outputFile, false, false);
    }

    @Test(groups = {"bucket"})
    public void checkGCSInput() throws Exception {
        File out = createTempFile("temp",".bam");
        String outputFile = out.getPath();
        testReadAndWrite(getCloudInput(), outputFile, true, false);
    }

    // this leaves output files around, for now.
    @Test(groups = {"bucket"})
    public void checkGCSOutput() throws Exception {
        File out = createTempFile("temp",".bam");
        String tempName = out.getName();
        String outputPath = getDataflowTestStaging() + tempName;
        testReadAndWrite(LOCAL_INPUT, outputPath, true, false);
    }

    // this leaves output files around, for now.
    @Test(groups = {"cloud"})
    public void checkCloud() throws Exception {
        File out = createTempFile("temp",".bam");
        String tempName = out.getName();
        String outputPath = getDataflowTestStaging() + tempName;
        testReadAndWrite(getCloudInput(), outputPath, true, true);
    }

    protected void testReadAndWrite(final String inputPath, final String outputPath, boolean enableGcs, boolean enableCloudExec) throws IOException {
        String localInput = inputPath;
        String localOutput = outputPath;
        logger.info("outputPath="+outputPath);
        final Pipeline pipeline = setupPipeline(enableGcs, enableCloudExec);
        if (BucketUtils.isCloudStorageUrl(inputPath)) {
            File out = createTempFile("temp-input",".bam");
            localInput = out.getPath();
            logger.info("Downloading the input to "+localInput+" (for comparison).");
            downloadFromGCS(pipeline.getOptions(), inputPath, localInput);
        }
        PCollection<Read> input = ingestReadsAndGrabHeader(pipeline, inputPath);
        SmallBamWriter.writeToFile(pipeline, input, header, outputPath);
        pipeline.run();
        if (BucketUtils.isCloudStorageUrl(outputPath)) {
            File out = createTempFile("temp-output",".bam");
            localOutput = out.getPath();
            logger.info("Downloading the output to " + localOutput+" (for comparison).");
            downloadFromGCS(pipeline.getOptions(), outputPath, localOutput);
        }
        IntegrationTestSpec.compareBamFiles(new File(localInput), new File(localOutput));
    }

    private void downloadFromGCS(PipelineOptions popts, String gcsPath, String localPath) throws IOException {
        try (
                InputStream in = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(popts).open(GcsPath.fromUri(gcsPath)));
                FileOutputStream fout = new FileOutputStream(localPath)) {
            final byte[] buf = new byte[1024 * 1024];
            int count;
            while ((count = in.read(buf)) > 0) {
                fout.write(buf, 0, count);
            }
        }
    }

    private Pipeline setupPipeline(boolean enableGcs, boolean enableCloudExec) {
        final GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
        if (enableCloudExec) {
            options.setStagingLocation(getDataflowTestStaging());
            options.setProject(getDataflowTestProject());
            options.setRunner(BlockingDataflowPipelineRunner.class);
        } else {
            options.setRunner(DirectPipelineRunner.class);
        }
        if (enableGcs) {
            options.setApiKey(getDataflowTestApiKey());
        }
        final Pipeline p = Pipeline.create(options);
        DataflowWorkarounds.registerGenomicsCoders(p);
        return p;
    }

    /** reads local disks or GCS -> header, and PCollection */
    private PCollection<Read> ingestReadsAndGrabHeader(final Pipeline pipeline, String filename) throws IOException {

        // input reads
        if (BucketUtils.isCloudStorageUrl(filename)) {
            // set up ingestion on the cloud
            // but read the header locally
            GcsPath path = GcsPath.fromUri(filename);
            InputStream inputstream = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(pipeline.getOptions())
                    .open(path));
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(inputstream));
            header = reader.getFileHeader();

            final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
            return new ReadsSource(filename, pipeline).getReadPCollection(intervals, ValidationStringency.SILENT);
        } else {
            // ingestion from local file
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(filename));
            header = reader.getFileHeader();
            List<Read> readLst = new ArrayList<>();
            for (SAMRecord sr : reader) {
                try {
                    Read e = ReadConverter.makeRead(sr);
                    readLst.add(e);
                } catch (SAMException x) {
                    logger.warn("Skipping read " + sr.getReadName() + " because we can't convert it.");
                } catch (NullPointerException y) {
                    logger.warn("Skipping read " + sr.getReadName() + " because we can't convert it. (null?)");
                }
            }
            return pipeline.apply(Create.of(readLst).withName("input ingest"));
        }
    }

    private String getCloudInput() {
        return getDataflowTestInputs() + THIS_TEST_FOLDER;
    }
}
