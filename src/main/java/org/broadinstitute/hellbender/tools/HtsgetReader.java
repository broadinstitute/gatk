package org.broadinstitute.hellbender.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.util.htsget.*;
import org.apache.commons.io.IOUtils;
import org.apache.http.client.utils.URIBuilder;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.OtherProgramGroup;

/**
 * A tool that downloads a file hosted on an htsget server to a local file
 * 
 * <h3>Usage example</h3>
 * <pre>
 * gatk HtsgetReader \
 *   --url htsget-server.org \
 *   --id A1.bam \
 *   --reference-name chr1
 *   -O output.bam
 * </pre>
 */

@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Download a file using htsget",
        oneLineSummary = "Download a file using htsget",
        programGroup = OtherProgramGroup.class
)
public class HtsgetReader extends CommandLineProgram {

    public static final String URL_LONG_NAME = "url";
    public static final String ID_LONG_NAME = "id";
    public static final String FORMAT_LONG_NAME = "format";
    public static final String CLASS_LONG_NAME = "class";
    public static final String FIELDS_LONG_NAME = "field";
    public static final String TAGS_LONG_NAME = "tag";
    public static final String NOTAGS_LONG_NAME = "notag";
    public static final String NUM_THREADS_LONG_NAME = "reader-threads";
    public static final String CHECK_MD5_LONG_NAME = "check-md5";

    @Argument(doc = "Output file.",
            fullName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    @Argument(doc = "URL of htsget endpoint.",
            fullName = URL_LONG_NAME,
            shortName = URL_LONG_NAME)
    private URI endpoint;

    @Argument(doc = "ID of record to request.",
            fullName = ID_LONG_NAME,
            shortName = ID_LONG_NAME)
    private String id;

    @Argument(doc = "Format to request record data in.",
            fullName = FORMAT_LONG_NAME,
            shortName = FORMAT_LONG_NAME,
            optional = true)
    private HtsgetFormat format;

    @Argument(doc = "Class of data to request.",
            fullName = CLASS_LONG_NAME,
            shortName = CLASS_LONG_NAME,
            optional = true)
    private HtsgetClass dataClass;

    @Argument(doc = "The interval and reference sequence to request",
            fullName = StandardArgumentDefinitions.INTERVALS_LONG_NAME,
            shortName = StandardArgumentDefinitions.INTERVALS_SHORT_NAME,
            optional = true)
    private SimpleInterval interval;

    @Argument(doc = "A field to include, default: all",
            fullName = FIELDS_LONG_NAME,
            shortName = FIELDS_LONG_NAME,
            optional = true)
    private List<HtsgetRequestField> fields;

    @Argument(doc = "A tag which should be included.",
            fullName = TAGS_LONG_NAME,
            shortName = TAGS_LONG_NAME,
            optional = true)
    private List<String> tags;

    @Argument(doc = "A tag which should be excluded.",
            fullName = NOTAGS_LONG_NAME,
            shortName = NOTAGS_LONG_NAME,
            optional = true)
    private List<String> notags;

    @Advanced
    @Argument(fullName = NUM_THREADS_LONG_NAME,
        shortName = NUM_THREADS_LONG_NAME,
        doc = "How many simultaneous threads to use when reading data from an htsget response;" +
            "higher values may improve performance when network latency is an issue.",
        optional = true,
        minValue = 1)
    private int readerThreads = 1;

    @Argument(fullName = CHECK_MD5_LONG_NAME, shortName = CHECK_MD5_LONG_NAME, doc = "Boolean determining whether to calculate the md5 digest of the assembled file "
            + "and validate it against the provided md5 hash, if it exists.", optional = true)
    private boolean checkMd5 = false;

    private ExecutorService executorService;

    private void checkMd5(final String expectedMd5) {
        try {
            final String actualMd5 = Utils.calculateFileMD5(this.outputFile);
            if (!actualMd5.equals(expectedMd5)) {
                throw new UserException("Expected md5: " + expectedMd5 + " and actual md5: " + actualMd5 + " do not match");
            }
        } catch (final IOException e) {
            throw new UserException("Could not calculate md5 checksum from downloaded file", e);
        }
    }

    @Override
    protected void onStartup() {
        if (this.readerThreads > 1) {
            this.executorService = Executors.newFixedThreadPool(this.readerThreads, new ThreadFactoryBuilder()
                .setNameFormat("htsget-reader-thread-%d")
                .setDaemon(true)
                .build());
        }
    }

    @Override
    protected void onShutdown() {
        if (this.executorService != null) {
            this.executorService.shutdownNow();
        }
    }

    @Override
    public Object doWork() {
        // Construct request from command line args and convert to URI
        final URI endpointWithId;
        try {
            final String endpointPath = endpoint.getPath();
            final String endpointPathWithSlash = endpointPath.endsWith("/") ? endpointPath : endpointPath.concat("/");
            endpointWithId = new URIBuilder(endpoint)
                .setPath(endpointPathWithSlash + id)
                .build();
        } catch (final URISyntaxException e) {
            throw new UserException(String.format("Error appending id: %s to provided endpoint: %s", id, endpoint), e);
        }

        final HtsgetRequest req = new HtsgetRequest(endpointWithId)
            .withFormat(format)
            .withDataClass(dataClass)
            .withInterval(interval)
            .withFields(fields)
            .withTags(tags)
            .withNotags(notags);

        final HtsgetResponse resp = req.getResponse();
        if (resp.getMd5() == null) {
            this.checkMd5 = false;
            logger.info("No md5 checksum received");
        }

        try (final OutputStream outputstream = new FileOutputStream(this.outputFile)) {
            if (this.readerThreads > 1) {
                final List<Future<InputStream>> futures = new ArrayList<>();
                resp.getBlocks().forEach(block -> futures.add(this.executorService.submit(() -> {
                    final Path tempFile = org.broadinstitute.hellbender.utils.io.IOUtils.createTempPath("htsget-temp", "");
                    try (final OutputStream ostream = Files.newOutputStream(tempFile)) {
                        org.apache.commons.io.IOUtils.copy(block.getData(), ostream);
                        return Files.newInputStream(tempFile);
                    } catch (final IOException e) {
                        throw new UserException("Error while downloading htsget block", e);
                    }
                })));
                try {
                    for (final Future<InputStream> future : futures) {
                        IOUtils.copy(future.get(), outputstream);
                    }
                } catch (final IOException e) {
                    throw new UserException("IOException while writing output file", e);
                } catch (final InterruptedException | ExecutionException e) {
                    throw new UserException("Interrupted while writing output file", e);
                }
            } else {
                try {
                    IOUtils.copy(resp.getDataStream(), outputstream);
                } catch (final IOException e) {
                    throw new UserException("IOException while writing output file", e);
                }
            }
        } catch (final IOException e) {
            throw new UserException("IOException during htsget download", e);
        }

        if (this.checkMd5) {
            this.checkMd5(resp.getMd5());
        }

        return null;
    }
}
