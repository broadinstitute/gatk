package org.broadinstitute.hellbender.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.MapperFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.json.JsonMapper;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import org.apache.commons.io.Charsets;
import org.apache.commons.io.IOUtils;
import org.apache.http.Header;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.util.EntityUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.htsgetreader.HtsgetClass;
import org.broadinstitute.hellbender.tools.htsgetreader.HtsgetErrorResponse;
import org.broadinstitute.hellbender.tools.htsgetreader.HtsgetFormat;
import org.broadinstitute.hellbender.tools.htsgetreader.HtsgetRequestBuilder;
import org.broadinstitute.hellbender.tools.htsgetreader.HtsgetRequestField;
import org.broadinstitute.hellbender.tools.htsgetreader.HtsgetResponse;
import org.broadinstitute.hellbender.utils.HttpUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

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
        programGroup = ExampleProgramGroup.class
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

    private CloseableHttpClient client;

    @Override
    public void onStartup() {
        if (this.readerThreads > 1) {
            logger.info("Initializing with " + this.readerThreads + " threads");
            final ThreadFactory threadFactory = new ThreadFactoryBuilder()
                .setNameFormat("htsgetReader-thread-%d")
                .setDaemon(true).build();
            this.executorService = Executors.newFixedThreadPool(readerThreads, threadFactory);
        }
        this.client = HttpUtils.getClient();
    }

    @Override
    public void onShutdown() {
        if (this.executorService != null) {
            this.executorService.shutdownNow();
        }
        super.onShutdown();
    }

    /**
     * Downloads data blocks provided by response to outputFile in serial
     */
    private void getData(final HtsgetResponse response) {
        try (final OutputStream ostream = new FileOutputStream(this.outputFile)) {
            response.getBlocks().forEach(b -> {
                try (final InputStream istream = b.getData()) {
                    IOUtils.copy(istream, ostream);
                } catch (final IOException e) {
                    throw new UserException("Failed to copy data block to output file", e);
                }
            });
        } catch (final IOException e) {
            throw new UserException("Could not create output file: " + outputFile, e);
        }
    }

    /**
     * Downloads data blocks provided by response to outputFile in parallel, using
     * the number of threads specified by user
     */
    private void getDataParallel(final HtsgetResponse response) {
        final List<Future<InputStream>> futures = new ArrayList<>(response.getBlocks().size());
        response.getBlocks().forEach(b -> futures.add(this.executorService.submit(b::getData)));

        try (final OutputStream ostream = new FileOutputStream(this.outputFile)) {
            futures.forEach(f -> {
                try (final InputStream istream = f.get()) {
                    IOUtils.copy(istream, ostream);
                } catch (final IOException e) {
                    throw new UserException("Error while copying data block to output file", e);
                } catch (final ExecutionException | InterruptedException e) {
                    throw new UserException("Error while waiting to download block", e);
                }
            });
        } catch (final IOException e) {
            throw new UserException("Could not create output file", e);
        }
    }

    /**
     * Checks md5 digest provided in response, if one exists, against calculated md5
     * hash of downloaded file, warning user if they differ
     */
    private void checkMd5(final HtsgetResponse resp) {
        final String expectedMd5 = resp.getMd5();
        if (expectedMd5 == null) {
            logger.warn("No md5 digest provided by response");
        } else {
            try {
                final String actualMd5 = Utils.calculateFileMD5(outputFile);
                if (!actualMd5.equals(expectedMd5)) {
                    throw new UserException("Expected md5: " + expectedMd5 + " did not match actual md5: " + actualMd5);
                }
            } catch (final IOException e) {
                throw new UserException("Unable to calculate md5 digest", e);
            }
        }
    }

    private JsonMapper getObjectMapper() {
        return JsonMapper.builder()
                .enable(DeserializationFeature.UNWRAP_ROOT_VALUE)
                .configure(MapperFeature.ACCEPT_CASE_INSENSITIVE_PROPERTIES, true)
                .build();
    }

    @Override
    public Object doWork() {
        // construct request from command line args and convert to URI
        final HtsgetRequestBuilder req = new HtsgetRequestBuilder(endpoint, id)
            .withFormat(format)
            .withDataClass(dataClass)
            .withInterval(interval)
            .withFields(fields)
            .withTags(tags)
            .withNotags(notags);
        final URI reqURI = req.toURI();

        final HttpGet getReq = new HttpGet(reqURI);
        try (final CloseableHttpResponse resp = this.client.execute(getReq)) {
            // get content of response
            final HttpEntity entity = resp.getEntity();
            final Header encodingHeader = entity.getContentEncoding();
            final Charset encoding = encodingHeader == null 
                ? StandardCharsets.UTF_8
                : Charsets.toCharset(encodingHeader.getValue());
            final String jsonBody = EntityUtils.toString(entity, encoding);

            final ObjectMapper mapper = this.getObjectMapper();

            if (resp.getStatusLine() == null) {
                throw new UserException(String.format("htsget server response did not contain status line for request %s", reqURI));
            }
            final int statusCode = resp.getStatusLine().getStatusCode();
            if (400 <= statusCode && statusCode < 500) {
                final HtsgetErrorResponse err = mapper.readValue(jsonBody, HtsgetErrorResponse.class);
                throw new UserException(String.format("Invalid request %s, received error code: %d, error type: %s, message: %s",
                        reqURI,
                        statusCode,
                        err.getError(),
                        err.getMessage()));
            } else if (statusCode == 200) {
                final HtsgetResponse response = mapper.readValue(jsonBody, HtsgetResponse.class);

                if (this.readerThreads > 1) {
                    this.getDataParallel(response);
                } else {
                    this.getData(response);
                }

                logger.info("Successfully wrote to: " + outputFile);

                if (checkMd5) {
                    this.checkMd5(response);
                }
            } else {
                throw new UserException(String.format("Unrecognized status code: %d for request %s", statusCode, reqURI));
            }
        } catch (final IOException e) {
            throw new UserException(String.format("IOException during htsget download for %s", reqURI), e);
        }
        return null;
    }
}
