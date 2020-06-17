package org.broadinstitute.hellbender.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.htsgetreader.*;
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
    public static final String PARALLEL_DOWNLOAD_LONG_NAME = "parallel";
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
    @Argument(doc = "Whether to try to download blocks in parallel",
        fullName = PARALLEL_DOWNLOAD_LONG_NAME,
        shortName = PARALLEL_DOWNLOAD_LONG_NAME,
        optional = true)
    private final boolean parallelDownload = false;

    @Argument(fullName = CHECK_MD5_LONG_NAME, shortName = CHECK_MD5_LONG_NAME, doc = "Boolean determining whether to calculate the md5 digest of the assembled file "
            + "and validate it against the provided md5 hash, if it exists.", optional = true)
    private boolean checkMd5 = false;

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
    public Object doWork() {
        // Construct request from command line args and convert to URI
        final HtsgetRequest req = new HtsgetRequest(endpoint, id)
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
            if (this.parallelDownload) {
                resp.getBlocks().parallelStream().map(block -> {
                        final Path tempFile = org.broadinstitute.hellbender.utils.io.IOUtils.createTempPath("htsget-temp", "");
                        try (final OutputStream ostream = Files.newOutputStream(tempFile)) {
                            org.apache.commons.io.IOUtils.copy(block.getData(), ostream);
                            return Files.newInputStream(tempFile);
                        } catch (final IOException e) {
                            throw new UserException("Error while downloading htsget block", e);
                        }
                    }).forEachOrdered(inputStream -> {
                        try {
                            IOUtils.copy(inputStream, outputstream);
                        } catch (final IOException e) {
                            throw new UserException("IOException while writing output file", e);
                        }
                    });
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

        if (this.checkMd5) this.checkMd5(resp.getMd5());

        return null;
    }
}
