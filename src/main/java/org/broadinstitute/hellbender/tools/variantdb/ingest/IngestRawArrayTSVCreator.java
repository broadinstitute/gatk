package org.broadinstitute.hellbender.tools.variantdb.ingest;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;


import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class IngestRawArrayTSVCreator extends IngestTSVCreator {
    static final Logger logger = LogManager.getLogger(IngestRawArrayTSVCreator.class);

    private SimpleXSVWriter rawArrayWriter = null;
    private final String sampleId;

    public IngestRawArrayTSVCreator(String sampleName, String sampleId, Path sampleDirectoryPath) {
        this.sampleId = sampleId;
        // If the raw directory inside it doesn't exist yet -- create it
        final String rawDirectoryName = "raw";
        final Path rawDirectoryPath = sampleDirectoryPath.resolve(rawDirectoryName);
        final File rawDirectory = new File(rawDirectoryPath.toString());
        if (!rawDirectory.exists()) {
            rawDirectory.mkdir();
        }
        try {
            // Create a raw file to go into the raw dir for _this_ sample
            final String rawOutputName = sampleName + rawDirectoryName + IngestVariantWalker.FILETYPE;
            final Path rawOutputPath = rawDirectoryPath.resolve(rawOutputName);
            // write header to it
            List<String> rawHeader = IngestRawArrayTSVCreator.getHeaders();
            rawArrayWriter = new SimpleXSVWriter(rawOutputPath, IngestVariantWalker.SEPARATOR);
            rawArrayWriter.setHeaderLine(rawHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create raw outputs", e);
        }


    }

    public List<String> createRow(final long start, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        row.add(String.valueOf(start));
        row.add(sampleId);
        for ( final RawArrayFieldEnum fieldEnum : RawArrayFieldEnum.values() ) {
            if (!fieldEnum.equals(RawArrayFieldEnum.position) && !fieldEnum.equals(RawArrayFieldEnum.sample_id)) {
                row.add(fieldEnum.getColumnValue(variant));
            }        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(RawArrayFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final List<String> TSVLinesToCreate = createRow(SchemaUtils.encodeLocation(variant.getContig(), variant.getStart()), variant, sampleId);

        // write the row to the XSV
        SimpleXSVWriter.LineBuilder rawLine = rawArrayWriter.getNewLineBuilder();
        rawLine.setRow(TSVLinesToCreate);
        rawLine.write();
    }

    @Override
    public void closeTool() {
        if (rawArrayWriter != null) {
            try {
                rawArrayWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }
    }
}
