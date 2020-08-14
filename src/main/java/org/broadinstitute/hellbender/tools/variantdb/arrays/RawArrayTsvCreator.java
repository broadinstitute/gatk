package org.broadinstitute.hellbender.tools.variantdb.arrays;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;


import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class RawArrayTsvCreator {
    static final Logger logger = LogManager.getLogger(RawArrayTsvCreator.class);

    private SimpleXSVWriter rawArrayWriter = null;
    private final String sampleId;
    private final Map<String, ProbeInfo> probeDataByName;
    private final boolean useCompressedData;
    private static String RAW_FILETYPE_PREFIX = "raw_";

    enum GT_encoding {
        AA("AA"),
        AB("AB"),
        BB("BB"),
        MISSING(".");

        String value;
        GT_encoding(String v) {
            value = v;
        }
        String getValue() {
            return value;
        }
    }

    public RawArrayTsvCreator(final String sampleName, final String sampleId, final String tableNumberPrefix, final Map<String, ProbeInfo> probeDataByName, boolean useCompressedData) {
        this.sampleId = sampleId;
        this.probeDataByName = probeDataByName;
        this.useCompressedData = useCompressedData;
        try {
            // Create a raw file to go into the raw dir for _this_ sample
            final String rawOutputName = RAW_FILETYPE_PREFIX + tableNumberPrefix + sampleName  + IngestConstants.FILETYPE;
            // write header to it
            List<String> rawHeader = RawArrayTsvCreator.getHeaders();
            rawArrayWriter = new SimpleXSVWriter(Paths.get(rawOutputName), IngestConstants.SEPARATOR);
            rawArrayWriter.setHeaderLine(rawHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create raw outputs", e);
        }
    }

    public List<String> createRow(final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        String rsid = variant.getID();
        if (rsid == null) {
            throw new IllegalStateException("Cannot be missing required value for site_name"); // TODO, should this be UserException too?
        }
        ProbeInfo probeInfo = probeDataByName.get(rsid);
        if (probeInfo == null) {
            logger.warn("no probe found for variant with ID: " + variant.getID() + "\t" + variant);
        } else {
            RawArrayFieldEnum[] fields = RawArrayFieldEnum.getCompressedRawArrayFieldEnums();
            if (!useCompressedData) {
                fields = RawArrayFieldEnum.getUncompressedRawArrayFieldEnums();
            }
            for (final RawArrayFieldEnum fieldEnum : fields) {
                row.add(fieldEnum.getColumnValue(variant, probeInfo, sampleId));
            }
        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(RawArrayFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (!variant.getFilters().contains("ZEROED_OUT_ASSAY")) {
            final List<String> rowData = createRow(variant, sampleId);

            // write the row to the XSV
            if (rowData.size() == RawArrayFieldEnum.values().length) {
                SimpleXSVWriter.LineBuilder rawLine = rawArrayWriter.getNewLineBuilder();
                rawLine.setRow(rowData);
                rawLine.write();
            }
        }
    }

    public void closeTool() {
        if (rawArrayWriter != null) {
            try {
                rawArrayWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close RAW array writer", e);
            }
        }
    }
}
