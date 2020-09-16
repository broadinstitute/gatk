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


import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class RawArrayTsvCreator {
    static final Logger logger = LogManager.getLogger(RawArrayTsvCreator.class);

    public static final String NORMX = "NORMX";
    public static final String NORMY = "NORMY";
    public static final String BAF = "BAF";
    public static final String LRR = "LRR";
    public static final GT_encoding value_to_drop = GT_encoding.HOM_REF;

    private SimpleXSVWriter rawArrayWriter = null;
    private final String sampleId;
    private final Map<String, ProbeInfo> probeDataByName;
    private static String RAW_FILETYPE_PREFIX = "raw_";

    enum GT_encoding {
        HOM_REF("R"),
        HET0_1("X"),
        HOM_VAR("A"),
        HET1_2("Y"),
        HOM_ALT2("B"),
        MISSING(".");

        String value;
        GT_encoding(String v) {
            value = v;
        }
        String getValue() {
            return value;
        }
    }

    public RawArrayTsvCreator(final String sampleName, final String sampleId, final String tableNumberPrefix, final Map<String, ProbeInfo> probeDataByName, final File outputDirectory) {
        this.sampleId = sampleId;
        this.probeDataByName = probeDataByName;
        try {
            // Create a raw file to go into the raw dir for _this_ sample
            final File rawOutputName = new File(outputDirectory, RAW_FILETYPE_PREFIX + tableNumberPrefix + sampleName  + IngestConstants.FILETYPE);
            // write header to it
            List<String> rawHeader = RawArrayTsvCreator.getHeaders();
            rawArrayWriter = new SimpleXSVWriter(rawOutputName.toPath(), IngestConstants.SEPARATOR);
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
            throw new IllegalStateException("Cannot be missing required probe ID for variant " + variant.getID() + "\t" + variant);
        }

        for (final RawArrayFieldEnum fieldEnum : RawArrayFieldEnum.getUncompressedRawArrayFieldEnums()) {
            switch (fieldEnum) {
                case sample_id:
                    row.add(sampleId);
                    break;
                case probe_id:
                    row.add(String.valueOf(probeInfo.probeId));
                    break;
                default:
                    row.add(fieldEnum.getColumnValue(variant));
            }
        }

        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(RawArrayFieldEnum.getUncompressedRawArrayFieldEnums()).map(String::valueOf).collect(Collectors.toList());
    }

    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (!variant.getFilters().contains("ZEROED_OUT_ASSAY")) {
            final List<String> rowData = createRow(variant, sampleId);

            int length = RawArrayFieldEnum.getUncompressedRawArrayFieldEnums().length;
            // write the row to the XSV
            if (rowData.size() == length) {
                SimpleXSVWriter.LineBuilder rawLine = rawArrayWriter.getNewLineBuilder();
                rawLine.setRow(rowData);
                rawLine.write();
            } else {
                throw new UserException("Length of row data didn't match length of expected row data.");
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
