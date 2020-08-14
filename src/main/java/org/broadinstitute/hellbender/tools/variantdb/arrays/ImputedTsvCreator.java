package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
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
import java.util.*;
import java.util.stream.Collectors;

public final class ImputedTsvCreator {
    static final Logger logger = LogManager.getLogger(ImputedTsvCreator.class);

    private SimpleXSVWriter arrayWriter = null;
    private final String runId;
    private final boolean useCompressedData;
    private static String IMPUTED_FILETYPE_PREFIX = "imputed_";
    private final Set<String> sampleNames;
    private final Map<String, Integer> sampleNameMap;


    public ImputedTsvCreator(final String tableNumberPrefix, final String runId, final Set<String> sampleNames, final Map<String, Integer> sampleNameMap, final boolean useCompressedData) {
        this.useCompressedData = useCompressedData;
        this.runId = runId;
        this.sampleNames = sampleNames;
        this.sampleNameMap = sampleNameMap;
        try {
            // Create a raw file to go into the raw dir for _this_ sample
            final String imputedOutputName = IMPUTED_FILETYPE_PREFIX + tableNumberPrefix + "_" + runId  + IngestConstants.FILETYPE;
            // write header to it
            List<String> imputedHeader = ImputedTsvCreator.getHeaders();
            arrayWriter = new SimpleXSVWriter(Paths.get(imputedOutputName), IngestConstants.SEPARATOR);
            arrayWriter.setHeaderLine(imputedHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create raw outputs", e);
        }
    }

    public List<String> createRow(final Genotype gt, final VariantContext variant) {
        List<String> row = new ArrayList<>();
        ImputedFieldEnum[] fields = ImputedFieldEnum.getUncompressedRawArrayFieldEnums();
//        if (useCompressedData) {
//            fields = ImputedFieldEnum.getCompressedRawArrayFieldEnums();
//        }
        for (final ImputedFieldEnum fieldEnum : fields) {
            row.add(fieldEnum.getColumnValue(variant, sampleNameMap.get(gt.getSampleName()).toString()));
        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(ImputedFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void apply(final VariantContext variant) { //, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        GenotypesContext gts = variant.getGenotypes(sampleNames);

        gts.forEach(genotype -> {
        final List<String> rowData = createRow(genotype, variant);

        // write the row to the XSV
        if (rowData.size() == RawArrayFieldEnum.values().length) {
            SimpleXSVWriter.LineBuilder rawLine = arrayWriter.getNewLineBuilder();
            rawLine.setRow(rowData);
            rawLine.write();
        }
        });
    }

    public void closeTool() {
        if (arrayWriter != null) {
            try {
                arrayWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close imputed array writer", e);
            }
        }
    }
}
