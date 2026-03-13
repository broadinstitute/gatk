package org.broadinstitute.hellbender.tools.gvs.ingest.tsv;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.ingest.VetFieldEnum;
import org.broadinstitute.hellbender.tools.gvs.ingest.VetWriter;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class VetTsvWriter extends SimpleXSVWriter implements VetWriter {

    private final boolean skipLoadingVqsrFields;
    private final boolean forceLoadingFromNonAlleleSpecific;

    public VetTsvWriter(File file, boolean skipLoadingVqsrFields, boolean forceLoadingFromNonAlleleSpecific) throws IOException {
        super(file.toPath(), IngestConstants.SEPARATOR);
        this.skipLoadingVqsrFields = skipLoadingVqsrFields;
        this.forceLoadingFromNonAlleleSpecific = forceLoadingFromNonAlleleSpecific;
        this.setHeaderLine(getHeaders());
    }

    private List<String> getHeaders() {
        return Arrays.stream(VetFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    @Override
    public void write(long location, VariantContext variant, long sampleId) throws IOException {
        final List<String> row = createRow(
                location,
                variant,
                String.valueOf(sampleId)
        );
        this.getNewLineBuilder().setRow(row).write();
    }


    public List<String> createRow(final long location, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        for ( final VetFieldEnum fieldEnum : VetFieldEnum.values() ) {
            if (fieldEnum.equals(VetFieldEnum.location)) {
                row.add(String.valueOf(location));
            } else if (fieldEnum.equals(VetFieldEnum.sample_id)) {
                row.add(sampleId);
            } else if (!(skipLoadingVqsrFields && fieldEnum.isVqsrSpecificField())) {
                row.add(fieldEnum.getColumnValue(variant, forceLoadingFromNonAlleleSpecific));
            }
        }
        return row;
    }
}
