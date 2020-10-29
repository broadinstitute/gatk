package org.broadinstitute.hellbender.tools.variantdb.nextgen;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class ExomeVetTsvCreator extends VetTsvCreator {

    static final Logger logger = LogManager.getLogger(ExomeVetTsvCreator.class);

    public ExomeVetTsvCreator(String sampleName, String sampleId, String tableNumberPrefix, final File outputDirectory) {
        super(sampleName, sampleId, tableNumberPrefix, outputDirectory, getHeaders());
    }

    @Override
    public List<String> createRow(final long location, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        for ( final ExomeFieldEnum fieldEnum : ExomeFieldEnum.values() ) {
            if (fieldEnum.equals(ExomeFieldEnum.location)) {
                row.add(String.valueOf(location));
            } else if (fieldEnum.equals(ExomeFieldEnum.sample_id)) {
                row.add(sampleId);
            } else {
                row.add(fieldEnum.getColumnValue(variant));
            }
        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(ExomeFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
