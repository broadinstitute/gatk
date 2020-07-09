package org.broadinstitute.hellbender.tools.variantdb.ingest.nextgen;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class GenomeVetTsvCreator extends VetTsvCreator {
    static final Logger logger = LogManager.getLogger(GenomeVetTsvCreator.class);

    public GenomeVetTsvCreator(String sampleName, String sampleId, String tableNumberPrefix) {
        super(sampleName, sampleId, tableNumberPrefix);
    }

    @Override
    public List<String> createRow(final long start, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        row.add(String.valueOf(SchemaUtils.encodeLocation(variant.getContig(), variant.getStart())));
        row.add(sampleId);
        for ( final GenomeFieldEnum fieldEnum : GenomeFieldEnum.values() ) {
            if (!fieldEnum.equals(GenomeFieldEnum.location) && !fieldEnum.equals(GenomeFieldEnum.sample)) {
                row.add(fieldEnum.getColumnValue(variant));
            }        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(GenomeFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
