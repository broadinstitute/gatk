package org.broadinstitute.hellbender.tools.variantdb.ingest;

import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class IngestPetCreation {

    /**
     * Expected headers for the Position Table (PET)
     */
    public enum HeaderFieldEnum {
        position,
        sample,
        state,
    }

    public enum GQStateEnum {
        VARIANT("v"),
        STAR("*"),
        ZERO("0"),
        TEN("1"),
        TWENTY("2"),
        THIRTY("3"),
        FORTY("4"),
        FIFTY("5"),
        SIXTY("6"),
        MISSING("m"),
        UNKNOWN("u");

        String value;

        GQStateEnum(String value) {
            this.value = value;
        }

    }

    public static List<List<String>> createPositionRows(final long start, final long end, final VariantContext variant, final String sampleId) {

        List<List<String>> rows = new ArrayList<>();

        if (!variant.isReferenceBlock()) {
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(start));
            row.add(sampleId);
            row.add(GQStateEnum.VARIANT.value);
            rows.add(row);

            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (long i = start + 1 ; i <= end; i++){
                row = new ArrayList<>();
                row.add(String.valueOf(i));
                row.add(sampleId);
                row.add(GQStateEnum.STAR.value);
                rows.add(row);
            }
        } else {
            // TODO check in the tool to make sure it's only one sample
            GQStateEnum state = getGQStateEnum(variant.getGenotype(0).getGQ());

            for (long position = start; position <= end; position++){ // break up ref blocks
                List<String> row = new ArrayList<>();

                row.add(String.valueOf(position));
                row.add(sampleId);
                row.add(state.value);
                rows.add(row);
            }
        }

        return rows;
    }

    public static List<List<String>> createArrayPositionRows(final long start, final long end, final VariantContext variant, final String sampleName) {

        List<List<String>> rows = new ArrayList<>();

        List<String> row = new ArrayList<>();
        row.add(String.valueOf(start));
        row.add(sampleName);
        row.add(GQStateEnum.UNKNOWN.value);
        rows.add(row);

        return rows;
    }

    public static List<List<String>> createSpanDelRows(final long start, final long end, final VariantContext variant, final String sampleName) {
        if (variant.isReferenceBlock()){
            throw new IllegalStateException("Cannot create span deletion rows for a reference block");
        }

        List<List<String>> rows = new ArrayList<>();

        for (long position = start; position <= end; position++){ // break up ref blocks
            List<String> row = new ArrayList<>();

            row.add(String.valueOf(position));
            row.add(sampleName);
            row.add(GQStateEnum.STAR.value);
            rows.add(row);
        }

        return rows;
    }

    public static List<List<String>> createMissingTSV(long start, long end, String sampleName) {
        List<List<String>> rows = new ArrayList<>();

        for (long position = start; position <= end; position ++){
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(position));
            row.add(sampleName);
            row.add(GQStateEnum.MISSING.value);
            rows.add(row);
        }

        return rows;
    }

    public static GQStateEnum getGQStateEnum(int GQ){
        if (GQ < 10) {
            return GQStateEnum.ZERO;
        } else if (GQ < 20) {
            return GQStateEnum.TEN;
        } else if (GQ < 30) {
            return GQStateEnum.TWENTY;
        } else if (GQ < 40) {
            return GQStateEnum.THIRTY;
        } else if (GQ < 50) {
            return GQStateEnum.FORTY;
        } else if (GQ < 60) {
            return GQStateEnum.FIFTY;
        } else {
            return GQStateEnum.SIXTY;
        }

    }

    public static List<String> getHeaders() {
        return Arrays.stream(HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
