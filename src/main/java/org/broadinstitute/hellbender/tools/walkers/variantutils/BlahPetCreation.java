package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class BlahPetCreation {

    /**
     * Expected headers for the Position Table (PET)
     */
    public enum HeaderFieldEnum {
        POSITION,
        SAMPLE,
        STATE,
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
        MISSING("m");

        String value;

        GQStateEnum(String value) {
            this.value = value;
        }

    }

    public static List<List<String>> createPositionRows(final int start, final VariantContext variant, final int end) {

        List<List<String>> rows = new ArrayList<>();
        String sampleName = variant.getSampleNamesOrderedByName().get(0);

        if (!variant.isReferenceBlock()) {
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(start));
            row.add(sampleName);
            row.add(GQStateEnum.VARIANT.value);
            rows.add(row);

            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (int i = start + 1 ; i <= end; i++){
                row = new ArrayList<>();
                row.add(String.valueOf(i));
                row.add(sampleName);
                row.add(GQStateEnum.STAR.value);
                rows.add(row);
            }
        } else {
            // TODO check in the tool to make sure it's only one sample
            GQStateEnum state = getGQStateEnum(variant.getGenotype(0).getGQ());

            for (int position = start; position <= end; position++){ // break up ref blocks
                List<String> row = new ArrayList<>();

                row.add(String.valueOf(position));
                row.add(sampleName);
                row.add(state.value);
                rows.add(row);
            }
        }

        return rows;
    }

    public static List<List<String>> createSpanDelRows(final int start, final VariantContext variant, final int end) {
        if (variant.isReferenceBlock()){
            throw new IllegalStateException("Cannot create span deletion rows for a reference block");
        }

        List<List<String>> rows = new ArrayList<>();
        String sampleName = variant.getSampleNamesOrderedByName().get(0);

        for (int position = start; position <= end; position++){ // break up ref blocks
            List<String> row = new ArrayList<>();

            row.add(String.valueOf(position));
            row.add(sampleName);
            row.add(GQStateEnum.STAR.value);
            rows.add(row);
        }

        return rows;
    }

    public static List<List<String>> createMissingTSV(int start, int end, String sampleName) {
        List<List<String>> rows = new ArrayList<>();

        for (int position = start; position <= end; position ++){
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
