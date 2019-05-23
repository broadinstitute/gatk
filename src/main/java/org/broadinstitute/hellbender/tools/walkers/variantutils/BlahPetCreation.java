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

    public static List<List<String>> createPositionRows(final VariantContext variant) {

        List<List<String>> rows = new ArrayList<>();
        String sampleName = variant.getSampleNamesOrderedByName().get(0);

        if (!variant.isReferenceBlock()) {
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(variant.getStart()));
            row.add(sampleName);
            row.add(GQStateEnum.VARIANT.value);
            rows.add(row);

            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (int i = variant.getStart() + 1 ; i <= variant.getEnd(); i++){
                row = new ArrayList<>();
                row.add(String.valueOf(i));
                row.add(sampleName);
                row.add(GQStateEnum.STAR.value);
                rows.add(row);
            }
        } else {
            // TODO check in the tool to make sure it's only one sample
            int genotypeQual = variant.getGenotype(0).getGQ();  // ok because we only have one sample
            GQStateEnum state;

            if (genotypeQual < 10) {
                state = GQStateEnum.ZERO;
            } else if (genotypeQual < 20) {
                state = GQStateEnum.TEN;
            } else if (genotypeQual < 30) {
                state = GQStateEnum.TWENTY;
            } else if (genotypeQual < 40) {
                state = GQStateEnum.THIRTY;
            } else if (genotypeQual < 50) {
                state = GQStateEnum.FORTY;
            } else if (genotypeQual < 60) {
                state = GQStateEnum.FIFTY;
            } else if (genotypeQual >= 60) {
                state = GQStateEnum.SIXTY;
            } else {
                throw new IllegalArgumentException("GQ is not in the range we expect");
            }

            for (int position = variant.getStart(); position <= variant.getEnd(); position++){ // break up ref blocks
                List<String> row = new ArrayList<>();
                row.add(String.valueOf(position));
                row.add(sampleName);
                row.add(state.value);
                rows.add(row);
            }
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

    public static List<String> getHeaders() {
        return Arrays.stream(HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
