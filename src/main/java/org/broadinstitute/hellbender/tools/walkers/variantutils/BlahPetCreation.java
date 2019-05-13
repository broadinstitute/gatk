package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.fasterxml.jackson.annotation.JsonValue;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;

import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties( // TODO -- how should this be edited?
        summary = "Creates position expanded table",
        oneLineSummary = "Creates position expanded table",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class BlahPetCreation {

    /**
     * Expected headers for the Variant Table (VET)
     */
    public enum GQStateEnum {
        VARIANT("v"),

        DELETION("s"), // TODO is this s or * ?

        ZERO("0"),

        TEN("10"),

        TWENTY("20"),

        THIRTY("30"),

        FORTY("40"),

        FIFTY("50"),

        SIXTY("60");

        private String value;

        GQStateEnum(String value) {
            this.value = value;
        }

        @Override
        @JsonValue
        public String toString() {
            return String.valueOf(value);
        }
    }

    public static List<String> createTSV(final VariantContext variant) {

        // TODO hardcode the dropped band for the GQ

        List<String> row = new ArrayList<>();
        String sampleName = variant.getSampleNamesOrderedByName().get(0);

        if (variant.isVariant()){
             row.add(variant.getStart() + "," + sampleName + "," + GQStateEnum.VARIANT.toString());
            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (int i = variant.getStart() + 1 ; i < variant.getEnd(); i++){
                row.add(i + "," + sampleName + "," + GQStateEnum.DELETION.toString());
            }
        } else {
            int genotypeQual = variant.getGenotype(0).getGQ();
            if (genotypeQual < 10) {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){ // break up ref blocks
                    row.add(i + "," + sampleName + "," + GQStateEnum.ZERO.toString());
                }
            } else if (genotypeQual < 20) {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){
                    row.add(i + "," + sampleName + "," + GQStateEnum.TEN.toString());
                }
            } else if (genotypeQual < 30) {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){
                    row.add(i + "," + sampleName + "," + GQStateEnum.TWENTY.toString());
                }
            } else if (genotypeQual < 40) {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){
                    row.add(i + "," + sampleName + "," + GQStateEnum.THIRTY.toString());
                }
            } else if (genotypeQual < 50) {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){
                    row.add(i + "," + sampleName + "," + GQStateEnum.FORTY.toString());
                }
            } else if (genotypeQual < 60) {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){
                    row.add(i + "," + sampleName + "," + GQStateEnum.FIFTY.toString());
                }
            } else {
                for (int i = variant.getStart(); i < variant.getEnd(); i++){
                    row.add(i + "," + sampleName + "," + GQStateEnum.SIXTY.toString());
                }
            }
        } // TODO is there an exception thrown if there's no GQ?

        return row;
    }
}
