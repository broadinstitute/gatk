package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.fasterxml.jackson.annotation.JsonValue;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Prints variants supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Prints variants with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class BlahVetCreation {

    /**
     * Expected headers for the Variant Table (VET)
     */
    public enum HeaderFieldEnum {
        // TODO is this where the validation step (required vs not) lives  -- fail if there is missing data for a required field
        // and just leave it empty if not required

        START_POSITION("start_position") { // Required
             public String getColumnValue(final VariantContext variant) throws IOException {
                 final String startPosition = String.valueOf(variant.getStart());
                 if (startPosition.equals(null)) {
                     throw new IllegalArgumentException("Cannot be missing required value for start_position");
                 }
                 return startPosition; // TODO see what James' writer does
            }
        },

        REFERENCE_BASES("reference_bases") { // Required
            public String getColumnValue(final VariantContext variant) throws IOException {
                final String referenceBase = variant.getReference().getBaseString();
                if (referenceBase.equals(null)) {
                    throw new IllegalArgumentException("Cannot be missing required value for reference_bases");
                }
                return referenceBase;
            }
        },

        ALTERNATE_BASES_ALT("alternate_bases.alt") {
            public String getColumnValue(final VariantContext variant) {
                return variant.getAlternateAllele(0).getBaseString();
            }
        },

        ALTERNATE_BASES_AS_RAW_MQ("alternate_bases.AS_RAW_MQ") { // Required
            public String getColumnValue(final VariantContext variant) throws IOException {
                if (variant.getAttribute("AS_RAW_MQ").equals(null)) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_RAW_MQ");
                }
                return variant.getAttributeAsString("AS_RAW_MQ", null);
            }
        },

        ALTERNATE_BASES_AS_MQ_DP("alternate_bases.AS_MQ_DP") { // Required // TODO check for what this should be
            public String getColumnValue(final VariantContext variant) throws IOException {
                if (variant.getAttribute("AS_MQ_DP").equals(null)) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_MQ_DP");
                }
                return variant.getAttributeAsString("AS_MQ_DP", null);
            }
        },

        ALTERNATE_BASES_AS_RAW_MQRankSum("alternate_bases.AS_RAW_MQRankSum") {
            public String getColumnValue(final VariantContext variant) {
                return variant.getAttributeAsString("AS_RAW_MQRankSum", null);
            }
        },

        ALTERNATE_BASES_AS_QUAL_APPROX("alternate_bases.AS_QUALapprox") { // Required
            public String getColumnValue(final VariantContext variant) throws IOException {
                if (variant.getAttribute("AS_QUALapprox").equals(null)) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_QUALapprox");
                }
                return variant.getAttributeAsString("AS_QUALapprox", null);
            }
        },

        ALTERNATE_BASES_AS_RAW_READPOSRANKSUM("alternate_bases.AS_RAW_ReadPosRankSum") {
            public String getColumnValue(final VariantContext variant) {
                return variant.getAttributeAsString("AS_RAW_ReadPosRankSum", null);
            }
        },

        ALTERNATE_BASES_AS_SB_TABLE("alternate_bases.AS_SB_TABLE") { // Required
            public String getColumnValue(final VariantContext variant) throws IOException {
                if (variant.getAttribute("AS_SB_TABLE").equals(null)) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_SB_TABLE");
                }
                return variant.getAttributeAsString("AS_SB_TABLE", null);
            }
        },

        ALTERNATE_BASES_AS_VARDP("alternate_bases.AS_VarDP") { // Required
            public String getColumnValue(final VariantContext variant) throws IOException {
                if (variant.getAttribute("AS_VarDP").equals(null)) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_VarDP");
                }
                return variant.getAttributeAsString("AS_VarDP", null);
            }
        },

        /*FILTER("filter") {
            public Set<String> getColumnValue(final VariantContext variant) {
                return variant.getFilters();
            }
        },*/

        CALL_NAME("call.name") {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getGenotypeString(); // TODO double check
            }
        },

        CALL_GENOTYPE("call.genotype") {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getGenotypeString(); // TODO double check
            }
        },

         /*CALL_PHASESET("call.phaseset") {
            public int[] getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).get .....// TODO lookup
            }
        },

        CALL_AD("call.AD") {
            public Object getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getAD();
            }
        },

        CALL_DP("call.DP") {
            public Object getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getDP();
            }
        },

        CALL_GQ("call.GQ") { // Required
            if (variant.getGenotype(0).getGQ().equals(null)) {
               throw new IllegalArgumentException("Cannot be missing required value for call.GQ");
            }
            public Object getColumnValue(final VariantContext variant) throws IOException {
                return variant.getGenotype(0).getGQ();
            }
        },

        CALL_PGT("call.PGT") {
            public Object getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getAnyAttribute("PGT");// TODO lookup -- should be in Format -- cast to int
            }
        },

        CALL_PID("call.PID") {
            public Object getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getAnyAttribute("PID");// TODO lookup -- should be in Format
            }
        },*/


        CALL_PL("call.PL") { // TODO make this a string -- Need to check in with Laura about how we want to concat the string together
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getPL().toString();
            }
        /*},

        AS_SB_TABLE("AS_SB_TABLE") {
            public String getColumnValue(final VariantContext variant) {
                return variant.get .....// TODO lookup
            }*/
        };

        private String value;

        HeaderFieldEnum(String value) {
            this.value = value;
        }

        public String getColumnValue(final VariantContext variant) throws IOException {
            throw new IllegalArgumentException("Not implemented");
        }

        @Override
        @JsonValue
        public String toString() {
            return String.valueOf(value);
        }
    }


    public static List<String> createTSV(final VariantContext variant) throws IOException {
        List<String> row = new ArrayList<>();

        for ( final HeaderFieldEnum fieldEnum : HeaderFieldEnum.values() ) {
            row.add(fieldEnum.getColumnValue(variant));
        }
        return row;
    }
}
