package org.broadinstitute.hellbender.tools.variantdb;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class IngestVetArrayCreation {

    /**
     * Expected headers for the Variant Table (VET)
     *     position, // req
     *     sample, // req
     *     reference_bases, // req
     *     alternate_bases_alt, // req
     *     call_genotype, // req
     *
     *     IGC // IGC < --- instead of the GQ
     *     LRR // Log R ratio --> intensity value instead of DP
     *     BAF // b allele fraction --> AD proxy
     *
     */
    public enum HeaderFieldEnum {
        // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
        // and just leave it empty if not required

        position, // Required-- start position for sample
        sample, // Required-- sample Id for sample


        ref { // Required
            public String getColumnValue(final VariantContext variant) {
                final String referenceBase = variant.getReference().getBaseString();
                if (referenceBase == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for reference_bases"); // TODO, should this be UserException too?
                }
                return referenceBase;
            }
        },

        alt { // remove "<NON_REF>"
            //TODO what if this field is null and if <NON_REF> is not there--throw an error
            public String getColumnValue(final VariantContext variant) {
                List<String> outList = new ArrayList<>();
                for(Allele a : variant.getAlternateAlleles()) {
                    if (!a.isNonRefAllele()) { // TODO unit test this
                        outList.add(a.getDisplayString());
                    }
                }
                return String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
            }
        },

        call_GT {
            public String getColumnValue(final VariantContext variant) {
                IndexedAlleleList<Allele> alleleList = new IndexedAlleleList<>(variant.getAlleles());
                ArrayList<Integer> allele_indices = new ArrayList<Integer>();

                for (Allele allele : variant.getGenotype(0).getAlleles()) {
                    allele_indices.add(alleleList.indexOfAllele(allele));
                }
                if (allele_indices.size() != 2){
                    throw new IllegalArgumentException("GT doesnt have two alleles");
                }
                String separator = variant.getGenotype(0).isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED;
                return StringUtils.join(allele_indices, separator);
            }
        },

        call_IGC { // Required instead of GQ variant.getGenotypes().get(0).getExtendedAttributes().get("IGC")
            public String getColumnValue(final VariantContext variant) {
                return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("IGC"));
            }
        },

        call_LRR { // variant.getGenotypes().get(0).getExtendedAttributes().get("LRR")
            public String getColumnValue(final VariantContext variant) {
                return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("LRR"));
            }
        },

        call_BAF { // variant.getGenotypes().get(0).getExtendedAttributes().get("BAF")
            public String getColumnValue(final VariantContext variant) {
                return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("BAF"));
            }
        };

        public String getColumnValue(final VariantContext variant) {
            throw new IllegalArgumentException("Not implemented");
        }

        private static String getAttribute(VariantContext vc, String key, String defaultValue){
            Object attr = vc.getAttribute(key);
            if ( attr == null ) return defaultValue;
            if ( attr instanceof String ) return (String)attr;
            if ( attr instanceof List) return StringUtils.join((List)attr, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
            return String.valueOf(attr); // throws an exception if this isn't a string
        }
    }

    public static List<String> createVariantRow(final long start, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        row.add(String.valueOf(start));
        row.add(sampleId);
        for ( final HeaderFieldEnum fieldEnum : HeaderFieldEnum.values() ) {
            if (!fieldEnum.equals(HeaderFieldEnum.position) && !fieldEnum.equals(HeaderFieldEnum.sample)) {
                row.add(fieldEnum.getColumnValue(variant));
            }        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
