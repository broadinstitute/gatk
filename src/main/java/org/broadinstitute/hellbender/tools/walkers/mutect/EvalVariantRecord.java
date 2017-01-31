package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

/**
 * Created by tsato on 2/8/17.
 */
public class EvalVariantRecord {
    public static final String CHROMOSOME_COLUMN_NAME = "CHROM";
    public static final String START_POSITION_COLUMN_NAME = "START";
    public static final String END_POSITION_COLUMN_NAME = "END";
    public static final String REF_ALLELE_COLUMN_NAME = "REF";
    public static final String ALT_ALLELE_COLUMN_NAME = "ALT";
    public static final String ALLELE_FRACTION_COLUMN_NAME = "TUMOR_ALT_AF";
    public static final String TRUTH_STATUS_COLUMN_NAME = "TRUTH_STATUS";
    public static final String[] VARIANT_TABLE_COLUMN_HEADERS = {CHROMOSOME_COLUMN_NAME, START_POSITION_COLUMN_NAME, END_POSITION_COLUMN_NAME,
            REF_ALLELE_COLUMN_NAME, ALT_ALLELE_COLUMN_NAME, ALLELE_FRACTION_COLUMN_NAME, TRUTH_STATUS_COLUMN_NAME};

    String truthStatus;
    VariantContext variantContext;
    // TODO: should be an array of doubles
    String tumorAlleleFraction;

    public EvalVariantRecord(final VariantContext variantContext, final String truthStatus, final String tumorSampleName){
        this.truthStatus = truthStatus;
        this.variantContext = variantContext;
        this.tumorAlleleFraction = (String) variantContext.getGenotype(tumorSampleName).getAnyAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);
    }

    public String getTruthStatus() {
        return truthStatus;
    }

    public VariantContext getVariantContext() {
        return variantContext;
    }

    public String getTumorAlleleFraction(){ return tumorAlleleFraction; }
}