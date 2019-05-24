package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * Implements fields for use in known locatables.
 *
 * This class should only be used in specialized cases.
 *
 * This class creates a funcotation with only the locatable inforamtion (CONTIG, START, and END) from a locatable.
 *
 * If you are trying to create funcotations based on the attributes of a variant context, please see
 * {@link FuncotatorUtils#createFuncotations(VariantContext, FuncotationMetadata, String)}
 *
 */
public class LocatableFuncotationCreator {

    @VisibleForTesting
    final static String CONTIG_FIELD_NAME = "CONTIG";

    @VisibleForTesting
    final static String START_FIELD_NAME = "START";

    @VisibleForTesting
    final static String END_FIELD_NAME = "END";

    @VisibleForTesting
    final static FuncotationMetadata METADATA =
            VcfFuncotationMetadata.create(
                    Arrays.asList(
                            new VCFInfoHeaderLine(CONTIG_FIELD_NAME,1, VCFHeaderLineType.String, "contig"),
                            new VCFInfoHeaderLine(START_FIELD_NAME,1, VCFHeaderLineType.Integer, "start position"),
                            new VCFInfoHeaderLine(END_FIELD_NAME,1, VCFHeaderLineType.Integer, "end position")

                    )
            );

    /**
     * Create a funcotation with only the locatable inforamtion (CONTIG, START, and END) from a locatable.
     *
     * If you are trying to create funcotations based on the attributes of a variant context, please see
     * {@link FuncotatorUtils#createFuncotations(VariantContext, FuncotationMetadata, String)}
     *
     * @param locatable Never {@code null}
     * @param altAllele Never {@code null}
     * @param dataSourceName Never {@code null}
     * @return a Funcotation with field names {@link LocatableFuncotationCreator#CONTIG_FIELD_NAME},
     *  {@link LocatableFuncotationCreator#START_FIELD_NAME}, and {@link LocatableFuncotationCreator#END_FIELD_NAME}.
     *  Never {@code null}.  Metadata is populated, see {@link LocatableFuncotationCreator#METADATA}
     */
    public static Funcotation create(final Locatable locatable, final Allele altAllele, final String dataSourceName) {
        Utils.nonNull(locatable);
        Utils.nonNull(altAllele);
        Utils.nonNull(dataSourceName);
        return TableFuncotation.create(Arrays.asList(CONTIG_FIELD_NAME, START_FIELD_NAME, END_FIELD_NAME),
                Arrays.asList(locatable.getContig(), String.valueOf(locatable.getStart()), String.valueOf(locatable.getEnd())),
                altAllele, dataSourceName, METADATA);
    }
}
