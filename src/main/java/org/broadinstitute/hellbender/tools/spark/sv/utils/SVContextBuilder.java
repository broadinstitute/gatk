package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.OptionalLong;

public class SVContextBuilder extends VariantContextBuilder {

    private OptionalLong startTimeInMillis = OptionalLong.empty();

    private GenotypesContext genotypesContext;

    public SVContextBuilder() {
    }

    public SVContextBuilder(final SVContext initial) {
        super(initial);
    }

    public void setGenotypingContextSizes(final int numberOfTemplates, final int numberOfHaplotypes) {
        attribute(GATKSVVCFConstants.TEMPLATE_COUNT, numberOfTemplates);
        attribute(GATKSVVCFConstants.HAPLOTYPE_COUNT, numberOfHaplotypes);
    }


    public void startRecordingProcessingTime() {
        startTimeInMillis = OptionalLong.of(System.currentTimeMillis());
    }

    public void stopRecordingProcessingTime() {
        if (startTimeInMillis.isPresent()) {
            final long stopTimeInMillis = System.currentTimeMillis();
            @SuppressWarnings("ConstantConditions") // the enclosing if makes sure that this will work.
            final long elapse = stopTimeInMillis - startTimeInMillis.getAsLong();
            attribute(GATKSVVCFConstants.GENOTYPING_PROCESSING_TIME, elapse);
            startTimeInMillis = OptionalLong.empty();
        } else {
            throw new IllegalStateException("start recording processing time wasn't called earlier");
        }
    }

    public SVContext make() {
        if (startTimeInMillis.isPresent()) {
            stopRecordingProcessingTime();
        }
        return SVContext.of(super.make());
    }
}
