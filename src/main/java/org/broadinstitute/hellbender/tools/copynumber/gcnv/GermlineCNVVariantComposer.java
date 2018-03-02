package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Base class for {@link GermlineCNVIntervalVariantComposer} and {@link GermlineCNVSegmentVariantComposer}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class GermlineCNVVariantComposer<DATA extends Locatable> {
    static final String VARIANT_PREFIX = "CNV";
    static final Allele REF_ALLELE = Allele.create("N", true);
    static final Allele DEL_ALLELE = Allele.create("<DEL>", false);
    static final Allele DUP_ALLELE = Allele.create("<DUP>", false);
    static final List<Allele> ALL_ALLELES = Arrays.asList(REF_ALLELE, DEL_ALLELE, DUP_ALLELE);

    protected final VariantContextWriter outputWriter;
    protected final String sampleName;

    GermlineCNVVariantComposer(final VariantContextWriter outputWriter,
                               final String sampleName) {
        this.outputWriter = Utils.nonNull(outputWriter);
        this.sampleName = Utils.nonNull(sampleName);
    }

    /**
     * Compose the header of the variant context.
     *
     * @param vcfDefaultToolHeaderLines default header lines of the VCF generation tool
     */
    abstract void composeVariantContextHeader(final Set<VCFHeaderLine> vcfDefaultToolHeaderLines);

    abstract VariantContext composeVariantContext(final DATA data);

    public final void writeAll(final List<DATA> dataList) {
        for (final DATA data : dataList) {
            outputWriter.add(composeVariantContext(data));
        }
    }
}
