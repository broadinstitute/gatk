package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.ArrayDeque;

/**
 * Remove indels that are close to another indel from a vcf.  This is a preprocessing step
 * for the CRSP sensitivity validation truth data.
 *
 * Example usage:
 * java -jat gatk.jar RemoveNearbyIndels -V input.vcf -O output.vcf -minIndelSpacing 20
 *
 * Created by David Benjamin on 1/30/17.
 */
@CommandLineProgramProperties(
        summary = "Remove indels that are close to each other from a vcf.  For any pair of indels that are within" +
                "some minimum allowed distance, both indels are removed, regardless of any intervening non-indel variants.",
        oneLineSummary = "Remove indels that are close to each other from a vcf",
        programGroup = VariantProgramGroup.class
)
public class RemoveNearbyIndels extends VariantWalker {

    public static final String MIN_INDEL_SPACING_NAME = "minIndelSpacing";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output filtered VCF file",
            optional = false)
    private final String outputVcf = null;

    @Argument(fullName = MIN_INDEL_SPACING_NAME,
            shortName = MIN_INDEL_SPACING_NAME,
            doc = "Minimum spacing between neighboring indels to be emitted",
            optional = false)
    private int minIndelSpacing = 1;

    private VariantContextWriter vcfWriter;

    private VariantBuffer variantBuffer = new VariantBuffer();

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final VCFHeader vcfHeader = new VCFHeader(inputHeader.getMetaDataInSortedOrder(), inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        variantBuffer.add(vc);
    }

    @Override
    public Object onTraversalSuccess() {
        variantBuffer.emitRemaining();
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    private final class VariantBuffer {

        private VariantContext lastIndel = null;

        // INVARIANT: this buffer will contain at most one indel
        private final ArrayDeque<VariantContext> buffer = new ArrayDeque<>();

        public VariantBuffer() { }

        public void add(final VariantContext vc) {
            if (lastIndel == null && vc.isIndel()) {
                buffer.add(vc);
                lastIndel = vc;
            } else if (nearby(lastIndel, vc)) {
                if (vc.isIndel()) {
                    emitAllNonIndels(); // throw out {@code lastIndel} and {@code vc}
                } else {
                    buffer.add(vc);
                }
            } else {
                emitAllVariants();
                buffer.add(vc);
            }

            lastIndel = vc.isIndel() ? vc : lastIndel;
        }

        private boolean nearby(final VariantContext left, final VariantContext right) {
            return left != null && left.getContig().equals(right.getContig())
                    && (right.getStart() - left.getEnd() < minIndelSpacing);
        }

        private void emitAllNonIndels() {
            buffer.stream().filter(vc -> !vc.isIndel()).forEach(vcfWriter::add);
            buffer.clear();
        }

        private void emitAllVariants() {
            buffer.forEach(vcfWriter::add);
            buffer.clear();
        }

        public void emitRemaining() {
            buffer.stream().filter(vc -> !(vc.isIndel() && nearby(lastIndel, vc)) || vc == lastIndel).forEach(vcfWriter::add);
            buffer.clear();
        }
    }
}
