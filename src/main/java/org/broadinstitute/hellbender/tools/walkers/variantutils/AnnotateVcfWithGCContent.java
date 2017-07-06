package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.walkers.annotator.GCContent;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.IntStream;

/**
 * Created by tsato on 7/6/17.
 */
@CommandLineProgramProperties(
        summary = "Annotate a vcf with the GC content of the reference around each variant locus",
        oneLineSummary = "(Internal) Annotate a vcf with the GC content of the reference around each variant locus",
        programGroup = VariantProgramGroup.class
)

@DocumentedFeature
public class AnnotateVcfWithGCContent extends VariantWalker {
    @Argument(fullName=StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Annotated VCF file",
            optional=false)
    private final String outputVcf = null;

    @Argument(fullName="KmerSize",
            shortName="K",
            optional=true)
    private int k = 100;

    private VariantContextWriter vcfWriter;

    private String GC_CONTENT_ANNOTATION_NAME = "GC";

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(new VCFInfoHeaderLine(GC_CONTENT_ANNOTATION_NAME, 1, VCFHeaderLineType.Float, "GC content in the reference around the variant"));
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        SimpleInterval oldWindow = refContext.getWindow();
        // TODO: make interval precise
        refContext.setWindow(k/2, k/2);
        final int C = 67;
        final int G = 71;
        final byte[] bases = refContext.getBases();
        final long gccount = IntStream.range(0, bases.length).map(i -> bases[i]).filter(j -> j == C || j == G).count();
        final double gcfraction = (double) gccount / bases.length;
        // TODO: restore window?

        vcfWriter.add(new VariantContextBuilder(vc).attribute(GC_CONTENT_ANNOTATION_NAME, gcfraction).make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}