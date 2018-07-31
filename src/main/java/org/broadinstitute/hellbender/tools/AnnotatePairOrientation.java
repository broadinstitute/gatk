package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.OxoGReadCounts;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "(Experimental) This adds fields normally emitted by M2 to a VCF.  There should never be a need to run this tool on a VCF that was produced by M2." +
                "\n  The output of this tool should be usable with FilterByOrientationBias." +
                "\n  The output of this tool only counts reads that fully overlap (and match) the variant or reference sequence (this is relevant for indels)." +
                "\n  IMPORTANT:  This tool does not produce the exact same F1R2/F2R1 as M2, due to the nature of how M2 calls variants (using read likelihoods, whereas this tool uses a base quality filter).",
        oneLineSummary = "(EXPERIMENTAL) Annotate a non-M2 VCF (using the associated tumor bam) with pair orientation fields (e.g. " + GATKVCFConstants.F1R2_KEY + " ).",
        programGroup = VariantEvaluationProgramGroup.class
)
@BetaFeature
public class AnnotatePairOrientation extends VariantWalker {

    @Argument(
            doc = "Output Somatic SNP/Indel VCF file with additional annotations.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    protected File outputFile;

    public final static String CUTOFF_SHORT_NAME = "cutoff";
    public final static String CUTOFF_LONG_NAME = "min-base-quality-cutoff";
    public final static int MIN_BASE_QUALITY_DEFAULT_CUTOFF = 7;
    @Argument(
            doc = "Cutoff for the min base quality value(s) to count the read.  These are for bases that overlap the variant.",
            shortName = CUTOFF_SHORT_NAME, fullName = CUTOFF_LONG_NAME, minValue = 0, maxRecommendedValue = 20,
            optional = true
    )
    private int minBaseQualityCutoff = MIN_BASE_QUALITY_DEFAULT_CUTOFF;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outputFile);
        vcfWriter.writeHeader(createVCFHeader(getHeaderForVariants(), getCommandLine()));
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        final ReadPileup readPileup = GATKProtectedVariantContextUtils.getPileup(variant, readsContext);
        final List<Genotype> updatedGenotypes = new ArrayList<>();

        final Map<String, ReadPileup> sampleToReadPileup = readPileup.splitBySample(getHeaderForReads(), null);

        for (Genotype g : variant.getGenotypes()) {
            final ReadPileup genotypeSamplePileup = sampleToReadPileup.get(g.getSampleName());
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            OxoGReadCounts.annotateSingleVariant(variant, gb, genotypeSamplePileup, minBaseQualityCutoff);
            updatedGenotypes.add(gb.make());
        }

        vcfWriter.add(new VariantContextBuilder(variant).genotypes(updatedGenotypes).make());
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    private static VCFHeader createVCFHeader(final VCFHeader inputVCFHeader, final String commandLine) {
        Utils.nonNull(inputVCFHeader);

        // Setup header for output file
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F1R2_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F2R1_KEY));
        headerLines.add(new VCFHeaderLine("command", commandLine));
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());
        return new VCFHeader(headerLines, samples.asSetOfSamples());
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}