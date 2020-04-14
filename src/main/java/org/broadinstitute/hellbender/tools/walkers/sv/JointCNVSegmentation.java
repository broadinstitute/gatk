package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVClusterEngine;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordWithEvidence;
import org.broadinstitute.hellbender.tools.sv.SVDepthOnlyCallDefragmenter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@BetaFeature
@CommandLineProgramProperties(
        summary = "Gathers single-sample segmented gCNV VCFs, harmonizes breakpoints, and outputs a cohort VCF with genotypes.",
        oneLineSummary = "Combined single-sample segmented gCNV VCFs.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class JointCNVSegmentation extends MultiVariantWalkerGroupedOnStart {

    private SortedSet<String> samples;
    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary dictionary;
    private SVDepthOnlyCallDefragmenter defragmenter;
    private SVClusterEngine clusterEngine;

    private String currentContig;

    public static final String MIN_QUALITY_LONG_NAME = "minimum-qs-score";

    @Argument(fullName = MIN_QUALITY_LONG_NAME, doc = "Minimum QS score to combine a variant segment")
    private int minQS = 20;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined output file", optional=false)
    private File outputFile;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary, 0.0);
        clusterEngine = new SVClusterEngine(dictionary);

        vcfWriter = getVCFWriter();
    }

    private VariantContextWriter getVCFWriter() {
        samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));

        VariantContextWriter writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /**
     * @param variantContexts  VariantContexts from driving variants with matching start positon
     *                         NOTE: This will never be empty
     * @param referenceContext ReferenceContext object covering the reference of the longest spanning VariantContext
     * @param readsContexts
     */
    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        if (currentContig == null) {
            currentContig = variantContexts.get(0).getContig(); //variantContexts should have identical start, so choose 0th arbitrarily
        } else if (!variantContexts.get(0).getContig().equals(currentContig)) {
            processClusters();
        }
        for (final VariantContext vc : variantContexts) {
            defragmenter.add(new SVCallRecordWithEvidence(SVCallRecord.createDepthOnlyFromGCNV(vc, minQS)));
        }
    }

    private void processClusters() {
        final List<SVCallRecordWithEvidence> defragmentedCalls = defragmenter.getOutput();
        defragmentedCalls.stream().forEachOrdered(clusterEngine::add);
        //defragmented calls may still be overlapping, so run the clustering engine to combine based on reciprocal overlap
        final List<SVCallRecordWithEvidence> clusteredCalls = clusterEngine.getOutput();
        write(clusteredCalls);
    }

    private void write(final List<SVCallRecordWithEvidence> calls) {
        calls.stream()
                .sorted(Comparator.comparing(c -> c.getStartAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(this::buildVariantContext)
                .forEachOrdered(vcfWriter::add);
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final Allele refAllele = Allele.REF_N;
        final VariantContextBuilder builder = new VariantContextBuilder("", call.getContig(), call.getStart(), call.getEnd(),
                Lists.newArrayList(refAllele, altAllele));
        builder.attribute(VCFConstants.END_KEY, call.getEnd());
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        builder.attribute(VCFConstants.SVTYPE, call.getType());
        final List<Genotype> genotypes = new ArrayList<>();
        for (final String sample : samples) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            if (call.getSamples().contains(sample)) {
                genotypeBuilder.alleles(Lists.newArrayList(refAllele, altAllele));
                final Genotype currentGenotype = call.getGenotypes().stream().filter(g -> g.getSampleName().equals(sample)).collect(Collectors.toList()).get(0);
                if (currentGenotype.hasAnyAttribute(GermlineCNVSegmentVariantComposer.CN)) {
                    genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, currentGenotype.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN));
                }

            } else {
                genotypeBuilder.alleles(Lists.newArrayList(refAllele, refAllele));
                genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, HomoSapiensConstants.DEFAULT_PLOIDY);
            }
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder.make();
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
