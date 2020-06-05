package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.MethylationProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
/**
 *
 * <p>Identifies methylated bases from bisulfite sequencing data. Given a bisulfite sequenced, methylation-aware
 * aligned BAM and a reference, it outputs methylation-site coverage to a specified output vcf file.</p>
 *
 * <h>Usage example</h>
 * <pre>
 *     gatk MethylationTypeCaller \
 *     -R "GRCm38_primary_assembly_genome.fa \
 *     -I bisulfite_input.bam \
 *     -O output.vcf
 * </pre>
 *
 * @author Benjamin Carlin
 */
@CommandLineProgramProperties(
        summary = "Tool that prints methylation-based coverage from supplied bisulfite sequenced, methylation-aware aligned BAM to the specified output vcf file",
        oneLineSummary = "Identify methylated bases from bisulfite sequenced, methylation-aware BAMs",
        programGroup = MethylationProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class MethylationTypeCaller extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output VCF file")
    private GATKPath outputFile;

    private VariantContextWriter vcfWriter;

    private static final int REFERENCE_CONTEXT_LENGTH = 2;


    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outputFile.toPath());
        vcfWriter.writeHeader(createMethylationHeader(getHeaderForReads(), getDefaultToolVCFHeaderLines()));
    }

    private static VCFHeader createMethylationHeader(SAMFileHeader header, Set<VCFHeaderLine> headerLines) {
        if(header == null) {
            throw new UserException.BadInput("Error writing header, getHeaderForReads() returns null");
        }

        final VCFInfoHeaderLine unconvertedCoverageLine = new VCFInfoHeaderLine(GATKVCFConstants.UNCONVERTED_BASE_COVERAGE_KEY, 1, VCFHeaderLineType.Integer, "Count of reads supporting methylation that are unconverted ");
        final VCFInfoHeaderLine coverageLine = new VCFInfoHeaderLine(GATKVCFConstants.CONVERTED_BASE_COVERAGE_KEY, 1, VCFHeaderLineType.Integer, "Count of reads supporting methylation that are converted ");
        final VCFInfoHeaderLine contextLine = new VCFInfoHeaderLine(GATKVCFConstants.METHYLATION_REFERENCE_CONTEXT_KEY, 1, VCFHeaderLineType.String, "Forward Strand Reference context");
        final VCFInfoHeaderLine readDepthLine = VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY);
        final VCFFormatHeaderLine gtLine = VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY);

        headerLines.add(unconvertedCoverageLine);
        headerLines.add(coverageLine);
        headerLines.add(contextLine);
        headerLines.add(readDepthLine);
        headerLines.add(gtLine);

        final List<String> samples = header.getReadGroups()
                .stream()
                .map(SAMReadGroupRecord::getSample)
                .sorted()
                .distinct()
                .collect(Collectors.toList());

        return new VCFHeader(headerLines, samples);
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final byte referenceBase = referenceContext.getBases()[0];
        final int unconvertedBases;
        final int convertedBases;
        final byte alt;
        byte[] context = null;

        // check the forward strand for methylated coverage
        if (referenceBase == (byte)'C') {
            alt = (byte)'T';
            final ReadPileup forwardBasePileup = alignmentContext.stratify(AlignmentContext.ReadOrientation.FORWARD).getBasePileup();
            // unconverted: C, index=1; converted: T, index=3
            final int[] forwardBaseCounts = forwardBasePileup.getBaseCounts();
            unconvertedBases = forwardBaseCounts[BaseUtils.simpleBaseToBaseIndex((byte)'C')];
            convertedBases = forwardBaseCounts[BaseUtils.simpleBaseToBaseIndex((byte)'T')];

            // if there is methylated coverage
            if (unconvertedBases + convertedBases > 0) {
                context = referenceContext.getBases(0, REFERENCE_CONTEXT_LENGTH);
            }
        }
        // check the reverse strand for methylated coverage
        else if (referenceBase == (byte)'G') {
            alt = (byte)'A';
            final ReadPileup reverseBasePileup = alignmentContext.stratify(AlignmentContext.ReadOrientation.REVERSE).getBasePileup();
            // unconverted: G, index=2; converted: A, index=0
            final int[] reverseBaseCounts = reverseBasePileup.getBaseCounts();
            unconvertedBases = reverseBaseCounts[BaseUtils.simpleBaseToBaseIndex((byte)'G')];
            convertedBases = reverseBaseCounts[BaseUtils.simpleBaseToBaseIndex((byte)'A')];

            // if there is methylated coverage
            if (unconvertedBases + convertedBases > 0) {
                // get the reverse complement for context b/c we are on the reverse strand
                context = BaseUtils.simpleReverseComplement(referenceContext.getBases(REFERENCE_CONTEXT_LENGTH,0));
            }
        }
        // if reference strand does not support methylation
        else {
            return;
        }

        // if there are reads that have methylated coverage
        if (context != null) {
            final LinkedHashSet<Allele> alleles = new LinkedHashSet<>();
            alleles.add(Allele.create(referenceBase, true));
            alleles.add(Allele.create(alt, false));

            final VariantContextBuilder vcb = new VariantContextBuilder();
            vcb.chr(alignmentContext.getContig());
            vcb.start(alignmentContext.getPosition());
            vcb.stop(alignmentContext.getPosition());
            vcb.noID();
            vcb.alleles(alleles);
            vcb.noGenotypes();
            vcb.unfiltered();
            vcb.attribute(GATKVCFConstants.UNCONVERTED_BASE_COVERAGE_KEY, unconvertedBases);
            vcb.attribute(GATKVCFConstants.CONVERTED_BASE_COVERAGE_KEY, convertedBases);
            vcb.attribute(GATKVCFConstants.METHYLATION_REFERENCE_CONTEXT_KEY, new String(context));
            vcb.attribute(VCFConstants.DEPTH_KEY, alignmentContext.size());

            // write to VCF
            final VariantContext vc = vcb.make();
            vcfWriter.add(vc);
        }
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }
}
