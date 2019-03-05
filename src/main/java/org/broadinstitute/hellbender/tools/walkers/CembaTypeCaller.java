package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Example tool that prints methtylation-based coverage from supplied read to the specified output vcf file
 * walker that prints methtylation-based coverage with contextual data
 *
 * @author Benjamin Carlin
 */
@CommandLineProgramProperties(
    summary = "Example tool that prints methtylation-based coverage from supplied read to the specified output vcf file",
    oneLineSummary = "walker that prints methtylation-based coverage with contextual data",
    programGroup = ExampleProgramGroup.class
)
public class CembaTypeCaller extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output VCF file")
    private File OUTPUT_FILE = null;

    private htsjdk.variant.variantcontext.writer.VariantContextWriter vcfWriter = null;

    @Override
    public void onTraversalStart() {
        try {
            vcfWriter = createVCFWriter(OUTPUT_FILE);
        }
        catch ( RuntimeException e ) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }
        vcfWriter.writeHeader(createMethylationHeader());
    }

    private VCFHeader createMethylationHeader() {
        VCFInfoHeaderLine mcLine = new VCFInfoHeaderLine("mC", 1, VCFHeaderLineType.Integer, "count of reads supporting methylation that are unconverted ");
        VCFInfoHeaderLine coverageLine = new VCFInfoHeaderLine("methCov", 1, VCFHeaderLineType.Integer, "count of reads supporting methylation");
        VCFInfoHeaderLine contextLine = new VCFInfoHeaderLine("CT", 1, VCFHeaderLineType.Integer, "reference context");
        VCFInfoHeaderLine readDepthLine = VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY);
        VCFFormatHeaderLine gtLine = VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY);

        LinkedHashSet<VCFHeaderLine> headerLines = new LinkedHashSet<>();
        headerLines.add(readDepthLine);
        headerLines.add(contextLine);
        headerLines.add(mcLine);
        headerLines.add(coverageLine);
        headerLines.add(gtLine);

        htsjdk.samtools.SAMFileHeader header = getHeaderForReads();
        if(header != null) {
            List<String> samples = header.getReadGroups()
                    .stream()
                    .map(SAMReadGroupRecord::getSample)
                    .sorted()
                    .distinct()
                    .collect(Collectors.toList());

            return new VCFHeader(headerLines, samples);
        }

        return null;
    }


    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if(referenceContext.hasBackingDataSource()) {
            String strand = new String(referenceContext.getBases());
            int unconverted_bases, converted_bases = unconverted_bases = 0;
            String alt, context = alt = null;

            // check the forward strand for methylated coverage
            if(strand.equals("C")) {
                ReadPileup forwardBasePileup = alignmentContext.stratify(AlignmentContext.ReadOrientation.FORWARD).getBasePileup();
                // unconverted: C; converted: T
                unconverted_bases = forwardBasePileup.getBaseCounts()[1];
                converted_bases = forwardBasePileup.getBaseCounts()[3];

                // if there is methylated coverage
                if (unconverted_bases + converted_bases > 0) {
                    alt = "T";
                    context = new String(referenceContext.getBases(0,2));
                }
            }
            // check the reverse strand for methylated coverage
            else if(strand.equals("G")) {
                ReadPileup reverseBasePileup = alignmentContext.stratify(AlignmentContext.ReadOrientation.REVERSE).getBasePileup();
                // unconverted: G; converted: A
                unconverted_bases = reverseBasePileup.getBaseCounts()[2];
                converted_bases = reverseBasePileup.getBaseCounts()[0];

                // if there is methylated coverage
                if (unconverted_bases + converted_bases > 0) {
                    alt = "A";
                    // get the reverse complement for context b/c we are on the reverse strand
                    context = new String(org.broadinstitute.hellbender.utils.BaseUtils.simpleReverseComplement(
                            referenceContext.getBases(2,0)));
                }
            }

            // if there are reads that have methylated coverage
            if (alt != null) {
                LinkedHashSet<Allele> alleles = new LinkedHashSet<>();
                alleles.add(Allele.create(strand, true));
                alleles.add(Allele.create(alt, false));

                VariantContextBuilder vcb = new VariantContextBuilder();
                vcb.chr(alignmentContext.getContig());
                vcb.start(alignmentContext.getPosition());
                vcb.stop(alignmentContext.getPosition());
                vcb.noID();
                vcb.alleles(alleles);
                vcb.noGenotypes();
                vcb.unfiltered();
                vcb.attribute("mC", unconverted_bases);
                vcb.attribute("methCov", converted_bases);
                vcb.attribute("CT", context);
                vcb.attribute(VCFConstants.DEPTH_KEY, alignmentContext.size());

                // write to VCF
                VariantContext vc = vcb.make();
                vcfWriter.add(vc);
            }
        }
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null )
           vcfWriter.close();
    }
}
