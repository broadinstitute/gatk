package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ContigChimericAlignmentIterativeInterpreter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection;

/**
 * (Internal) Examines aligned contigs from local assemblies and calls structural variants
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It packages structural
 * variant calling as a separate tool, independent of the generation of local assemblies.
 * Most researchers will run StructuralVariationDiscoveryPipelineSpark, which both generates local assemblies
 * of interesting genomic regions, and then calls structural variants from these assemblies.</p>
 * <p>This tool takes a SAM/BAM/CRAM containing the alignments of assembled contigs from local assemblies
 * and searches it for split alignments indicating the presence of structural variations. To do so the tool parses
 * primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV,
 * two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than
 * or equal to min-alignment-length. Imprecise variants with approximate locations are also called.</p>
 * <p>The input file is typically the output file produced by FindBreakpointEvidenceSpark.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of assembled contigs or long reads aligned to reference.</li>
 *     <li>The reference to which the contigs have been aligned.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A vcf file describing the discovered structural variants.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk DiscoverVariantsFromContigAlignmentsSAMSpark \
 *     -I assemblies.sam \
 *     -R reference.2bit \
 *     -O structural_variants.vcf
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 *
 * <h3>Notes</h3>
 * <p>The reference is broadcast by Spark, and must therefore be a .2bit file due to current restrictions.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Examines aligned contigs from local assemblies and calls structural variants",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It packages structural" +
        " variant calling as a separate tool, independent of the generation of local assemblies." +
        " Most researchers will run StructuralVariationDiscoveryPipelineSpark, which both generates local assemblies" +
        " of interesting genomic regions, and then calls structural variants from these assemblies." +
        " This tool takes a SAM/BAM/CRAM containing the alignments of assembled contigs from local assemblies" +
        " and searches it for split alignments indicating the presence of structural variations. To do so the tool parses" +
        " primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV," +
        " two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than" +
        " or equal to min-alignment-length. Imprecise variants with approximate locations are also called.\n" +
        " The input file is typically the output file produced by FindBreakpointEvidenceSpark.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class DiscoverVariantsFromContigAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(DiscoverVariantsFromContigAlignmentsSAMSpark.class);

    @ArgumentCollection
    private final DiscoverVariantsFromContigAlignmentsSparkArgumentCollection discoverStageArgs =
            new DiscoverVariantsFromContigAlignmentsSparkArgumentCollection();

    @Argument(doc = "prefix for discovery (non-genotyped) VCF; sample name will be appended after the provided argument, followed by \"_inv_del_ins.vcf\"",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String prefixForOutput;


    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.MAPPED);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        validateParams();

        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast =
                StructuralVariationDiscoveryPipelineSpark.broadcastCNVCalls(ctx, getHeaderForReads(),
                        discoverStageArgs.cnvCallsFile);

        final String vcfOutputPath = getVcfOutputPath();

        final SvDiscoveryInputMetaData svDiscoveryInputMetaData =
                new SvDiscoveryInputMetaData(ctx, discoverStageArgs, null, vcfOutputPath,
                        null, null, null,
                        cnvCallsBroadcast,
                        getHeaderForReads(), getReference(), getDefaultToolVCFHeaderLines(), localLogger);

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new SvDiscoverFromLocalAssemblyContigAlignmentsSpark
                        .SAMFormattedContigAlignmentParser(getReads(),
                                                           svDiscoveryInputMetaData.getSampleSpecificData().getHeaderBroadcast().getValue(), true)
                        .getAlignedContigs();

        // assembly-based breakpoints
        List<VariantContext> annotatedVariants =
                ContigChimericAlignmentIterativeInterpreter
                        .discoverVariantsFromChimeras(svDiscoveryInputMetaData, parsedContigAlignments);

        final SAMSequenceDictionary refSeqDictionary = svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast().getValue();
        SVVCFWriter.writeVCF(annotatedVariants, vcfOutputPath, refSeqDictionary, getDefaultToolVCFHeaderLines(), localLogger);
    }


    private void validateParams() {
        discoverStageArgs.validate();
    }

    private String getVcfOutputPath() {
        if ( Files.exists(Paths.get(prefixForOutput)) ) {
            if ( Files.isDirectory(Paths.get(prefixForOutput)) ) // existing directory
                return prefixForOutput + (prefixForOutput.endsWith("/") ? "" : "/") + SVUtils.getSampleId(getHeaderForReads()) + "_inv_del_ins.vcf";
            else
                throw new UserException("Provided prefix for output is pointing to an existing file: " + prefixForOutput); // to avoid accidental override of a file
        } else { // prefixForOutput doesn't point to an existing file or directory
            return prefixForOutput + (prefixForOutput.endsWith("/") ? "" : "_") + SVUtils.getSampleId(getHeaderForReads()) + "_inv_del_ins.vcf";
        }
    }
}
