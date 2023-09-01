package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.AssemblerOffRamp;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.PostAssemblerOnRamp;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.PostFilterOnRamp;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.PreFilterOffRamp;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;

/**
 * This is a specialized HaplotypeCaller tool, designed to allow for breaking the monolithic haplotype
 * caller process into smaller discrete steps. This version allows for stopping the process at after (or before)
 * a specific step and saving a state file (an off-ramp) that can later be used to restart the process from the
 * same step (an on-ramp).
 *
 * The major steps of the haplotype caller are: assembly, filtering and genotyping
 *
 * At this off ramps are implemented before and after the assembler and before the filtering
 * On ramps are implemented after the assembler and after the filtering.
 *
 * For a description of the HaplotypeCaller tool, see {@link HaplotypeCaller}
 *
 * <h3>Usage examples</h3>
 *
 * <p>This tool implements the same arguments as the the base {@link HaplotypeCaller} while adding a small
 * number of argument for ramp control (see below). Usage examples are based on HaplotypeCaller usage examples</p>
 *
 * <br />
 * <h4>Off-ramp example</h4>
 * Note that pre_assembler.zip is an output file
 * <pre>
 * gatk --java-options "-Xmx4g" HaplotypeCaller  \
 *   // standard HaplotypeCaller arguments and their values
 *   --off-ramp-type PRE_ASSEMBLER_OFF
 *   --off-ramp-file pre_assembler.zip
 * </pre>
 *
 *
 * <br />
 * <h4>On-ramp example</h4>
 * Note that post_assembler.zip is an input file
 * <pre>
 * gatk --java-options "-Xmx4g" HaplotypeCaller  \
 *   // standard HaplotypeCaller arguments and their values
 *   --on-ramp-type POST_ASSEMBLER_ON
 *   --on-ramp-file post_assembler.zip
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "specialized HaplotypeCaller tool, designed to allow for breaking the monolithic haplotype caller " +
                "process into smaller discrete steps. NOTE: this is a debugging tool!",
        oneLineSummary = "Call germline SNPs and indels via local re-assembly of haplotypes (ramped version)",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class RampedHaplotypeCaller extends HaplotypeCaller {

    @ArgumentCollection
    private RampedHaplotypeCallerArgumentCollection rpArgs = new RampedHaplotypeCallerArgumentCollection();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
    }

    @Override
    protected HaplotypeCallerEngine buildHaplotypeCallerEngine(final HaplotypeCallerArgumentCollection hcArgs, final AssemblyRegionArgumentCollection assemblyRegionArgs, final boolean createOutputBamIndex, final boolean createOutputBamMD5, final SAMFileHeader headerForReads, final CachingIndexedFastaSequenceFile referenceReader, final VariantAnnotatorEngine variantAnnotatorEngine) {
        return new RampedHaplotypeCallerEngine(hcArgs, assemblyRegionArgs, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), getReferenceReader(referenceArguments), variantAnnotatorEngine, rpArgs);
    }

    @Override
    public boolean nonRandomDownsamplingMode() {
        return true;
    }


}
