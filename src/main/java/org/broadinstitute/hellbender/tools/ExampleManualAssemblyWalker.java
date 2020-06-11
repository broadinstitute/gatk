package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Assembly example tool",
        oneLineSummary = "Assembly example tool",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExampleManualAssemblyWalker extends ReadWalker {

    private static final Logger logger = LogManager.getLogger(ExampleManualAssemblyWalker.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public GATKPath output;

    private SAMFileGATKReadWriter outputWriter;

    private SampleList sampleList;
    private ReferenceSequenceFile referenceReader;
    private ReadThreadingAssembler assembler;
    private SmithWatermanAligner aligner;
    private ReadLikelihoodCalculationEngine likelihoodEngine;

    // These arguments are not exposed on the CLI by the tool, since they are not annotated as an ArgumentCollection
    private final HaplotypeCallerArgumentCollection assemblySettings = new HaplotypeCallerArgumentCollection();

    private final List<GATKRead> readsForAssembly = new ArrayList<>();
    private GATKRead previousRead;

    @Override
    public boolean requiresReference() {
        return true;
    }

    // You can add filters and/or read transformers via the methods below:
    
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return super.getDefaultReadFilters();
    }

    @Override
    public ReadTransformer makePreReadFilterTransformer() {
        return super.makePreReadFilterTransformer();
    }

    @Override
    public ReadTransformer makePostReadFilterTransformer() {
        return super.makePostReadFilterTransformer();
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);

        sampleList = new IndexedSampleList(ReadUtils.getSamplesFromHeader(getHeaderForReads()));
        referenceReader = new CachingIndexedFastaSequenceFile(referenceArguments.getReferenceSpecifier());
        assembler = assemblySettings.assemblerArgs.makeReadThreadingAssembler();
        aligner = SmithWatermanAligner.getAligner(assemblySettings.smithWatermanImplementation);
        likelihoodEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(assemblySettings.likelihoodArgs);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        if ( startNewRegion(previousRead, read) ) {
            if ( ! readsForAssembly.isEmpty() ) {
                final SimpleInterval regionSpan = IntervalUtils.getSpanningInterval(readsForAssembly);
                final AssemblyRegion assemblyRegion = new AssemblyRegion(regionSpan, regionSpan, true, getHeaderForReads());
                assemblyRegion.addAll(readsForAssembly);

                processRegion(assemblyRegion);
            }

            readsForAssembly.clear();
        }

        readsForAssembly.add(read);
        previousRead = read;
    }

    @Override
    public Object onTraversalSuccess() {
        // Any final actions go here...
        return null;
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }

    private boolean startNewRegion(final GATKRead previousRead, final GATKRead currentRead) {
        return false; // Placeholder
    }

    private void processRegion(final AssemblyRegion region) {
        final AssemblyResultSet assemblyResult = AssemblyBasedCallerUtils.assembleReads(region, Collections.emptyList(), assemblySettings, getHeaderForReads(), sampleList, logger, referenceReader, assembler, aligner, ! assemblySettings.doNotCorrectOverlappingBaseQualities);

        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();

        System.out.println("\nThere were " + haplotypes.size() + " haplotypes found in region " + region + ". Here they are:");
        for (String haplotype : haplotypes.stream().map(haplotype -> haplotype.toString()).sorted().collect(Collectors.toList())) {
            System.out.println(haplotype);
        }

        final Map<String,List<GATKRead>> readsBySample = AssemblyBasedCallerUtils.splitReadsBySample(sampleList, getHeaderForReads(), region.getReads());
        final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods =
                likelihoodEngine.computeReadLikelihoods(assemblyResult, sampleList, readsBySample);

        // Get the best haplotype for each read:
        final Collection<AlleleLikelihoods<GATKRead, Haplotype>.BestAllele> bestAlleles = readLikelihoods.bestAllelesBreakingTies(AssemblyBasedCallerUtils.HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY);

        // The actual HaplotypeCaller/Mutect2 then go on to call into their respective genotyping engines.
        // See, eg., HaplotypeCallerEngine.callRegion()

        // For an alternate assembly option not used by the HC/Mutect2 but used by the SV tools,
        // see the FermiLiteAssembler
    }

}
