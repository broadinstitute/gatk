package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class PairHMMErrorCorrectionEngine extends ErrorCorrectionEngineBase {
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private M2ArgumentCollection mtac = new M2ArgumentCollection();
    private HaplotypeCallerArgumentCollection hcac = new HaplotypeCallerArgumentCollection();

    PairHMMErrorCorrectionEngine(SAMFileHeader header, ReferenceInputArgumentCollection referenceArguments, SAMFileGATKReadWriter outputWriter) {
        super(header, referenceArguments, outputWriter);
        hcac.standardArgs.genotypeArgs.samplePloidy = 1;
        mtac.likelihoodArgs.pairHMM = PairHMM.Implementation.LOGLESS_CACHING;
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(mtac.likelihoodArgs);
    }

    /**
     * Reassemble reads and choose the best haplotype as the consensus read, the way haplotypcaller
     * finds variants with ploidy=1
     ***/
    @Override
    protected GATKRead getIndelConsensusRead(List<GATKRead> reads, ReferenceContext referenceContext) {
        final SampleList indexedSampleList = new IndexedSampleList(new ArrayList<>(samplesList));
        /** {@link ReadThreadingAssembler.findBestPaths() } requires a graph. **/
        // PairHMM is different.
        // Learning of the PCR error likelihood should be done by deep learning. This should be the self-motivated project.
        // Find the consensus read
        final boolean assembleReads = false;
        final AssemblyResultSet assemblyResult = new AssemblyResultSet();

        // Can I just go from reads to haplotypes to events? Yes.
        final int regionStart = reads.stream().mapToInt(GATKRead::getStart).min().getAsInt();
        final int regionEnd = reads.stream().mapToInt(GATKRead::getEnd).max().getAsInt();
        final int REFERENCE_PADDING = 50; // We need padding because in the event of insertion we might need bases in the reference that falls outside of the tight window
        // must take the min of end of contig and regionEnd + padding
        final SimpleInterval refInterval = new SimpleInterval(referenceContext.getContig(), Math.max(regionStart - REFERENCE_PADDING, 0), regionEnd + REFERENCE_PADDING);
        final SimpleInterval justCurious = referenceContext.getInterval();
        // final Set<Haplotype> haplotypes = readsToHaplotypeSet(reads, refInterval);
        final List<ReadAndHaplotype> readAndHaplotypes = getReadAndHaplotypes(reads, refInterval);
        final List<Haplotype> haplotypes = ReadAndHaplotype.getHaplotypeList(readAndHaplotypes);
        // Must create an empty assembly result set with haplotypes in order to call computeReadLikelihoods(), which
        // does not use any other information in the assembly result set.
        readAndHaplotypes.forEach(rah -> assemblyResult.add(rah.getHaplotype()));
        assemblyResult.setPaddedReferenceLoc(refInterval);
        final ReferenceSequenceFile referenceReader = AssemblyBasedCallerUtils.createReferenceReader(Utils.nonNull(referenceArguments.getReferenceFileName()));
        // final byte[] referenceBases = AssemblyRegion.getReference(referenceReader, 0, refInterval);
        // assemblyResult.setFullReferenceWithPadding(referenceBases);

        // assemblyResult = doAssembly(reads, indexedSampleList);
        // final SortedSet<VariantContext> set = assemblyResult.getVariationEvents(1);

        return null;
    }

    private AssemblyResultSet doAssembly(final List<GATKRead> reads, final SampleList indexedSampleList){
        final int maxDeletionLength = 100;
        final int start = reads.stream().mapToInt(GATKRead::getStart).min().getAsInt();
        final int end = reads.stream().mapToInt(GATKRead::getEnd).max().getAsInt();
        final SimpleInterval interval = new SimpleInterval(reads.get(0).getContig(), start, end); // is this wrong? should I get the min of the start positions?
        final AssemblyRegion assemblyRegion = new AssemblyRegion(interval, 0, header);
        assemblyRegion.addAll(reads);
        final ReferenceSequenceFile referenceReader = AssemblyBasedCallerUtils.createReferenceReader(Utils.nonNull(referenceArguments.getReferenceFileName()));
        final M2ArgumentCollection mtac = new M2ArgumentCollection(); // AssemblyBasedCallerArgumentCollect suffices, too.
        final ReadThreadingAssembler assemblyEngine = mtac.createReadThreadingAssembler();
        final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(mtac.smithWatermanImplementation);

        // ts: logger being null is not good.
        final AssemblyResultSet assemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyRegion, Collections.emptyList(),
                mtac, header, indexedSampleList, null, referenceReader, assemblyEngine, aligner, false);
        return assemblyResult;
    }

    /** 12/4/19 AF-based filter, and unnecessarily complicated code that handles multi-ploidy case, etc,
     * opting for a simple implementation using the whole haplotype
     * **/
    public void likelihoodToGenotypesHC(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                        final AssemblyResultSet assemblyResult){
        final SampleList indexedSampleList = new IndexedSampleList(new ArrayList<>(samplesList));
        final HaplotypeCallerGenotypingEngine hcGenotypingEngine = new HaplotypeCallerGenotypingEngine(hcac, indexedSampleList, ! hcac.doNotRunPhysicalPhasing);
        final CalledHaplotypes calledHaplotypes = hcGenotypingEngine.assignGenotypeLikelihoods(assemblyResult.getHaplotypeList(), readLikelihoods, new HashMap<>(),
                assemblyResult.getFullReferenceWithPadding(), // ts: ref bases. Where should it start/end?
                assemblyResult.getPaddedReferenceLoc(), // ts: what about reference loc?
                assemblyResult.getRegionForGenotyping().getSpan(),
                null, new ArrayList<>(), false, 1, header, false);
        // Genotypes are empty. How can I go from calls to genotype likelihoods?
    }

    private void getEventsWithinDuplciateSetSketch(){
        // EventMap.buildEventMapsForHaplotypes(new ArrayList<>(haplotypes), referenceContext.getBases(), referenceContext.getInterval(), false, 1);
        // final SortedSet<VariantContext> events = AssemblyResultSet.getAllVariantContexts(new ArrayList<>(haplotypes));
        // final List<VariantContext> vc = AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(); // ts: this is events at a particular locus
    }


}
