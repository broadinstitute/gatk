package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class ConsensusReadsByRealignment extends ReadWalker {
    final List<GATKRead> reads = new ArrayList<>(30);

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        reads.add(read);
    }

    @Override
    public Object onTraversalSuccess() {
        // call region in m2engine
        final SimpleInterval interval = new SimpleInterval(reads.get(0).getContig(), reads.get(0).getStart(), reads.get(0).getEnd());
        final AssemblyRegion assemblyRegion = new AssemblyRegion(interval,0, getHeaderForReads());
        assemblyRegion.addAll(reads);
        final SampleList sample = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(getHeaderForReads())));
        final ReferenceSequenceFile referenceReader = AssemblyBasedCallerUtils.createReferenceReader(Utils.nonNull(referenceArguments.getReferenceFileName()));
        final M2ArgumentCollection mtac = new M2ArgumentCollection(); // AssemblyBasedCallerArgumentCollect suffices, too.
        final ReadThreadingAssembler assemblyEngine = mtac.createReadThreadingAssembler();
        final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(mtac.smithWatermanImplementation);

        final AssemblyResultSet assemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyRegion, Collections.emptyList(),
                mtac, getHeaderForReads(), sample, logger, referenceReader, assemblyEngine, aligner, false);
        int d = 3;

        return "SUCCESS";
    }
}
