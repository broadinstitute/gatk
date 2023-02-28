package org.broadinstitute.hellbender.tools.walkers.groundtruth;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AlleleLikelihoodWriter;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedHMMEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LikelihoodEngineArgumentCollection;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;

import java.util.*;

@CommandLineProgramProperties(
        summary = "Align Reads to Haplotypes using FlowBasedPairHMM",
        oneLineSummary = "Produces readxhaplotype matrix with likelihoods of read / haplotype",
        programGroup = FlowBasedProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public class FlowPairHMMAlignReadsToHaplotypes extends ReadWalker {

    private static final Logger logger = LogManager.getLogger(FlowPairHMMAlignReadsToHaplotypes.class);

    @Argument(fullName = "haplotypes", shortName = "R", doc="Fasta file with haplotypes")
    public GATKPath haplotypes_fa;

    @Argument(fullName = "input", shortName = "I", doc="Input BAM/CRAM file with reads")
    public GATKPath reads;

    @Argument(fullName = "output", shortName = "O", doc="Readxhaplotype log-likelihood matrix")
    public String output;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
    LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();
    SAMFileHeader sequenceHeader;
    final FlowBasedHMMEngine aligner = new FlowBasedHMMEngine(fbargs, (byte) likelihoodArgs.gcpHMM,
            10000, likelihoodArgs.expectedErrorRatePerBase, likelihoodArgs.pcrErrorModel,
            likelihoodArgs.dontUseDragstrPairHMMScores ? null : DragstrParamUtils.parse(likelihoodArgs.dragstrParams),
            likelihoodArgs.enableDynamicReadDisqualification, likelihoodArgs.readDisqualificationThresholdConstant,
            likelihoodArgs.minUsableIndelScoreToUse, (byte) likelihoodArgs.flatDeletionPenalty,
            (byte) likelihoodArgs.flatInsertionPenatly);

    AlleleLikelihoodWriter outputWriter;

    List<Haplotype> haplotypeList;

    List<GATKRead> readBuffer;
    final int BUFFER_SIZE_LIMIT = 1000;
    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        sequenceHeader = getHeaderForReads();
        outputWriter = new AlleleLikelihoodWriter(IOUtils.getPath(output), null);
        haplotypeList = new ArrayList<>();
        FastaSequenceFile haplotypeFile = new FastaSequenceFile(haplotypes_fa.toPath(), false);
        ReferenceSequence seq = haplotypeFile.nextSequence();
        byte[] flowOrder = FlowBasedReadUtils.getReadFlowOrder(sequenceHeader, null);
        readBuffer = new ArrayList<>();
        while (seq != null ){
            Haplotype hap = new Haplotype(seq.getBases());
            haplotypeList.add(hap);
            seq = haplotypeFile.nextSequence();
        }
    }
    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext){
        readBuffer.add(read);
        if (readBuffer.size() == BUFFER_SIZE_LIMIT){
            AlleleLikelihoods<GATKRead, Haplotype> readByHaplotypeMatrix = calculateLikelihoods(readBuffer);
            outputWriter.writeAlleleLikelihoods(readByHaplotypeMatrix);
            readBuffer.clear();
        }
    }

    private AlleleLikelihoods<GATKRead,Haplotype> calculateLikelihoods(final List<GATKRead> reads){
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));
        Map<String, List<GATKRead>> tmpMap = new HashMap<>();
        tmpMap.put("sm1", reads);
        return aligner.computeReadLikelihoods(haplotypeList, sequenceHeader, samples, tmpMap, false);
    }
}
