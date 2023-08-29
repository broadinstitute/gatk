package org.broadinstitute.hellbender.tools.walkers.featuremapping;


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
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;

import java.util.*;
/**
 * A tool to align list of reads to a list of haplotypes. The alignment score is calculated based on assumption
 * that the reads were generated from one of the haplotypes and only sequencing errors. Thus, the alignment score
 * is exactly the likelihood of the read given haplotype that is calculated in HaplotypeCaller.
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> BAM/CRAM file</li>
 *     <li> Fasta of haplotypes (fa) </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> TSV file that is either in the format of read x haplotype matrix or in
 *          a "concise" format.
 *     </li>
 * </ul>
 *
 *    Since the tool was designed for alignment of the flow-based reads, it currently supports two alignment engines:
 *    FlowPairHMM and FlowBasedAlignment (FBA), but can be easily extended.
 *    At present, there are two output formats that can be specified using parameter --output-format: extended and concise.
 *    The extended format contains a readxhaplotype matrix that shows alignment score of each read versus each haplotype.
 *    Condensed format will contain the following columns for each processed read:
 *    likelihood score, the best haplotype, the second best haplotype and the difference of alignment scores
 *    between the best and the second best haplotype. In addition, as in many cases most of the reads are
 *    coming from the "reference" haplotype we can also output the distance from the (marked) reference haplotype
 *
 * <h3>Usage examples</h3>
 * <pre>
*             gatk FlowPairHMMAlignReadsToHaplotypes \
 *            -H ~{haplotype_list} -O ~{base_file_name}.matches.tsv \
 *            -I ~{input_bam} --flow-use-t0-tag -E FBA \
 *            --flow-fill-empty-bins-value 0.00001 --flow-probability-threshold 0.00001 \
 *            --flow-likelihood-optimized-comp
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */

@CommandLineProgramProperties(
        summary = "Align Reads to Haplotypes using FlowBasedPairHMM",
        oneLineSummary = "Produces readxhaplotype matrix with likelihoods of read / haplotype",
        programGroup = FlowBasedProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public final class FlowPairHMMAlignReadsToHaplotypes extends ReadWalker {

    private static final Logger logger = LogManager.getLogger(FlowPairHMMAlignReadsToHaplotypes.class);

    @Argument(fullName = "haplotypes", shortName = "H", doc="Fasta file with haplotypes")
    public GATKPath haplotypesFa;

    @Argument(fullName = "ref-haplotype", doc="Fasta file with haplotypes", optional = true)
    public String refHaplotypeName;

    @Argument(fullName = "output", shortName = "O", doc="Read x haplotype log-likelihood matrix")
    public String output;

    @Argument(fullName = "concise-output-format", doc="concise or expanded output format: " +
            "expanded - output full read x haplotype, concise - output for each read best haplotype and " +
            "score differences from the next best and the reference haplotype, default: false (expanded format)", optional=true)
    public boolean conciseOutputFormat = false;

    @Argument(fullName = "aligner", shortName = "E", doc="Aligner: FlowBasedHMM or FlowBasedAligner (FlowBased)")
    public ReadLikelihoodCalculationEngine.Implementation alignerName = ReadLikelihoodCalculationEngine.Implementation.FlowBased;

    @ArgumentCollection
    public FlowBasedAlignmentArgumentCollection fbargs = new FlowBasedAlignmentArgumentCollection();

    @ArgumentCollection
    private static final LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();
    SAMFileHeader sequenceHeader;

    ReadLikelihoodCalculationEngine aligner;

    AlleleLikelihoodWriter outputWriter;

    List<Haplotype> haplotypeList;
    Map<String, String> haplotypeToName;
    List<GATKRead> readBuffer;
    boolean writeHeader=true;
    int readCount = 0;
    final int BUFFER_SIZE_LIMIT = 50;
    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        sequenceHeader = getHeaderForReads();
        if (!conciseOutputFormat) {
            outputWriter = new AlleleLikelihoodWriter(IOUtils.getPath(output), null);
        } else {
            outputWriter = new ConciseAlleleLikelihoodWriter(IOUtils.getPath(output));
        }
        haplotypeList = new ArrayList<>();
        fbargs.keepBoundaryFlows = true;
        fbargs.trimToHaplotype = false;
        fbargs.exactMatching = true;
        FastaSequenceFile haplotypeFile = new FastaSequenceFile(haplotypesFa.toPath(), false);
        ReferenceSequence seq = haplotypeFile.nextSequence();
        readBuffer = new ArrayList<>();
        haplotypeToName = new HashMap<>();
        while (seq != null ){
            Haplotype hap = new Haplotype(seq.getBases(), seq.getName().equals(refHaplotypeName));
            haplotypeList.add(hap);
            haplotypeToName.put(new String(seq.getBases()), seq.getName());
            seq = haplotypeFile.nextSequence();
        }

        if (alignerName == ReadLikelihoodCalculationEngine.Implementation.FlowBasedHMM) {
            aligner = new FlowBasedHMMEngine(fbargs, (byte) likelihoodArgs.gcpHMM,
                    10000, likelihoodArgs.expectedErrorRatePerBase, likelihoodArgs.pcrErrorModel,
                    null,
                    likelihoodArgs.enableDynamicReadDisqualification, likelihoodArgs.readDisqualificationThresholdConstant,
                    likelihoodArgs.minUsableIndelScoreToUse, (byte) likelihoodArgs.flatDeletionPenalty,
                    (byte) likelihoodArgs.flatInsertionPenatly);
        } else if (alignerName == ReadLikelihoodCalculationEngine.Implementation.FlowBased) {
            aligner = new FlowBasedAlignmentLikelihoodEngine(fbargs,
                    10000, 1,
                    false, 0);
        } else {
            throw new RuntimeException("Accepted engines are FlowBasedHMM or FlowBased");
        }
    }
    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext){

        readBuffer.add(read);
        readCount++;

        if (readBuffer.size() == BUFFER_SIZE_LIMIT){
            AlleleLikelihoods<GATKRead, Haplotype> readByHaplotypeMatrix = calculateLikelihoods(readBuffer);
            outputWriter.writeAlleleLikelihoodsAsMatrix(readByHaplotypeMatrix, haplotypeToName, writeHeader, readCount - BUFFER_SIZE_LIMIT);
            writeHeader=false;
            readBuffer.clear();
        }
    }

    @Override
    public Object onTraversalSuccess(){
        AlleleLikelihoods<GATKRead, Haplotype> readByHaplotypeMatrix = calculateLikelihoods(readBuffer);
        outputWriter.writeAlleleLikelihoodsAsMatrix(readByHaplotypeMatrix, haplotypeToName, writeHeader, readCount - BUFFER_SIZE_LIMIT);
        writeHeader=false;
        readBuffer.clear();
        return null;
    }
    private AlleleLikelihoods<GATKRead,Haplotype> calculateLikelihoods(final List<GATKRead> reads){
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));
        Map<String, List<GATKRead>> tmpMap = new HashMap<>();
        tmpMap.put("sm1", reads);
        if (alignerName == ReadLikelihoodCalculationEngine.Implementation.FlowBasedHMM) {
            return ((FlowBasedHMMEngine)aligner).computeReadLikelihoods(haplotypeList, sequenceHeader, samples, tmpMap, false, false);
        } else if (alignerName == ReadLikelihoodCalculationEngine.Implementation.FlowBased){
            return ((FlowBasedAlignmentLikelihoodEngine)aligner).computeReadLikelihoods(haplotypeList, sequenceHeader, samples, tmpMap, false, false);
        } else {
            return null;
        }
    }
}
