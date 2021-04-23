package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.stream.Collectors;

public class SVDepthOnlyCallDefragmenter extends LocatableClusterEngine<SVCallRecordWithEvidence> {

    private final double minSampleOverlap;
    private final double paddingFraction;
    private static final double DEFAULT_MIN_SAMPLE_OVERLAP = 0.9;

    @VisibleForTesting
    protected static final double DEFAULT_PADDING_FRACTION = 0.5;

    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary) {
        this(dictionary, DEFAULT_PADDING_FRACTION, DEFAULT_MIN_SAMPLE_OVERLAP, null);
    }

    //for single-sample clustering case
    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary, double paddingFraction,
                                       double minSampleOverlap, List<GenomeLoc> coverageIntervals) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE, coverageIntervals);
        this.minSampleOverlap = minSampleOverlap;
        this.paddingFraction = paddingFraction;
    }

    /**
     * Find a single call representative of all the calls in the cluster
     * @param cluster   the events that are clustered together
     * @return  a call encompassing all the cluster's events and containing all the algorithms and genotypes
     */
    @Override
    protected SVCallRecordWithEvidence flattenCluster(final Collection<SVCallRecordWithEvidence> cluster) {
        final int newStart =  cluster.stream().mapToInt(SVCallRecordWithEvidence::getStart).min().getAsInt();
        final int newEnd = cluster.stream().mapToInt(SVCallRecordWithEvidence::getEnd).max().getAsInt();
        final SVCallRecordWithEvidence exampleCall = cluster.iterator().next();
        final int length = newEnd - newStart + 1;  //+1 because GATK intervals are inclusive
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList()); //should be depth only
        final List<Genotype> clusterGenotypes = deduplicateGenotypes(cluster.stream().flatMap(v -> v.getGenotypes().stream()).collect(Collectors.toList()));
        return new SVCallRecordWithEvidence(exampleCall.getContig(), newStart, exampleCall.getStartStrand(),
                exampleCall.getEndContig(), newEnd, exampleCall.getEndStrand(), exampleCall.getType(), length, algorithms, clusterGenotypes,
                exampleCall.getStartSplitReadSites(), exampleCall.getEndSplitReadSites(), exampleCall.getDiscordantPairs());
    }

    protected List<Genotype> deduplicateGenotypes(final List<Genotype> clusterGenotypes) {
        final Set<String> samples = clusterGenotypes.stream()
                .filter(Genotype::isCalled)
                .map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
        if (samples.size() == clusterGenotypes.size()) {
            return clusterGenotypes;
        }
        final List<Genotype> mergedGenotypes = new ArrayList<>();
        final Map<String, List<Genotype>> moreGenotypes = new LinkedHashMap<>();
        for (final Genotype g : clusterGenotypes) {
            if (!moreGenotypes.containsKey(g.getSampleName())) {
                final List<Genotype> newList = new ArrayList<>();
                newList.add(g);
                moreGenotypes.put(g.getSampleName(), newList);
            } else {
                (moreGenotypes.get(g.getSampleName())).add(g);
            }
        }
        for (final List<Genotype> gList : moreGenotypes.values()) {
            if (gList.size() > 1) {
                mergedGenotypes.add(defragmentGenotypes(gList));
            } else if (gList.size() == 1) {
                mergedGenotypes.add(gList.get(0));
            }
        }
        return mergedGenotypes;
    }

    /**
     *
     * @param genotypesForSameSample
     * @return
     */
    private Genotype defragmentGenotypes(final List<Genotype> genotypesForSameSample) {
        final String sampleName = genotypesForSameSample.get(0).getSampleName();
        if (!genotypesForSameSample.stream().allMatch(g -> g.getSampleName().equals(sampleName))) {
            throw new IllegalArgumentException("This method expects a list of genotypes from the same sample, " +
                    "but not all input genotypes represent sample " + sampleName + ".");
        }
        final GenotypeBuilder gb = new GenotypeBuilder(genotypesForSameSample.get(0));
        //For now just make sure genotypes have the same copy number -- qualities will be recalculated elsewhere
        gb.noAttributes();
        final int copyNumber = Integer.parseInt(genotypesForSameSample.get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString());
        if (genotypesForSameSample.stream().allMatch(g -> Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) == copyNumber)) {
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
            return gb.make();
        } else {
            throw new IllegalArgumentException("This method will only merge genotypes with the same copy number. Expected all genotypes to be copy number " + copyNumber + ".");
        }
    }

    /**
     * Determine if two calls should cluster based on their padded intervals and genotyped samples
     * @param a
     * @param b
     * @return true if the two calls should be in the same cluster
     */
    @Override
    protected boolean clusterTogether(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!isDepthOnlyCall(a) || !isDepthOnlyCall(b)) return false;
        Utils.validate(a.getContig().equals(a.getEndContig()), "Call A is depth-only but interchromosomal");
        Utils.validate(b.getContig().equals(b.getEndContig()), "Call B is depth-only but interchromosomal");
        if (!a.getType().equals(b.getType())) return false;
        final Set<String> sharedSamples = new LinkedHashSet<>(a.getSamples());
        sharedSamples.retainAll(b.getSamples());
        final double sampleOverlap = Math.min(sharedSamples.size() / (double) a.getSamples().size(), sharedSamples.size() / (double) b.getSamples().size());
        if (sampleOverlap < minSampleOverlap) return false;
        //in the single-sample case, match copy number strictly if we're looking at the same sample
        boolean copyNumbersAgree = true;
        if (a.getGenotypes().size() == 1 && b.getGenotypes().size() == 1 && a.getGenotypes().get(0).getSampleName().equals(b.getGenotypes().get(0).getSampleName())) {
            if (a.getGenotypes().get(0).hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) && b.getGenotypes().get(0).hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) &&
                !(a.getGenotypes().get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).equals(b.getGenotypes().get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)))) {
                copyNumbersAgree = false;
            }
        }
        return getClusteringInterval(a, null)
                .overlaps(getClusteringInterval(b, null)) && copyNumbersAgree;
    }


    /**
     * Determine an overlap interval for clustering using padding specified at object construction
     * Returned interval represents the interval in which the start position of a new event must fall in order to be
     * added to the cluster (including the new event)
     * @param call  new event to be clustered
     * @param currentClusterInterval    the cluster of interest, may be null
     * @return  an interval describing the cluster after the new event is added
     */
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecordWithEvidence call, final SimpleInterval currentClusterInterval) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getCallInterval(call);
        final int paddedCallStart, paddedCallEnd;
        if (genomicToBinMap != null) {
            final GenomeLoc callStart = parser.createGenomeLoc(call.getContig(), call.getStart(), call.getStart());
            final GenomeLoc callEnd = parser.createGenomeLoc(call.getContig(), call.getEnd(), call.getEnd());
            //first interval that is equal to or "greater than" the call start, such that the start of the bin should match the call start, with a little wiggle room
            final Map.Entry<GenomeLoc, Integer> startBin = genomicToBinMap.ceilingEntry(callStart);
            if (startBin == null) {
                throw new UserException.BadInput("Call start " + callStart + " for  call " + call.prettyPrint() + " not found in model call intervals.");
            }
            final int callStartIndex = startBin.getValue();
            //last interval that is equal to or "less than" the call start, such that the end of the bin should match the call end
            final Map.Entry<GenomeLoc, Integer> endBin = genomicToBinMap.floorEntry(callEnd);
            if (endBin == null) {
                throw new UserException.BadInput("Call end " + callEnd + " for call " + call.prettyPrint() + " not found in model call intervals.");
            }
            final int callEndIndex = endBin.getValue();
            final int callBinLength = callEndIndex - callStartIndex + 1;
            if (callBinLength <= 0) {
                throw new UserException.BadInput("Copy number call at " + call.getContig() + ":" + call.getStart() + "-"
                        + call.getEnd() + " does not align with supplied model calling intervals. Use the filtered intervals input from GermlineCNVCaller for this cohort/model.");
            }
            final int paddedStartIndex = Math.max(callStartIndex - (int)Math.round(callBinLength * paddingFraction), 0);
            if (coverageIntervals.get(paddedStartIndex).getContig().equals(callStart.getContig())) {
                paddedCallStart = coverageIntervals.get(paddedStartIndex).getStart();
            } else {
                paddedCallStart = callStart.getStart();
            }
            final int paddedEndIndex = Math.min(callEndIndex + (int)Math.round(callBinLength * paddingFraction), genomicToBinMap.size() - 1);
            if (coverageIntervals.get(paddedEndIndex).getContig().equals(callEnd.getContig())) {
                paddedCallEnd = coverageIntervals.get(paddedEndIndex).getEnd();
            } else {
                paddedCallEnd = callEnd.getEnd();
            }
        } else {
            paddedCallStart = (int) (callInterval.getStart() - paddingFraction * callInterval.getLengthOnReference());
            paddedCallEnd = (int) (callInterval.getEnd() + paddingFraction * callInterval.getLengthOnReference());
        }
        final int contigLength = dictionary.getSequence(call.getContig()).getSequenceLength();
        if (currentClusterInterval == null) {
            return IntervalUtils.trimIntervalToContig(call.getContig(), paddedCallStart, paddedCallEnd, contigLength);
        }
        //NOTE: this is an approximation -- padding should be based on the length of the call plus currentClusterIntervals
        final int newMinStart = Math.min(paddedCallStart, currentClusterInterval.getStart());
        final int newMaxEnd = Math.max(paddedCallEnd, currentClusterInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(call.getContig(), newMinStart, newMaxEnd, contigLength);
    }

    // Not used for single-linkage clustering
    @Override
    protected boolean itemsAreIdentical(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        throw new GATKException.ShouldNeverReachHereException("Deduplication should not be called for single-linkage clustering in depth-only defragmentation.");
    }

    // Not used for single-linkage clustering
    @Override
    protected SVCallRecordWithEvidence deduplicateIdenticalItems(final Collection<SVCallRecordWithEvidence> items) {
        return null;
    }

    private SimpleInterval getCallInterval(final SVCallRecordWithEvidence call) {
        return new SimpleInterval(call.getContig(), call.getStart(), call.getEnd());
    }

    public static boolean isDepthOnlyCall(final SVCallRecord call) {
        if (call.getAlgorithms().isEmpty()) return false;
        for (final String alg : call.getAlgorithms()) {
            if (!alg.equals(GATKSVVCFConstants.DEPTH_ALGORITHM)) return false;
        }
        return true;
    }

    @VisibleForTesting
    public static double getDefaultPaddingFraction() {
        return DEFAULT_PADDING_FRACTION;
    }
}
