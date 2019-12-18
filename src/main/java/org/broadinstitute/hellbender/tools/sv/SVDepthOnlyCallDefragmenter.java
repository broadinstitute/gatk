package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SVDepthOnlyCallDefragmenter extends LocatableClusterEngine<SVCallRecord> {

    private final double minSampleOverlap;
    private static final double PADDING_FRACTION = 0.5;

    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary) {
        this(dictionary, 0.9, null);
    }

    //for single-sample clustering case
    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary, double minSampleOverlap, List<GenomeLoc> coverageIntervals) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE, coverageIntervals);
        this.minSampleOverlap = minSampleOverlap;
    }

    /**
     * Find a single call representative of all the calls in the {@param cluster}
     * @param cluster   the events that are clustered together
     * @return  a call encompassing all the cluster's events and containing all the algorithms and genotypes
     */
    @Override
    protected SVCallRecord flattenCluster(final Collection<SVCallRecord> cluster) {
        final int newStart =  cluster.stream().mapToInt(SVCallRecord::getPositionA).min().getAsInt();
        final int newEnd = cluster.stream().mapToInt(SVCallRecord::getPositionB).max().getAsInt();
        final SVCallRecord exampleCall = cluster.iterator().next();
        final int length = newEnd - newStart + 1;  //+1 because GATK intervals are inclusive
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList()); //should be depth only
        final List<Genotype> clusterGenotypes = deduplicateGenotypes(cluster.stream().flatMap(v -> v.getGenotypes().stream()).collect(Collectors.toList()));
        return new SVCallRecord(exampleCall.getId(), exampleCall.getContigA(), newStart, exampleCall.getStrandA(),
                exampleCall.getContigB(), newEnd, exampleCall.getStrandB(), exampleCall.getType(), length, algorithms, clusterGenotypes);
    }

    @Override
    protected SVDeduplicator<SVCallRecord> getDeduplicator() {
        final Function<Collection<SVCallRecord>,SVCallRecord> collapser = SVCallRecordUtils::deduplicateWithRawCallAttribute;
        return new SVCallRecordDeduplicator<>(collapser, dictionary);
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
        Utils.nonEmpty(genotypesForSameSample);
        final String sampleName = genotypesForSameSample.get(0).getSampleName();
        if (!genotypesForSameSample.stream().allMatch(g -> g.getSampleName().equals(sampleName))) {
            throw new IllegalArgumentException("This method expects a list of genotypes from the same sample, " +
                    "but not all input genotypes represent sample " + sampleName + ".");
        }
        final Set<Integer> copyNumbers = genotypesForSameSample.stream()
                .filter(g -> g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))
                .map(g -> Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()))
                .collect(Collectors.toSet());
        final GenotypeBuilder gb = new GenotypeBuilder(genotypesForSameSample.get(0));
        if (copyNumbers.isEmpty()) {
            return gb.make();
        } else if (copyNumbers.size() == 1) {
            //For now just make sure genotypes have the same copy number -- qualities will be recalculated elsewhere
            gb.noAttributes();
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumbers.iterator().next());
            return gb.make();
        } else {
            throw new IllegalArgumentException("This method will only merge genotypes with the same copy number. Found " + copyNumbers.size() + " different copy numbers.");
        }
    }

    /**
     * Determine if two calls should cluster based on their padded intervals and genotyped samples
     * @param a
     * @param b
     * @return true if the two calls should be in the same cluster
     */
    @Override
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        if (!isDepthOnlyCall(a) || !isDepthOnlyCall(b)) return false;
        Utils.validate(a.getContigA().equals(a.getContigB()), "Call A is depth-only but interchromosomal");
        Utils.validate(b.getContigA().equals(b.getContigB()), "Call B is depth-only but interchromosomal");
        if (!a.getType().equals(b.getType())) return false;
        final Set<String> sharedSamples = new LinkedHashSet<>(a.getCalledSamples());
        sharedSamples.retainAll(b.getCalledSamples());
        final double sampleOverlap = Math.min(sharedSamples.size() / (double) a.getCalledSamples().size(), sharedSamples.size() / (double) b.getCalledSamples().size());
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
     * Determine an overlap interval for clustering using {@value #PADDING_FRACTION} padding
     * Returned interval represents the interval in which the start position of a new event must fall in order to be added to the cluster (including {@param call})
     * @param call  new event to be clustered
     * @param currentClusterInterval    the cluster of interest, may be null
     * @return  an interval describing the cluster after {@param call} is added
     */
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecord call, final SimpleInterval currentClusterInterval) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getCallInterval(call);
        final int paddedCallStart, paddedCallEnd;
        if (genomicToBinMap != null) {
            final GenomeLoc callStart = parser.createGenomeLoc(call.getContigA(), call.getPositionA(), call.getPositionA());
            final GenomeLoc callEnd = parser.createGenomeLoc(call.getContigA(), call.getPositionB(), call.getPositionB());
            //first interval that is equal to or "greater than" the call start, such that the start of the bin should match the call start, with a little wiggle room
            final Map.Entry<GenomeLoc, Integer> startBin = genomicToBinMap.ceilingEntry(callStart);
            if (startBin == null) {
                throw new UserException.BadInput("Call start " + callStart + " for  call " + call.getId() + " not found in model call intervals.");
            }
            final int callStartIndex = startBin.getValue();
            //last interval that is equal to or "less than" the call start, such that the end of the bin should match the call end
            final Map.Entry<GenomeLoc, Integer> endBin = genomicToBinMap.floorEntry(callEnd);
            if (endBin == null) {
                throw new UserException.BadInput("Call end " + callEnd + " for call " + call.getId() + " not found in model call intervals.");
            }
            final int callEndIndex = endBin.getValue();
            final int callBinLength = callEndIndex - callStartIndex + 1;
            if (callBinLength <= 0) {
                throw new UserException.BadInput("Copy number call at " + call.getContigA() + ":" + call.getPositionA() + "-"
                        + call.getPositionB() + " does not align with supplied model calling intervals. Use the filtered intervals input from GermlineCNVCaller for this cohort/model.");
            }
            final int paddedStartIndex = Math.max(callStartIndex - (int)Math.round(callBinLength * PADDING_FRACTION), 0);
            if (coverageIntervals.get(paddedStartIndex).getContig().equals(callStart.getContig())) {
                paddedCallStart = coverageIntervals.get(paddedStartIndex).getStart();
            } else {
                paddedCallStart = callStart.getStart();
            }
            final int paddedEndIndex = Math.min(callEndIndex + (int)Math.round(callBinLength * PADDING_FRACTION), genomicToBinMap.size() - 1);
            if (coverageIntervals.get(paddedEndIndex).getContig().equals(callEnd.getContig())) {
                paddedCallEnd = coverageIntervals.get(paddedEndIndex).getEnd();
            } else {
                paddedCallEnd = callEnd.getEnd();
            }
        } else {
            paddedCallStart = (int) (callInterval.getStart() - PADDING_FRACTION * callInterval.getLengthOnReference());
            paddedCallEnd = (int) (callInterval.getEnd() + PADDING_FRACTION * callInterval.getLengthOnReference());
        }
        final int contigLength = dictionary.getSequence(call.getContigA()).getSequenceLength();
        if (currentClusterInterval == null) {
            return IntervalUtils.trimIntervalToContig(call.getContigA(), paddedCallStart, paddedCallEnd, contigLength);
        }
        //NOTE: this is an approximation -- padding should be based on the length of the call plus currentClusterIntervals
        final int newMinStart = Math.min(paddedCallStart, currentClusterInterval.getStart());
        final int newMaxEnd = Math.max(paddedCallEnd, currentClusterInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(call.getContigA(), newMinStart, newMaxEnd, contigLength);
    }

    private SimpleInterval getCallInterval(final SVCallRecord call) {
        return new SimpleInterval(call.getContigA(), call.getPositionA(), call.getPositionB());
    }

    public static boolean isDepthOnlyCall(final SVCallRecord call) {
        if (call.getAlgorithms().isEmpty()) return false;
        for (final String alg : call.getAlgorithms()) {
            if (!alg.equals(GATKSVVCFConstants.DEPTH_ALGORITHM)) return false;
        }
        return true;
    }

    @VisibleForTesting
    public static double getPaddingFraction() {
        return PADDING_FRACTION;
    }
}
