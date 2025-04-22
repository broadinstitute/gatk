package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.walkers.sv.SVStratify;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Handles stratification of variants for genotype concordance analysis of SVs, similar to the functionality
 * implemented in {@link org.broadinstitute.hellbender.tools.walkers.sv.GroupedSVCluster}.
 *
 * Truth and eval variants can be added with their respective methods and are stratified using the stratification
 * machinery passed in through this class's constructor. Concordance results are buffered in the {@link ItemTracker}
 * subclass, which keeps track of which stratifications each eval variant has been submitted to and returned. The
 * stratification group with the highest priority, as determined by the ordering returned by the
 * {@link SVStratificationEngine#getStrata()} method,  and matched truth variant is returned to the final output buffer
 * once a variant has been returned from all stratification groups. The "default" group is always added and of lowest
 * priority. If an eval variant fails to match in all groups, then there is no guarantee on which group's output variant
 * is returned, but they are all identical (i.e. a copy of the input variant with a "FP" status field).
 *
 * Output records must be retrieved using the {@link #flush} method and are not sorted.
 *
 * TODO: This implementation will result in proliferation of deep copies of the input variant, with up to one copy per
 *   possible stratification group. Avoiding the deep copies, particularly of genotypes, would substantially improve
 *   performance on large inputs.
 *
 */
public class StratifiedConcordanceEngine {

    protected final Map<String, ClosestSVFinder> clusterEngineMap;
    protected final SVStratificationEngine stratificationEngine;
    protected final Map<Long, ItemTracker> variantStatusMap = new HashMap<>();
    protected List<VariantContext> outputBuffer;
    protected ClosestSVFinder defaultEngine;
    protected SAMSequenceDictionary dictionary;
    protected final SVStratificationEngineArgumentsCollection stratArgs;
    protected Long nextItemId = 0L;

    public StratifiedConcordanceEngine(final Map<String, ClosestSVFinder> clusterEngineMap,
                                       final SVStratificationEngine stratificationEngine,
                                       final SVStratificationEngineArgumentsCollection stratArgs,
                                       final SVClusterEngineArgumentsCollection defaultClusteringArgs,
                                       final SVConcordanceAnnotator defaultCollapser,
                                       final SAMSequenceDictionary dictionary) {
        Utils.validate(stratificationEngine.getStrata().size() == clusterEngineMap.size(),
                "Stratification and clustering configurations have a different number of groups.");
        for (final SVStratificationEngine.Stratum stratum : stratificationEngine.getStrata()) {
            Utils.validate(clusterEngineMap.containsKey(stratum.getName()),
                    "Could not find group " + stratum.getName() + " in clustering configuration.");
        }

        this.clusterEngineMap = clusterEngineMap;
        this.stratificationEngine = stratificationEngine;
        this.stratArgs = stratArgs;
        this.outputBuffer = new ArrayList<>();
        this.dictionary = dictionary;

        final SVConcordanceLinkage defaultLinkage = new SVConcordanceLinkage(dictionary);
        defaultLinkage.setDepthOnlyParams(defaultClusteringArgs.getDepthParameters());
        defaultLinkage.setMixedParams(defaultClusteringArgs.getMixedParameters());
        defaultLinkage.setEvidenceParams(defaultClusteringArgs.getPESRParameters());
        this.defaultEngine = new ClosestSVFinder(defaultLinkage, defaultCollapser::annotate, false, dictionary);
    }

    /**
     * Retrieve finished eval variants with concordance annotations
     * @param force returns all variants regardless of the position of the last added variant
     * @return final output records
     */
    public Collection<VariantContext> flush(boolean force) {
        for (final String name : clusterEngineMap.keySet()) {
            flushEngineToBuffer(clusterEngineMap.get(name), force, name);
        }
        flushEngineToBuffer(defaultEngine, force, SVStratify.DEFAULT_STRATUM);
        final Collection<VariantContext> result = outputBuffer;
        outputBuffer = new ArrayList<>();
        return result;
    }

    /**
     * Checks if there are any active eval variants or unflushed output variants
     */
    public boolean isEmpty() {
        return variantStatusMap.isEmpty() && outputBuffer.isEmpty();
    }

    /**
     * Adds a "truth" variant
     */
    public void addTruthVariant(final SVCallRecord record) {
        final Long id = nextItemId++;
        // truth variants can be matched across all stratification groups
        for (final String name : clusterEngineMap.keySet()) {
            addToEngine(record, id, true, clusterEngineMap.get(name), name);
        }
        addToEngine(record, id, true, defaultEngine, SVStratify.DEFAULT_STRATUM);
    }

    /**
     * Adds an "eval" variant
     */
    public void addEvalVariant(final SVCallRecord record) {
        final Long id = nextItemId++;
        // eval variants get stratified
        final List<String> stratifications = stratificationEngine.getMatches(record,
                        stratArgs.overlapFraction, stratArgs.numBreakpointOverlaps, stratArgs.numBreakpointOverlapsInterchrom)
                .stream().map(SVStratificationEngine.Stratum::getName).collect(Collectors.toUnmodifiableList());
        for (final String stratum : stratifications) {
            Utils.validate(clusterEngineMap.containsKey(stratum), "Group undefined: " + stratum);
            addToEngine(record, id, false, clusterEngineMap.get(stratum), stratum);
        }
        // default stratum
        addToEngine(record, id, false, defaultEngine, SVStratify.DEFAULT_STRATUM);
        // Keep track of which groups it was added to
        Utils.validate(!variantStatusMap.containsKey(id), "Attempted to add duplicate id");
        final List<String> stratList = new ArrayList<>(stratifications.size() + 1);
        stratList.addAll(stratifications);
        stratList.add(SVStratify.DEFAULT_STRATUM);
        variantStatusMap.put(id, new ItemTracker(stratList));
    }

    /**
     * Adds a record to the given concordance engine, flushing active variants to the buffer if it hit a new contig
     */
    protected void addToEngine(final SVCallRecord record, final Long id, final boolean isTruth, final ClosestSVFinder engine, final String name) {
        if (engine.getLastItemContig() != null && !record.getContigA().equals(engine.getLastItemContig())) {
            flushEngineToBuffer(engine, true, name);
        }
        engine.add(record, id, isTruth);
        flushEngineToBuffer(engine, false, name);
    }

    /**
     * Flushes active eval variants in a specific group's engine to the output buffer
     */
    protected void flushEngineToBuffer(final ClosestSVFinder engine, final boolean force, final String name) {
        for (final ClosestSVFinder.LinkageConcordanceRecord record: engine.flush(force)) {
            final ItemTracker tracker = variantStatusMap.get(record.id());
            Utils.validate(tracker != null, "Unregistered variant id: " + record.id());
            tracker.eject(name, record.record());
            if (tracker.allEjected()) {
                Utils.nonNull(tracker.getOutput(), "No tracker output");
                variantStatusMap.remove(record.id());
                final List<String> strata = tracker.getGroups().stream().sorted().collect(Collectors.toUnmodifiableList());
                final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(tracker.getOutput())
                        .attribute(GATKSVVCFConstants.STRATUM_INFO_KEY, strata);
                outputBuffer.add(builder.make());
            }
        }
    }

    /**
     * Tracks a given eval variant's stratification groups and its status in each. Prioritizes the final output variant
     * as the one produced by the highest priority group that also matched a truth variant.
     */
    protected static class ItemTracker {

        private final Map<String, Integer> groupPriorityMap;
        private SVCallRecord out;
        private Integer outPriority;
        private final List<String> groups;

        /**
         * Constructor with stratification groups to which this record belongs. Note that the input record is not
         * tracked because only the output record(s) are needed.
         * @param submittedTo stratification group names, sorted in order of decreasing priority
         */
        public ItemTracker(final List<String> submittedTo) {
            this.groups = submittedTo;
            this.groupPriorityMap = new HashMap<>(SVUtils.hashMapCapacity(submittedTo.size()));
            for (int i = 0; i < submittedTo.size(); i++) {
                this.groupPriorityMap.put(submittedTo.get(i), i);
            }
        }

        /**
         * To be called when a record is returned from one of the engines. May only be called once per group.
         * @param name name of the group
         * @param record record returned by that group's engine
         */
        public void eject(final String name, final SVCallRecord record) {
            final Integer priority = groupPriorityMap.remove(name);
            Utils.nonNull(priority, "Unregistered name: " + name);
            // assign to output if higher than current output's priority
            if (outPriority == null || (isTruePositive(record)) && (priority < outPriority || (priority > outPriority && !isTruePositive(out)))) {
                out = record;
                outPriority = priority;
            }
        }

        private boolean isTruePositive(final SVCallRecord record) {
            return ConcordanceState.TRUE_POSITIVE.getAbbreviation().equals(record.getAttributes().get(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE));
        }

        /**
         * Returns true if this tracker is done
         */
        public boolean allEjected() {
            return groupPriorityMap.isEmpty();
        }

        /**
         * Get the highest priority output record. Should only be called if {@link #allEjected()} is true.
         */
        public SVCallRecord getOutput() {
            return out;
        }

        public List<String> getGroups() {
            return groups;
        }
    }
}
