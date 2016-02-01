package org.broadinstitute.hellbender.utils.hdf5;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.exome.ReductionResult;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A PoN class that is stored in RAM.
 */
public class RamPoN implements PoN {

    private final List<String> targetNames;
    private final List<String> rawTargetNames;

    private final List<String> sampleNames;

    private final double version;

    private final List<String> panelTargetNames;
    private final List<String> panelSampleNames;

    private final List<Target> targets;
    private final List<Target> rawTargets;
    private final List<Target> panelTargets;

    private RealMatrix targetFactors;
    private final RealMatrix normalizedCounts;
    private final RealMatrix logNormalizedCounts;
    private final RealMatrix logNormalizedPInverseCounts;
    private final RealMatrix reducedPanelCounts;
    private final RealMatrix reducedPanelPInverseCounts;

    private final double[] targetVariances;

    /**
     * Instantiates a RamPoN that is a duplicate of the given PoN
     * @param inputPoN Never {@code null}
     */
    public RamPoN(final PoN inputPoN) {
        this.targetNames = inputPoN.getTargetNames();
        this.sampleNames = inputPoN.getSampleNames();
        this.version = inputPoN.getVersion();
        this.panelTargetNames = inputPoN.getPanelTargetNames();
        this.panelSampleNames = inputPoN.getPanelSampleNames();
        this.normalizedCounts = inputPoN.getNormalizedCounts();
        this.logNormalizedCounts = inputPoN.getLogNormalizedCounts();
        this.logNormalizedPInverseCounts = inputPoN.getLogNormalizedPInverseCounts();
        this.reducedPanelCounts = inputPoN.getReducedPanelCounts();
        this.reducedPanelPInverseCounts = inputPoN.getReducedPanelPInverseCounts();
        this.targets = new ArrayList<>(inputPoN.getTargets());
        this.rawTargets = new ArrayList<>(inputPoN.getRawTargets());
        this.rawTargetNames = new ArrayList<>(inputPoN.getRawTargetNames());
        this.panelTargets = new ArrayList<>(inputPoN.getPanelTargets());
        this.targetFactors = inputPoN.getTargetFactors();
        this.targetVariances = inputPoN.getTargetVariances();
    }

    /**
     *  For use when overriding values of an input PoN with a separate reduction result.
     *
     * @param inputPoN Never {@code null}
     * @param ponReduction Never {@code null}
     */
    public RamPoN(final PoN inputPoN, final ReductionResult ponReduction) {
        this.targetNames = inputPoN.getTargetNames();
        this.sampleNames = inputPoN.getSampleNames();
        this.version = inputPoN.getVersion();
        this.panelTargetNames = inputPoN.getPanelTargetNames();
        this.panelSampleNames = inputPoN.getPanelSampleNames();
        this.normalizedCounts = inputPoN.getNormalizedCounts();
        this.logNormalizedCounts = inputPoN.getLogNormalizedCounts();
        this.logNormalizedPInverseCounts = inputPoN.getLogNormalizedPInverseCounts();
        this.targets = new ArrayList<>(inputPoN.getTargets());
        this.rawTargets = new ArrayList<>(inputPoN.getRawTargets());
        this.rawTargetNames = new ArrayList<>(inputPoN.getRawTargetNames());
        this.targetFactors = inputPoN.getTargetFactors();

        // Note that if the PoN reduction step ever filters targets, this next line will no longer be correct.
        this.panelTargets = new ArrayList<>(inputPoN.getPanelTargets());

        // Use the reduction result for the reduction fields.
        this.reducedPanelCounts = ponReduction.getReducedCounts();
        this.reducedPanelPInverseCounts = ponReduction.getReducedInverse();

        this.targetVariances = inputPoN.getTargetVariances();
    }

    @Override
    public List<String> getTargetNames() {
        return Collections.unmodifiableList(new ArrayList<>(targetNames));
    }

    @Override
    public List<String> getRawTargetNames() {
        return Collections.unmodifiableList(new ArrayList<>(rawTargetNames));
    }

    @Override
    public List<String> getSampleNames() {
        return Collections.unmodifiableList(new ArrayList<>(sampleNames));
    }

    @Override
    public double getVersion() {
        return version;
    }

    @Override
    public List<String> getPanelTargetNames() {
        return Collections.unmodifiableList(new ArrayList<>(panelTargetNames));
    }

    @Override
    public List<String> getPanelSampleNames() {
        return Collections.unmodifiableList(new ArrayList<>(panelSampleNames));
    }

    @Override
    public RealMatrix getTargetFactors() {
        return targetFactors.copy();
    }

    @Override
    public List<Target> getTargets() {
        return targets.stream().map(t -> new Target(t.getName(), new SimpleInterval(t.getContig(), t.getStart(), t.getEnd()))).collect(Collectors.toList());
    }

    @Override
    public List<Target> getRawTargets() {
        return rawTargets.stream().map(t -> new Target(t.getName(), new SimpleInterval(t.getContig(), t.getStart(), t.getEnd()))).collect(Collectors.toList());
    }

    @Override
    public List<Target> getPanelTargets() {
        return panelTargets.stream().map(t -> new Target(t.getName(), new SimpleInterval(t.getContig(), t.getStart(), t.getEnd()))).collect(Collectors.toList());
    }

    @Override
    public double[] getTargetVariances() {
        return Arrays.copyOf(targetVariances, targetVariances.length);
    }

    @Override
    public void setTargetFactors(final RealMatrix targetFactors) {
        Utils.nonNull(targetFactors);
        if (targetFactors.getColumnDimension() != 1) {
            throw new IllegalArgumentException("the number of columns in the target factors matrix must be 1 but it is: " + targetFactors.getColumnDimension());
        }
        this.targetFactors = targetFactors;
    }

    @Override
    public RealMatrix getNormalizedCounts() {
        return normalizedCounts.copy();
    }

    @Override
    public RealMatrix getLogNormalizedCounts() {
        return logNormalizedCounts.copy();
    }

    @Override
    public RealMatrix getLogNormalizedPInverseCounts() {
        return logNormalizedPInverseCounts.copy();
    }

    @Override
    public RealMatrix getReducedPanelCounts() {
        return reducedPanelCounts.copy();
    }

    @Override
    public RealMatrix getReducedPanelPInverseCounts() {
        return reducedPanelPInverseCounts.copy();
    }
}
