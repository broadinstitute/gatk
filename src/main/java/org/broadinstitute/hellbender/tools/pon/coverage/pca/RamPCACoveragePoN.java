package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A PoN class that is stored in RAM.
 */
public final class RamPCACoveragePoN implements PCACoveragePoN {
    private final double version;

    private final List<Target> targets;
    private final List<Target> rawTargets;
    private final List<Target> panelTargets;

    private final double[] targetFactors;
    private final double[] targetVariances;

    private final RealMatrix normalizedCounts;
    private final RealMatrix logNormalizedCounts;
    private final RealMatrix logNormalizedPInverseCounts;
    private final RealMatrix reducedPanelCounts;
    private final RealMatrix reducedPanelPInverseCounts;

    private final List<String> targetNames;
    private final List<String> rawTargetNames;
    private final List<String> panelTargetNames;

    private final List<String> sampleNames;
    private final List<String> panelSampleNames;

    /**
     * Instantiates a RamPoN that is a duplicate of the given PoN
     * @param inputPoN Never {@code null}
     */
    public RamPCACoveragePoN(final PCACoveragePoN inputPoN) {
        //make defensive copies of lists
        this.version = inputPoN.getVersion();
        this.targets = new ArrayList<>(inputPoN.getTargets());
        this.rawTargets = new ArrayList<>(inputPoN.getRawTargets());
        this.panelTargets = new ArrayList<>(inputPoN.getPanelTargets());
        this.targetFactors = Arrays.copyOf(inputPoN.getTargetFactors(), inputPoN.getTargetFactors().length);
        this.targetVariances = Arrays.copyOf(inputPoN.getTargetVariances(), inputPoN.getTargetVariances().length);
        this.normalizedCounts = inputPoN.getNormalizedCounts();
        this.logNormalizedCounts = inputPoN.getLogNormalizedCounts();
        this.logNormalizedPInverseCounts = inputPoN.getLogNormalizedPInverseCounts();
        this.reducedPanelCounts = inputPoN.getReducedPanelCounts();
        this.reducedPanelPInverseCounts = inputPoN.getReducedPanelPInverseCounts();
        this.targetNames = new ArrayList<>(inputPoN.getTargetNames());
        this.rawTargetNames = new ArrayList<>(inputPoN.getRawTargetNames());
        this.panelTargetNames = new ArrayList<>(inputPoN.getPanelTargetNames());
        this.sampleNames = new ArrayList<>(inputPoN.getSampleNames());
        this.panelSampleNames = new ArrayList<>(inputPoN.getPanelSampleNames());
    }

    @Override
    public double getVersion() {
        return version;
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
    public double[] getTargetFactors() {
        return targetFactors;
    }

    @Override
    public double[] getTargetVariances() {
        return targetVariances;
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

    @Override
    public List<String> getTargetNames() {
        return Collections.unmodifiableList(targetNames);
    }

    @Override
    public List<String> getRawTargetNames() {
        return Collections.unmodifiableList(rawTargetNames);
    }

    @Override
    public List<String> getPanelTargetNames() {
        return Collections.unmodifiableList(panelTargetNames);
    }

    @Override
    public List<String> getSampleNames() {
        return Collections.unmodifiableList(sampleNames);
    }

    @Override
    public List<String> getPanelSampleNames() {
        return Collections.unmodifiableList(panelSampleNames);
    }
}
