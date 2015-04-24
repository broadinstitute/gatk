package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.StringUtil;

/**
 * A read name encoder conforming to the standard described by Illumina Casava 1.8.
 * 
 * @see <a href="http://biowulf.nih.gov/apps/CASAVA1_8_Changes.pdf">Casava 1.8 update</a>
 * @author mccowan
 */
public final class Casava18ReadNameEncoder implements ReadNameEncoder {
    final static int CONTROL_FIELD_VALUE = 0;
    final String runId, instrumentName, flowcellId;
    
    static enum IsFilteredLabel {
        Y, N;
        static IsFilteredLabel get(final boolean passesFilter) {
            return passesFilter ? N : Y;
        }
    }
    
    public Casava18ReadNameEncoder(final String instrumentName, final String runId, final String flowcellId) {
        this.runId = runId;
        this.instrumentName = instrumentName;
        this.flowcellId = flowcellId;
    }

    @Override
    public String generateReadName(final ClusterData cluster, final Integer pairNumber) {
        return String.format(
                "%s:%s:%s:%d:%d:%d:%d %s:%s:%d:%s",
                instrumentName,
                runId,
                flowcellId,
                cluster.getLane(),
                cluster.getTile(),
                cluster.getX(),
                cluster.getY(),
                StringUtil.asEmptyIfNull(pairNumber),
                IsFilteredLabel.get(cluster.isPf()),
                CONTROL_FIELD_VALUE,
                StringUtil.asEmptyIfNull(cluster.getMatchedBarcode())
        );
    }
}
