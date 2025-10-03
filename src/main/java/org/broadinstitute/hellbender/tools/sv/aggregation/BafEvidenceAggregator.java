package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class BafEvidenceAggregator extends SVEvidenceAggregator<BafEvidence> {

    private final double paddingFraction;

    /**
     * Constructor.
     * @param source        PE feature source
     * @param dictionary    sequence dictionary
     * @param paddingFraction  evidence window as a fraction of variant size
     */
    public BafEvidenceAggregator(final FeatureDataSource<BafEvidence> source,
                                 final SAMSequenceDictionary dictionary,
                                 final double paddingFraction) {
        super(source, dictionary); // there is no fixed upper bound on the window since it's a fraction of length
        this.paddingFraction = paddingFraction;
    }

    @Override
    public SimpleInterval getEvidenceQueryInterval(final SVCallRecord call) {
        final int padding = (int) Math.ceil(call.getLength() * paddingFraction);
        final SimpleInterval result = new SimpleInterval(call.getContigA(), call.getPositionA(), call.getPositionB())
                .expandWithinContig(padding, dictionary);
        Utils.nonNull(result, "Error generating padded interval for variant " + call.getId() + "; check that its coordinates are valid");
        return result;
    }

    public double getPaddingFraction() {
        return paddingFraction;
    }

    @Override
    public boolean evidenceFilter(final SVCallRecord record, final BafEvidence evidence) {
        // No additional properties to match on, other than assumed overlap
        return true;
    }
}
