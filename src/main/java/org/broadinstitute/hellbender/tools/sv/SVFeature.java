package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.PriorityQueue;
import java.util.Set;

public interface SVFeature extends Feature {
    SVFeature extractSamples( final Set<String> sampleNames, final Object header );
}
