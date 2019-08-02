package org.broadinstitute.hellbender.engine;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;

/**
 * A Dummy / Placeholder class that can be used where a {@link GATKTool} is required.
 * DO NOT USE THIS FOR ANYTHING OTHER THAN TESTING.
 * THIS MUST BE IN THE ENGINE PACKAGE DUE TO SCOPE ON `features`!
 * Created by jonn on 9/19/18.
 */
@CommandLineProgramProperties(
        summary = "A dummy GATKTool to help test Funcotator.",
        oneLineSummary = "Dummy dumb dumb tool for testing.",
        programGroup = TestProgramGroup.class
)
public final class DummyPlaceholderGatkTool extends GATKTool {

    public DummyPlaceholderGatkTool() {
        parseArgs(new String[]{});
        onStartup();
    }

    @Override
    public void traverse() {

    }

    @Override
    void initializeFeatures(){
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                getGenomicsDBOptions());
    }
}
