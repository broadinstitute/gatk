package org.broadinstitute.hellbender.utils.test;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class FuncotatorTestUtils {
    private FuncotatorTestUtils() {}


    /**
     * Since funcotation factories need an instance of {@link FeatureContext} to funcotate, this convenience method can
     *  create a new instance for test methods.
     *
     * @param funcotationFactories {@link List} of {@link DataSourceFuncotationFactory} that should be used to generate the
     *                                         {@link FeatureContext}.  Never {@code null}, but empty list is acceptable.
     * @param dummyToolInstanceName A name to use for the "tool".  Any string will work here.  Never {@code null}.
     * @param interval genomic interval for the result.  Typically, this would be the interval of the variant.  Never {@link null}.
     * @param featureQueryLookahead When querying FeatureDataSources, cache this many extra bases of context beyond
     *                              the end of query intervals in anticipation of future queries. Must be >= 0.  If uncertain, use zero.
     * @param cloudPrefetchBuffer See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, Path)}  If uncertain, use zero.
     * @param cloudIndexPrefetchBuffer See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, Path)}  If uncertain, use zero.
     * @param reference See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, Path)}  If uncertain, use {@code null}.
     * @return a {@link FeatureContext} ready for querying the funcotation factories on the given interval.  Never {@code null}.
     */
    public static FeatureContext createFeatureContext(final List<DataSourceFuncotationFactory> funcotationFactories, final String dummyToolInstanceName,
                                                      final SimpleInterval interval, final int featureQueryLookahead, final int cloudPrefetchBuffer,
                                                      final int cloudIndexPrefetchBuffer, final Path reference) {
        Utils.nonNull(funcotationFactories);
        Utils.nonNull(dummyToolInstanceName);
        Utils.nonNull(interval);

        final Map<FeatureInput<? extends Feature>, Class<? extends Feature>> featureInputsWithType =
                funcotationFactories.stream()
                        .collect(Collectors.toMap(ff -> ff.getMainSourceFileAsFeatureInput(), ff -> ff.getAnnotationFeatureClass()));

        return FeatureContext.create(featureInputsWithType, dummyToolInstanceName, interval,
                featureQueryLookahead, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, reference);
    }
}
