package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Utility class used to create {@link TargetCollection} from common target data sources.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCollectionUtils {

    /**
     * Disabled instantiation of this helper class.
     */
    private TargetCollectionUtils() {

    }

    /**
     * Creates a target collection out from an interval list where each element corresponds to a target.
     *
     * @param intervals the target intervals.
     * @throws IllegalArgumentException if {@code intervals} is {@code null}, it contains the {@code null} element or
     *   overlapping intervals.
     * @return never {@code null}.
     */
    public static TargetCollection<SimpleInterval> fromSimpleIntervalList(final List<SimpleInterval> intervals) {
        return new HashedListTargetCollection<SimpleInterval>(intervals) {

            @Override
            public String name(final SimpleInterval target) {
                return Utils.nonNull(target,"the input target cannot be null").toString();
            }

            @Override
            public SimpleInterval location(final SimpleInterval target) {
                return Utils.nonNull(target,"the input target cannot be null");
            }
        };
    }

     /**
     * Creates a target collection out from a BED feature list.
     * <p>
     *     Each feature in the input bed feature list will correspond to a target whose name is the name of the BED
     *     feature as returned by {@link BEDFeature#getName()}.
     * </p>
     *
     * @param features the target features.
     * @param <T> the specific {@link BEDFeature} subclass.
     * @throws IllegalArgumentException if {@code features} is {@code null}, it contains the {@code null} element or
     *     overlapping features.
     *
     * @return never {@code null}.
     */
    public static <T extends BEDFeature> TargetCollection<T> fromBEDFeatureList(final List<T> features) {
        return new HashedListTargetCollection<T>(Utils.nonNull(features,"the input feature list cannot be null")) {
            @Override
            public String name(final T target) {
                final String name = Utils.nonNull(target,"the input target cannot be null").getName();
                return (name == null || name.isEmpty()) ? null : name;
            }

            @Override
            public SimpleInterval location(final T target) {
                return new SimpleInterval(Utils.nonNull(target,"the input target cannot be null"));
            }
        };
    }

    /**
     * Creates a target collection out from a BED feature containing file.
     *
     * @param file the target feature file.
     * @param <T> the specific {@link BEDFeature} subclass.
     * @throws IllegalArgumentException if {@code features} is {@code null}, it contains the {@code null} element or
     *     overlapping features.
     *
     * @return never {@code null}.
     */
    public static <T extends BEDFeature> TargetCollection<T> fromBEDFeatureFile(final File file, final FeatureCodec<T, ?> codec) {
        Utils.nonNull(file,"the input file cannot be null");
        Utils.nonNull(codec,"the codec cannot be null");
        final FeatureDataSource<T> source = new FeatureDataSource<>(file,codec);
        final List<T> features = StreamSupport.stream(source.spliterator(),false).collect(Collectors.toList());
        return fromBEDFeatureList(features);
    }
}
