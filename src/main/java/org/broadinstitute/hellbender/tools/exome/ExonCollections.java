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
 * Utility class used to create {@link ExonCollections} from common exon data sources.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ExonCollections {

    /**
     * Disabled instantiation of this helper class.
     */
    private ExonCollections() {

    }

    /**
     * Creates an exon collection out from an interval list where each element corresponds to an exon.
     *
     * @param intervals the exon intervals.
     * @throws IllegalArgumentException if {@code intervals} is {@code null}, it contains the {@code null} element or
     *   overlapping intervals.
     * @return never {@code null}.
     */
    public static ExonCollection<SimpleInterval> fromSimpleIntervalList(final List<SimpleInterval> intervals) {
        return new HashedListExonCollection<SimpleInterval>(intervals) {

            @Override
            public String name(final SimpleInterval exon) {
                return Utils.nonNull(exon,"the input exon cannot be null").toString();
            }

            @Override
            public SimpleInterval location(final SimpleInterval exon) {
                return Utils.nonNull(exon,"the input exon cannot be null");
            }
        };
    }

     /**
     * Creates an exon collection out from a BED feature list.
     * <p>
     *     Each feature in the input bed feature list will correspond to an exon whose name is the name of the BED
     *     feature as returned by {@link BEDFeature#getName()}.
     * </p>
     *
     * @param features the exon features.
     * @param <T> the specific {@link BEDFeature} subclass.
     * @throws IllegalArgumentException if {@code features} is {@code null}, it contains the {@code null} element or
     *     overlapping features.
     *
     * @return never {@code null}.
     */
    public static <T extends BEDFeature> ExonCollection<T> fromBEDFeatureList(final List<T> features) {
        return new HashedListExonCollection<T>(Utils.nonNull(features,"the input feature list cannot be null")) {
            @Override
            public String name(final T exon) {
                final String name = Utils.nonNull(exon,"the input exon cannot be null").getName();
                return (name == null || name.isEmpty()) ? null : name;
            }

            @Override
            public SimpleInterval location(final T exon) {
                return new SimpleInterval(Utils.nonNull(exon,"the input exon cannot be null"));
            }
        };
    }

    /**
     * Creates an exon collection out from a BED feature containing file.
     *
     * @param file the exon feature file.
     * @param <T> the specific {@link BEDFeature} subclass.
     * @throws IllegalArgumentException if {@code features} is {@code null}, it contains the {@code null} element or
     *     overlapping features.
     *
     * @return never {@code null}.
     */
    public static <T extends BEDFeature> ExonCollection<T> fromBEDFeatureFile(final File file, final FeatureCodec<T, ?> codec) {
        Utils.nonNull(file,"the input file cannot be null");
        Utils.nonNull(codec,"the codec cannot be null");
        final FeatureDataSource<T> source = new FeatureDataSource<>(file,codec);
        final List<T> features = StreamSupport.stream(source.spliterator(),false).collect(Collectors.toList());
        return fromBEDFeatureList(features);
    }
}
