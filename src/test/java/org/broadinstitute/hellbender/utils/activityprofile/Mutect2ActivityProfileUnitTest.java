package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState.Type.*;

public class Mutect2ActivityProfileUnitTest extends GATKBaseTest {
    private static final String CONTIG = "1";
    private static final int START = 1001;
    private static final SAMFileHeader HEADER = new SAMFileHeader(new SAMSequenceDictionary(List.of(new SAMSequenceRecord(CONTIG, 100_000))));

    // default max assembly region size is 300 and default max phasing distance is 50
    private static final int DEFAULT_PHASING_DISTANCE = new AssemblyRegionArgumentCollection().maxProbPropagationDistance;

    /**
     * Active states (probability 1.0) of given types at given offsets from START.  All other states are inactive.
     */
    private static final class Sites {
        final Map<Integer, ActivityProfileState.Type> types = new HashMap<>();

        Sites add(final ActivityProfileState.Type type, final int... offsets) {
            Arrays.stream(offsets).forEach(offset -> types.put(offset, type));
            return this;
        }
    }

    /**
     * Add states at offsets 0 through numStates - 1 from START, pop regions, and describe each as "firstOffset-lastOffset:isActive"
     */
    private static List<String> popRegions(final int numStates, final Sites sites, final double germlineFraction,
                                           final int maxPhasingDistance, final boolean forceConversion) {
        final AssemblyRegionArgumentCollection args = new AssemblyRegionArgumentCollection();
        args.genotypeGermlineSitesFraction = germlineFraction;
        args.maxProbPropagationDistance = maxPhasingDistance;
        final Mutect2ActivityProfile profile = new Mutect2ActivityProfile(args, HEADER);

        for (int offset = 0; offset < numStates; offset++) {
            final ActivityProfileState.Type type = sites.types.getOrDefault(offset, NONE);
            final SimpleInterval loc = new SimpleInterval(CONTIG, START + offset, START + offset);
            profile.add(new ActivityProfileState(loc, type == NONE ? 0.0 : 1.0, type, null));
        }

        return profile.popReadyAssemblyRegions(0, args.minAssemblyRegionSize, args.maxAssemblyRegionSize, forceConversion).stream()
                .map(region -> (region.getStart() - START) + "-" + (region.getEnd() - START) + ":" + region.isActive())
                .collect(Collectors.toList());
    }

    @DataProvider(name = "regions")
    public Object[][] regions() {
        return new Object[][] {
                // the example from the Mutect2ActivityProfile javadoc, with phasing distance 100: somatic variants at 200,
                // 350, and 600 and germline variants at 280, 650, and 800.  Germline variants within the phasing distance of
                // somatic variants join their active regions, but the lone germline variant at 800 is not active.
                {"javadoc example", 1000, new Sites().add(SOMATIC, 200, 350, 600).add(GERMLINE, 280, 650, 800), 0.0, 100, true,
                        List.of("0-184:false", "185-365:true", "366-584:false", "585-665:true", "666-784:false", "785-815:false",
                                "816-984:false", "985-999:false")},

                // issue 4: at the end of an interval, a leftover stretch of PADDING or fewer states with no active state is
                // inactive, even when the germline fraction is 1
                {"reference-only interval end", 100, new Sites(), 1.0, DEFAULT_PHASING_DISTANCE, true,
                        List.of("0-84:false", "85-99:false")},
                {"reference-only remainder after active region", 20, new Sites().add(GERMLINE, 0), 1.0, DEFAULT_PHASING_DISTANCE, true,
                        List.of("0-15:true", "16-19:false")},

                // issue 1: a germline variant within the phasing distance of a somatic variant at the end of an interval is
                // in the same region
                {"cluster at interval end", 70, new Sites().add(GERMLINE, 0).add(SOMATIC, 40), 0.0, DEFAULT_PHASING_DISTANCE, true,
                        List.of("0-55:true", "56-69:false")},
                // issue 1: a cluster that fits within the max region size is not split just because its end is within the
                // phasing distance of the max region size
                {"cluster ending near max size", 400, new Sites().add(GERMLINE, 0).add(SOMATIC, 40, 70, 100, 130, 160, 190, 220, 250, 280),
                        0.0, DEFAULT_PHASING_DISTANCE, false, List.of("0-295:true")},

                // issue 2: when a cluster is too large and must be split, the region before the split is not active due to a
                // somatic variant after the split.  Here the cluster is split after offset 0, at its longest gap.
                {"split cluster activity", 500, new Sites().add(GERMLINE, 0).add(GERMLINE, IntStream.rangeClosed(5, 30).map(n -> 10 * n).toArray())
                        .add(SOMATIC, 100).add(GERMLINE, 340), 0.0, DEFAULT_PHASING_DISTANCE, false,
                        List.of("0-15:false", "16-34:false", "35-349:true")},

                // issue 3: splitting a cluster with no gaps or with equal gaps yields large regions, not 16-base fragments
                {"contiguous active states", 400, new Sites().add(SOMATIC, IntStream.range(0, 400).toArray()), 0.0, DEFAULT_PHASING_DISTANCE, false,
                        List.of("0-314:true")},
                {"evenly spaced active states", 600, new Sites().add(SOMATIC, IntStream.rangeClosed(0, 10).map(n -> 49 * n).toArray()),
                        0.0, DEFAULT_PHASING_DISTANCE, false, List.of("0-309:true")},

                // the germline fraction applies to clusters with germline states, but panel of normals sites alone are never active
                {"panel of normals cluster", 100, new Sites().add(ARTIFACT, 20, 40), 1.0, DEFAULT_PHASING_DISTANCE, true,
                        List.of("0-4:false", "5-55:false", "56-84:false", "85-99:false")},
                {"germline and panel of normals cluster", 100, new Sites().add(GERMLINE, 20).add(ARTIFACT, 40), 1.0, DEFAULT_PHASING_DISTANCE, true,
                        List.of("0-4:false", "5-55:true", "56-84:false", "85-99:false")},

                // a region's trailing padding never reaches the next cluster, even when the phasing distance is shorter than
                // the padding.  The somatic variant at 32 is within the phasing distance of the germline variant at 40 but
                // not the one at 20.
                {"short phasing distance", 100, new Sites().add(GERMLINE, 20, 40).add(SOMATIC, 32), 0.0, 10, true,
                        List.of("0-4:false", "5-30:false", "31-50:true", "51-84:false", "85-99:false")},
                {"default phasing distance", 100, new Sites().add(GERMLINE, 20, 40).add(SOMATIC, 32), 0.0, DEFAULT_PHASING_DISTANCE, true,
                        List.of("0-4:false", "5-55:true", "56-84:false", "85-99:false")}
        };
    }

    @Test(dataProvider = "regions")
    public void testRegions(final String description, final int numStates, final Sites sites, final double germlineFraction,
                            final int maxPhasingDistance, final boolean forceConversion, final List<String> expected) {
        Assert.assertEquals(popRegions(numStates, sites, germlineFraction, maxPhasingDistance, forceConversion), expected, description);
    }

    // somatic variants at 200 and 600 and a germline variant at 230, with the alt alleles Mutect2Engine::isActive reports
    private static Map<Integer, ActivityProfileState> activeStates() {
        return Map.of(200, activeState(200, SOMATIC, ActivityProfileState.AltAlleleType.G),
                230, activeState(230, GERMLINE, ActivityProfileState.AltAlleleType.T),
                600, activeState(600, SOMATIC, ActivityProfileState.AltAlleleType.INSERTION));
    }

    private static ActivityProfileState activeState(final int offset, final ActivityProfileState.Type type, final ActivityProfileState.AltAlleleType altAlleleType) {
        return new ActivityProfileState(new SimpleInterval(CONTIG, START + offset, START + offset), 1.0, type, null, Map.of(altAlleleType, 10.0));
    }

    private static List<AssemblyRegion> popRegionsWithActiveStates(final ActivityProfile profile, final AssemblyRegionArgumentCollection args,
                                                                    final Map<Integer, ActivityProfileState> activeStates) {
        for (int offset = 0; offset < 1000; offset++) {
            final SimpleInterval loc = new SimpleInterval(CONTIG, START + offset, START + offset);
            profile.add(activeStates.getOrDefault(offset, new ActivityProfileState(loc, 0.0, NONE, null)));
        }
        return profile.popReadyAssemblyRegions(0, args.minAssemblyRegionSize, args.maxAssemblyRegionSize, true);
    }

    // Each region records the states of the active loci that formed it, along with their alt alleles, which lets Mutect2
    // tell when a region is active because of a lone SNV
    @Test
    public void testRegionsRecordActiveStates() {
        final AssemblyRegionArgumentCollection args = new AssemblyRegionArgumentCollection();
        args.genotypeGermlineSitesFraction = 0.0;
        final Map<Integer, ActivityProfileState> activeStates = activeStates();
        final List<AssemblyRegion> regions = popRegionsWithActiveStates(new Mutect2ActivityProfile(args, HEADER), args, activeStates);

        final List<AssemblyRegion> activeRegions = regions.stream().filter(AssemblyRegion::isActive).collect(Collectors.toList());
        Assert.assertEquals(activeRegions.size(), 2);
        assertSameStates(activeRegions.get(0).getActiveStates(), List.of(activeStates.get(200), activeStates.get(230)));
        assertSameStates(activeRegions.get(1).getActiveStates(), List.of(activeStates.get(600)));
        Assert.assertEquals(activeRegions.get(1).getActiveStates().get(0).getAltAlleleTypes(), Set.of(ActivityProfileState.AltAlleleType.INSERTION));

        regions.stream().filter(region -> !region.isActive()).forEach(region -> Assert.assertTrue(region.getActiveStates().isEmpty()));
    }

    // the band-pass profile's states are smoothed copies of those added, so its regions record none
    @Test
    public void testBandPassRegionsRecordNoActiveStates() {
        final AssemblyRegionArgumentCollection args = new AssemblyRegionArgumentCollection();
        final BandPassActivityProfile profile = new BandPassActivityProfile(args.maxProbPropagationDistance, args.activeProbThreshold,
                BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, HEADER);
        final List<AssemblyRegion> regions = popRegionsWithActiveStates(profile, args, activeStates());

        Assert.assertTrue(regions.stream().anyMatch(AssemblyRegion::isActive));
        regions.forEach(region -> Assert.assertTrue(region.getActiveStates().isEmpty()));
    }

    private static void assertSameStates(final List<ActivityProfileState> actual, final List<ActivityProfileState> expected) {
        Assert.assertEquals(actual.size(), expected.size());
        IntStream.range(0, actual.size()).forEach(n -> Assert.assertSame(actual.get(n), expected.get(n)));
    }
}
