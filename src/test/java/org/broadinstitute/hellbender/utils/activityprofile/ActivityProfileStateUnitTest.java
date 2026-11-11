package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public final class ActivityProfileStateUnitTest {
    private SAMSequenceDictionary sequenceDictionary;

    @BeforeClass
    public void init() {
        // sequence
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100);
        sequenceDictionary =header.getSequenceDictionary();
    }

    @DataProvider(name = "ActiveProfileResultProvider")
    public Object[][] makeActiveProfileResultProvider() {
        final List<Object[]> tests = new LinkedList<>();

        final String chr = sequenceDictionary.getSequence(0).getSequenceName();
        for ( final SimpleInterval loc : Arrays.asList(
                new SimpleInterval(chr, 10, 10),
                new SimpleInterval(chr, 100, 100))) {
            for ( final double prob : Arrays.asList(0.0, 0.5, 1.0) ) {
                for ( final ActivityProfileState.Type state : ActivityProfileState.Type.values() ) {
                    for ( final Number value : Arrays.asList(1, 2, 4) ) {
                        tests.add(new Object[]{ loc, prob, state, value});
                    }
                }
                tests.add(new Object[]{ loc, prob, null, null});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ActiveProfileResultProvider")
    public void testActiveProfileResultProvider(SimpleInterval loc, final double prob, ActivityProfileState.Type maybeState, final Number maybeNumber) {
        final ActivityProfileState apr = maybeState == null
                ? new ActivityProfileState(loc, prob)
                : new ActivityProfileState(loc, prob, maybeState, maybeNumber);

        Assert.assertEquals(apr.getLoc(), loc);
        Assert.assertEquals(apr.getOffset(loc), 0);
        Assert.assertNotNull(apr.toString());
        Assert.assertEquals(apr.isActiveProb(), prob);
        Assert.assertEquals(apr.getResultState(), maybeState == null ? ActivityProfileState.Type.NONE : maybeState);
        Assert.assertEquals(apr.getResultValue(), maybeState == null ? null : maybeNumber);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testError1(){
        final String chr = sequenceDictionary.getSequence(0).getSequenceName();
        new ActivityProfileState(new SimpleInterval(chr, 10, 10), 0.1, null, -1.0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testError2(){
        final String chr = sequenceDictionary.getSequence(0).getSequenceName();
        new ActivityProfileState(new SimpleInterval(chr, 10, 11), 0.1, null, 1.0);
    }

    @Test
    public void testAltAlleleLogOdds() {
        final SimpleInterval loc = new SimpleInterval(sequenceDictionary.getSequence(0).getSequenceName(), 10, 10);

        // states that don't report alt alleles have none
        Assert.assertTrue(new ActivityProfileState(loc, 1.0).getAltAlleleLogOdds().isEmpty());
        Assert.assertTrue(new ActivityProfileState(loc, 1.0, ActivityProfileState.Type.SOMATIC, null).getAltAlleleTypes().isEmpty());

        // the state keeps its own copy of the reported alt alleles
        final Map<ActivityProfileState.AltAlleleType, Double> logOdds = new EnumMap<>(ActivityProfileState.AltAlleleType.class);
        logOdds.put(ActivityProfileState.AltAlleleType.G, 12.0);
        logOdds.put(ActivityProfileState.AltAlleleType.DELETION, 5.0);
        final ActivityProfileState state = new ActivityProfileState(loc, 1.0, ActivityProfileState.Type.SOMATIC, null, logOdds);
        logOdds.put(ActivityProfileState.AltAlleleType.T, 7.0);

        Assert.assertEquals(state.getAltAlleleLogOdds(), Map.of(ActivityProfileState.AltAlleleType.G, 12.0, ActivityProfileState.AltAlleleType.DELETION, 5.0));
        Assert.assertEquals(state.getAltAlleleTypes(), EnumSet.of(ActivityProfileState.AltAlleleType.G, ActivityProfileState.AltAlleleType.DELETION));
    }

    @Test
    public void testSubstitutedBase() {
        Assert.assertEquals(ActivityProfileState.AltAlleleType.A.substitutedBase(), (byte) 'A');
        Assert.assertEquals(ActivityProfileState.AltAlleleType.C.substitutedBase(), (byte) 'C');
        Assert.assertEquals(ActivityProfileState.AltAlleleType.G.substitutedBase(), (byte) 'G');
        Assert.assertEquals(ActivityProfileState.AltAlleleType.T.substitutedBase(), (byte) 'T');
        Assert.assertEquals(ActivityProfileState.AltAlleleType.INSERTION.substitutedBase(), (byte) 0);
        Assert.assertEquals(ActivityProfileState.AltAlleleType.DELETION.substitutedBase(), (byte) 0);
        Assert.assertEquals(ActivityProfileState.AltAlleleType.OTHER.substitutedBase(), (byte) 0);
    }
}
