package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class LongHomopolymerHaplotypeCollapsingEngineUnitTest extends GATKBaseTest {

    private static final Logger logger = LogManager.getLogger(LongHomopolymerHaplotypeCollapsingEngineUnitTest.class);

    @DataProvider(name = "needsCollapsingDataProvider")
    public Object[][] needsCollapsingDataProvider() {
        return new Object[][] {
                { "CCAATTGG", 12, false },
                { "CCAAAAAAAAAAAATTGG", 12, false},
                { "CCAAAAAAAAAAAAATTGG", 12, true},
        };
    }

    @Test(dataProvider = "needsCollapsingDataProvider")
    public void testNeedsCollapsing(final String bases, final int threshold, final boolean expected) {

        final boolean result = LongHomopolymerHaplotypeCollapsingEngine.needsCollapsing(bases.getBytes(), threshold, logger);

        Assert.assertEquals(result, expected);
    }

    @DataProvider(name = "identicalBySequenceDataProvider")
    public Object[] identicalBySequenceDataProvider() {
        return new Object[][] {
                {
                    new String[] {"AAAA", "AAAA/r", "CCCC", "AAAA"},        // /r marks a reference
                        new String[] {"1-0,1,3", "2-2"}                   // result map

                }
        };
    }

    @Test(dataProvider = "identicalBySequenceDataProvider")
    public void testIdenticalBySequence(final String[] haplotypesBases, final String[] expectedMapSpecs) {

        // this code is a bit overly complex by the need to actually verify collection inclusion without
        // using equal() - as soom of the input are equal()

        // create haplotypes
        List<Haplotype>         haplotypes = new LinkedList<>();;
        for ( final String haplotypeBases : haplotypesBases ) {
            final Haplotype     haplotype = createHaplotype(haplotypeBases);
            haplotypes.add(haplotype);
        }

        // execute the call
        final Map<Haplotype, List<Haplotype>> result = LongHomopolymerHaplotypeCollapsingEngine.identicalBySequence(haplotypes);

        // create expected (index) map
        final Map<Integer, Set<Integer>>      expected = createIndexMap(expectedMapSpecs);

        // compare maps
        Assert.assertEquals(result.size(), expected.size());
        expected.forEach((keyIndex, valuesIndexs) -> {

            // access expected key haplotype
            final Haplotype       keyHaplotype = haplotypes.get(keyIndex);

            // the key must be present in the result map
            final List<List<Haplotype>>   valueHaplotypesList = result.entrySet().stream()
                    .filter(e ->  e.getKey() == keyHaplotype)
                    .map(e -> e.getValue())
                    .collect(Collectors.toList());
            Assert.assertEquals(valueHaplotypesList.size(), 1);
            final List<Haplotype>  valueHaplotypes = valueHaplotypesList.get(0);

            // check values
            Assert.assertEquals(valueHaplotypes.size(), valuesIndexs.size());
            valuesIndexs.forEach(valueIndex -> {

                // value must be present in the set
                final Haplotype     valueHaplotype = haplotypes.get(valueIndex);

                Assert.assertEquals(valueHaplotypes.stream()
                        .filter(haplotype -> haplotype == valueHaplotype)
                        .collect(Collectors.toList()).size(), 1);
            });

        });
    }

    private Haplotype createHaplotype(final String bases) {
        if ( bases.endsWith("/r") ) {
            return new Haplotype(bases.substring(0, bases.length() - 2).getBytes(), true);
        } else {
            return new Haplotype(bases.getBytes(), false);
        }
    }

    private Map<Integer, Set<Integer>> createIndexMap(final String[] specs)
    {
        Map<Integer, Set<Integer>>     map = new LinkedHashMap<>();

        for ( final String spec : specs ) {
            final String[]      toks = spec.split("-");
            final Set<Integer>  values = new LinkedHashSet<>();
            for ( final String value : toks[1].split(",") ) {
                values.add(Integer.parseInt(value));
            }
            map.put(Integer.parseInt(toks[0]), values);
        }

        return map;
    }
}
