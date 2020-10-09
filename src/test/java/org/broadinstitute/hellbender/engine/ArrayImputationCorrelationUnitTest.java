package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.validation.ArrayImputationCorrelation;

import java.util.Iterator;
import java.util.List;

public class ArrayImputationCorrelationUnitTest extends GATKBaseTest {


    public void testPearsonCorrelationAggregator(final List<Double> xList, final List<Double> yList) {
        final ArrayImputationCorrelation.PearsonCorrelationAggregator correlationAggregator = new ArrayImputationCorrelation.PearsonCorrelationAggregator();
        final Iterator<Double> xIterator = xList.iterator();
        final Iterator<Double> yIterator = yList.iterator();
        while (xIterator.hasNext() && yIterator.hasNext()) {
            correlationAggregator.addEntry(xIterator.next(), yIterator.next());
        }
    }
}