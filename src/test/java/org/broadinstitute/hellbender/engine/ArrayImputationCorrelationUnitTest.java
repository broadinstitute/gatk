package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.GATKBaseTest;

import java.util.Iterator;
import java.util.List;

import static org.testng.Assert.*;

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