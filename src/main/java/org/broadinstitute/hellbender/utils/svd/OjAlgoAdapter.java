package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.ojalgo.access.Access2D;
import org.ojalgo.matrix.store.ElementsConsumer;
import org.ojalgo.matrix.store.ElementsSupplier;
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.matrix.store.RawStore;

public class OjAlgoAdapter implements Access2D<Double>, ElementsSupplier<Double> {

    public static OjAlgoAdapter of(final RealMatrix adaptee) {
        return new OjAlgoAdapter(adaptee);
    }

    private final RealMatrix myAdaptee;

    OjAlgoAdapter(final RealMatrix adaptee) {
        super();
        myAdaptee = adaptee;
    }

    public long countColumns() {
        return myAdaptee.getColumnDimension();
    }

    public long countRows() {
        return myAdaptee.getRowDimension();
    }

    public double doubleValue(final long row, final long column) {
        return myAdaptee.getEntry((int) row, (int) column);
    }

    public PhysicalStore.Factory<Double, RawStore> factory() {
        // return PrimitiveDenseStore.FACTORY;
        return RawStore.FACTORY;
    }

    public Double get(final long row, final long column) {
        return this.doubleValue(row, column);
    }

    public void supplyTo(final ElementsConsumer<Double> consumer) {
        final int tmpRowDim = myAdaptee.getRowDimension();
        final int tmpColDim = myAdaptee.getColumnDimension();
        for (int j = 0; j < tmpColDim; j++) {
            for (int i = 0; i < tmpRowDim; i++) {
                consumer.set(i, j, myAdaptee.getEntry(i, j));
            }
        }
    }

}
