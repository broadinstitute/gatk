package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.util.Collections;
import java.util.List;

/**
 * Empty meaning no reads thus elements.
 */
class EmptyIntervalPileup implements IntervalPileup {

    private final ReferenceBases referenceBases;
    EmptyIntervalPileup(final ReferenceBases referenceBases) {
        this.referenceBases = referenceBases;
    }

    @Override
    @SuppressWarnings("unchecked")
    public List<GATKRead> reads() {
        return Collections.EMPTY_LIST;
    }

    @Override
    public ReferenceBases reference() {
        return referenceBases;
    }

    @Override
    public int width() {
        return referenceBases.getInterval().size();
    }

    @Override
    public int height() {
        return 0;
    }

    @Override
    public byte baseAt(final int row, final int column) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public byte qualAt(final int row, final int column) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public Insert insertAt(int row, int column) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public GATKRead readAt(final int row, final int column) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public List<GATKRead> readsAt(final int row, final int column) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public Element element(final GATKRead read) {
        return null;
    }

    @Override
    public boolean hasInsertAt(int i, int j) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public byte[] insertBasesAt(int i, int j) {
        throw new IndexOutOfBoundsException();
    }

    @Override
    public byte[] insertQualsAt(int i, int j) {
        throw new IndexOutOfBoundsException();
    }
}
