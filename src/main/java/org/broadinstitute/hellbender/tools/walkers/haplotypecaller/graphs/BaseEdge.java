package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Comparator;

/**
 * Simple edge class for connecting nodes in the graph.
 *
 * Works equally well for all graph types (kmer or sequence)
 */
public class BaseEdge {
    private int multiplicity;
    private boolean isRef;

    /**
     * Create a new BaseEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param multiplicity the number of observations of this edge
     */
    public BaseEdge(final boolean isRef, final int multiplicity) {
        if ( multiplicity < 0 ) {
            throw new IllegalArgumentException("multiplicity must be >= 0 but got " + multiplicity);
        }

        this.multiplicity = multiplicity;
        this.isRef = isRef;
    }

    /**
     * Create a new copy of this BaseEdge
     */
    public BaseEdge copy() {
        return new BaseEdge(isRef(), getMultiplicity());
    }

    /**
     * Get the number of observations of paths connecting two vertices
     * @return a positive integer >= 0
     */
    public final int getMultiplicity() {
        return multiplicity;
    }

    /**
     * Get the DOT format label for this edge, to be displayed when printing this edge to a DOT file
     * @return a non-null string
     */
    public String getDotLabel() {
        return Integer.toString(getMultiplicity());
    }

    /**
     * Increase the multiplicity of this edge by incr
     * @param incr the change in this multiplicity, must be >= 0
     */
    public void incMultiplicity(final int incr) {
        if ( incr < 0 ) {
            throw new IllegalArgumentException("incr must be >= 0 but got " + incr);
        }
        multiplicity += incr;
    }

    /**
     * A special assessor that returns the multiplicity that should be used by pruning algorithm
     *
     * @return the multiplicity value that should be used for pruning
     */
    public int getPruningMultiplicity() {
        return getMultiplicity();
    }

    /**
     * Set the multiplicity of this edge to value
     * @param value an integer >= 0
     */
    public final void setMultiplicity( final int value ) {
        if ( multiplicity < 0 ) {
            throw new IllegalArgumentException("multiplicity must be >= 0");
        }
        multiplicity = value;
    }

    /**
     * Does this edge indicate a path through the reference graph?
     * @return true if so
     */
    public final boolean isRef() {
        return isRef;
    }

    /**
     * Indicate that this edge follows the reference sequence, or not
     * @param isRef true if this is a reference edge
     */
    public final void setIsRef( final boolean isRef ) {
        this.isRef = isRef;
    }

    /**
     * Sorts a collection of BaseEdges in decreasing order of weight, so that the most
     * heavily weighted is at the start of the list
     */
    public static final Comparator<BaseEdge> EDGE_MULTIPLICITY_ORDER = Comparator.comparingInt(BaseEdge::getMultiplicity).reversed();

    /**
     * Add edge to this edge, updating isRef and multiplicity as appropriate
     *
     * isRef is simply the or of this and edge
     * multiplicity is the sum
     *
     * @param edge the edge to add
     * @return this
     */
    public final BaseEdge add(final BaseEdge edge) {
        Utils.nonNull(edge, "edge cannot be null");
        this.multiplicity += edge.getMultiplicity();
        this.isRef = this.isRef || edge.isRef();
        return this;
    }

    /**
     * Create a new BaseEdge with the given multiplicity. The resulting edge is a reference edge if any of the argument edges are reference.
     *
     * @param edges a collection of edges to or their isRef values
     * @param multiplicity our desired multiplicity
     * @return a newly allocated BaseEdge
     */
    public static BaseEdge makeOREdge(final Collection<BaseEdge> edges, final int multiplicity) {
        Utils.nonNull(edges);
        final boolean anyRef = edges.stream().anyMatch(e -> e.isRef());
        return new BaseEdge(anyRef, multiplicity);
    }

    @Override
    public final String toString() {
        return String.format("BaseEdge{multiplicity=%d, isRef=%b}", multiplicity, isRef);
    }
}
