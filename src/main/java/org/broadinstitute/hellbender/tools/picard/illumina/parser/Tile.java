package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CollectionUtil;

import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.getSoleElement;
import static htsjdk.samtools.util.CollectionUtil.partition;
import static java.util.Arrays.asList;
import static java.util.Collections.unmodifiableMap;

/**
 * Represents a tile from TileMetricsOut.bin. Stores information on location (lane & tile #, density, number of clusters and the
 * phasing/prephasing values associated with this tile
 *
 * @author jgentry
 */
public final class Tile {
    private final int lane, tile;
    private final float density, clusters;

    private final Map<TileTemplateRead, Float> phasingMap;
    private final Map<TileTemplateRead, Float> prePhasingMap;

    /**
     * @param tilePhasingValues Either one or two TilePhasingValues, corresponding to the FIRST and potentially SECOND template reads
     */
    public Tile(final int lane, final int tile, final float density, final float clusters, final TilePhasingValue... tilePhasingValues) {
        this.lane = lane;
        this.tile = tile;
        this.density = density;
        this.clusters = clusters;

        final Collection<TilePhasingValue> phasingValues = ensureSoleTilePhasingValuesPerRead(asList(tilePhasingValues));

        final Map<TileTemplateRead, Float> phasingMap = new HashMap<TileTemplateRead, Float>();
        final Map<TileTemplateRead, Float> prePhasingMap = new HashMap<TileTemplateRead, Float>();

        /** For each of the TileReads, assign their phasing & prephasing values to the respective maps, which we will
         * use later to calculate the medians
         */
        for (final TilePhasingValue phasingValue : phasingValues) {
            phasingMap.put(phasingValue.getTileTemplateRead(), phasingValue.getPhasingValue());
            prePhasingMap.put(phasingValue.getTileTemplateRead(), phasingValue.getPrePhasingValue());
        }

        this.phasingMap = unmodifiableMap(phasingMap);
        this.prePhasingMap = unmodifiableMap(prePhasingMap);
    }

    /**
     * Returns the number of this tile's parent lane.
     */
    public int getLaneNumber() {
        return lane;
    }

    /**
     * Returns the number/name of this tile.
     */
    public int getTileNumber() {
        return tile;
    }

    /**
     * Returns the cluster density of this tile, in units of [cluster/mm^2].
     */
    public float getClusterDensity() {
        return density;
    }

    /**
     * Returns the number of on this tile.
     */
    public float getClusterCount() {
        return clusters;
    }

    public Map<TileTemplateRead, Float> getPhasingMap() {
        return phasingMap;
    }

    public Map<TileTemplateRead, Float> getPrePhasingMap() {
        return prePhasingMap;
    }

    /**
     * For any given TileTemplateRead, we want to make sure that there is only a single TilePhasingValue
     */
    private static Collection<TilePhasingValue> ensureSoleTilePhasingValuesPerRead(final Collection<TilePhasingValue> tilePhasingValues) {
        final Map<TileTemplateRead, Collection<TilePhasingValue>> partitionedMap = partition(tilePhasingValues,
                new CollectionUtil.Partitioner<TilePhasingValue, TileTemplateRead>() {
                    @Override
                    public TileTemplateRead getPartition(final TilePhasingValue phasingValue) {
                        return phasingValue.getTileTemplateRead();
                    }
                });

        final Collection<TilePhasingValue> newTilePhasingValues = new LinkedList<TilePhasingValue>();
        for (final TileTemplateRead read : partitionedMap.keySet()) {
            newTilePhasingValues.add(getSoleElement(partitionedMap.get(read)));
        }

        return newTilePhasingValues;
    }
}
