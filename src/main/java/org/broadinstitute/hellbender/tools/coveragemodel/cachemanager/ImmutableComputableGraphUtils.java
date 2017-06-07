package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A utilities class for {@link ImmutableComputableGraph}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ImmutableComputableGraphUtils {

    private ImmutableComputableGraphUtils() {}

    /**
     * A simple builder class for {@link ImmutableComputableGraph}.
     *
     * @implNote Node addition methods must perform node key uniqueness checks. Otherwise, some of the nodes
     * with the same key will be lost; see {@link CacheNode#equals(Object)}.
     */
    public static class ImmutableComputableGraphBuilder {
        private final Set<CacheNode> nodes;
        private final Set<CacheNode.NodeKey> keys;
        private boolean cacheAutoUpdate;

        ImmutableComputableGraphBuilder() {
            nodes = new HashSet<>();
            keys = new HashSet<>();
            cacheAutoUpdate = false;
        }

        /**
         * Add a primitive node with specified value
         */
        public ImmutableComputableGraphBuilder primitiveNode(@Nonnull final CacheNode.NodeKey key,
                                                             @Nonnull final CacheNode.NodeTag[] tags,
                                                             @Nonnull Duplicable value) {
            Utils.nonNull(key);
            Utils.nonNull(tags);
            Utils.nonNull(value);
            assertKeyUniqueness(key);
            nodes.add(new PrimitiveCacheNode(key, Arrays.stream(tags).collect(Collectors.toList()), value));
            keys.add(key);
            return this;
        }

        /**
         * Add an uninitialized primitive node holding an empty {@link DuplicableNDArray} and no tags
         */
        public ImmutableComputableGraphBuilder primitiveNodeWithEmptyNDArray(@Nonnull final CacheNode.NodeKey key) {
            return primitiveNode(key, new CacheNode.NodeTag[] {}, new DuplicableNDArray());
        }

        /**
         * Add a computable node
         */
        public ImmutableComputableGraphBuilder computableNode(@Nonnull final CacheNode.NodeKey key,
                                                              @Nonnull final CacheNode.NodeTag[] tags,
                                                              @Nonnull final CacheNode.NodeKey[] parents,
                                                              @Nullable final ComputableNodeFunction func,
                                                              final boolean cacheEvals) {
            Utils.nonNull(key);
            Utils.nonNull(tags);
            Utils.nonNull(parents);
            assertKeyUniqueness(key);
            nodes.add(new ComputableCacheNode(key,
                    Arrays.stream(tags).collect(Collectors.toList()),
                    Arrays.stream(parents).collect(Collectors.toList()),
                    func, cacheEvals));
            keys.add(key);
            return this;
        }

        /**
         * Add an "untracked" externally computable node (i.e. no parents and no tags). Note: it is the user's
         * responsibility to keep track of the value of these nodes and updating them if necessary. Since these nodes
         * have no parents, {@link ImmutableComputableGraph} will <b>not</b> perform any bookkeeping on them.
         */
        public ImmutableComputableGraphBuilder untrackedExternallyComputableNode(@Nonnull final CacheNode.NodeKey key) {
            return computableNode(key, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {}, null, true);
        }

        /**
         * Enable cache auto-update for all computable nodes in the graph
         */
        public ImmutableComputableGraphBuilder withCacheAutoUpdate() {
            cacheAutoUpdate = true;
            return this;
        }

        /**
         * Disable cache auto-update for all computable nodes in the graph
         */
        public ImmutableComputableGraphBuilder withoutCacheAutoUpdate() {
            cacheAutoUpdate = false;
            return this;
        }

        private void assertKeyUniqueness(@Nonnull final CacheNode.NodeKey key) {
            if (keys.contains(key)) {
                throw new DuplicateNodeKeyException("A node with key " + quote(key.toString()) + " already exists");
            }
        }

        public ImmutableComputableGraph build() {
            if (nodes.size() == 0) {
                throw new IllegalStateException("Can not make an empty cache node collection");
            } else {
                return new ImmutableComputableGraph(nodes, cacheAutoUpdate);
            }
        }

        /**
         * This exception will be thrown if a node with the same key is already added to the builder
         */
        static final class DuplicateNodeKeyException extends RuntimeException {
            private static final long serialVersionUID = 2016242121833170379L;

            DuplicateNodeKeyException(String s) {
                super(s);
            }
        }
    }

    static String quote(final String str) {
        return "\"" + str + "\"";
    }
}