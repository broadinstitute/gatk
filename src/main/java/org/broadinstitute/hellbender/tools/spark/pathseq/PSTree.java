package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a taxonomic tree with nodes assigned a name and taxonomic rank (e.g. order, family, genus, species, etc.)
 * Nodes are stored in a HashMap and keyed by a String id. Tree keeps track of each node's children as well as its
 * parent for efficient traversals top-down and bottom-up.
 * <p>
 * Designed to be constructed on-the-fly while parsing the NCBI taxonomy dump, which specifies each node and its parent.
 * As a result, invalid trees may be constructed (e.g. multiple roots or cycles). Use checkStructure() to confirm its validity.
 * <p>
 * Note the tree root is initialized in the constructor and cannot be modified (except adding children with addNode()).
 */
@DefaultSerializer(PSTree.Serializer.class)
public class PSTree {

    private final String root;
    private Map<String, PSTreeNode> tree;

    public PSTree(final String rootId) {
        tree = new HashMap<>();
        root = rootId;
        tree.put(root, new PSTreeNode());
        tree.get(root).setName("root");
        tree.get(root).setRank("root");
        tree.get(root).setParent(null);
    }

    @SuppressWarnings("unchecked")
    protected PSTree(final Kryo kryo, final Input input) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        root = input.readString();
        final int treeSize = input.readInt();
        tree = new HashMap<>(treeSize);
        for (int i = 0; i < treeSize; i++) {
            final String key = input.readString();
            final PSTreeNode value = kryo.readObject(input, PSTreeNode.class);
            tree.put(key, value);
        }

        kryo.setReferences(oldReferences);
    }

    /**
     * Returns short String of 20 arbitrarily chosen nodes
     */
    private static String getAbbreviatedNodeListString(final Set<String> nodes) {
        return "[" + nodes.stream().limit(20).collect(Collectors.joining(", ")) + "]";
    }

    protected void serialize(final Kryo kryo, final Output output) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        output.writeString(root);
        output.writeInt(tree.size());
        for (final String key : tree.keySet()) {
            output.writeString(key);
            kryo.writeObject(output, tree.get(key));
        }

        kryo.setReferences(oldReferences);
    }

    /**
     * Adds a node to the tree. ID must be unique, and all arguments must be non-null.
     * If the node exists, its properties are modified
     */
    public void addNode(final String id, final String name, final String parent, final long length, final String rank) {
        Utils.validateArg((id != null) && (name != null) && (parent != null) && (rank != null), "Passed a null argument to addNode()");
        Utils.validateArg(!root.equals(id), "Tried to set root attributes using addNode()");

        //If the node already exists, keep its current children and set everything else
        if (!tree.containsKey(id)) {
            tree.put(id, new PSTreeNode());
        }
        final PSTreeNode node = tree.get(id);
        node.setName(name);
        node.setParent(parent);
        node.setLength(length);
        node.setRank(rank);

        //If the parent doesn't exist, create a new node for it
        //We are assuming its attributes will be set later using addNode()
        if (!tree.containsKey(parent)) {
            tree.put(parent, new PSTreeNode());
        }
        tree.get(parent).addChild(id);
    }

    /**
     * Deletes nodes unreachable from the root and returns the set of deleted nodes.
     */
    public Set<String> removeUnreachableNodes() {
        final Set<String> reachable = traverse();
        final Set<String> unreachable = new HashSet<>(tree.keySet());
        unreachable.removeAll(reachable);
        retainNodes(reachable);
        return unreachable;
    }

    /**
     * Because of the piecemeal way we allow the tree to be constructed, we can end up with invalid tree structures.
     * Check that tree structure contains valid pointers and is connected.
     */
    public void checkStructure() {

        //Check child-parent pointers are consistent
        for (final String id : tree.keySet()) {
            final PSTreeNode n = tree.get(id);
            for (final String child : n.getChildren()) {
                if (!tree.containsKey(child)) {
                    throw new UserException.BadInput("Malformed tree detected. Node " + id + " has non-existent child " + child);
                }
                if (!tree.get(child).getParent().equals(id)) {
                    throw new UserException.BadInput("Malformed tree detected. Node " + id + " has child " + child + ", whose parent is " + tree.get(child).getParent() + " instead of " + id);
                }
            }
            final String parent = n.getParent();
            if (parent != null) {
                if (!tree.containsKey(parent)) {
                    throw new UserException.BadInput("Malformed tree detected. Node " + id + " has non-existent parent " + parent);
                }
                if (!tree.get(parent).getChildren().contains(id)) {
                    throw new UserException.BadInput("Malformed tree detected. Node " + id + " has parent " + parent + ", which does not have child " + id);
                }
            }
        }

        //Check disconnected sets of nodes using a traversal from the root
        final Set<String> unreachable = removeUnreachableNodes();
        if (!unreachable.isEmpty()) {
            final String nodesList = getAbbreviatedNodeListString(unreachable);
            throw new UserException.BadInput("Malformed tree detected. Tree is disconnected. " + unreachable.size() + " of " + tree.size() + " nodes were unreachable: " + nodesList);
        }
    }

    /**
     * Find all nodes reachable from the root
     */
    private Set<String> traverse() {
        final Queue<String> queue = new LinkedList<>();
        final Set<String> visited = new HashSet<>(tree.size());
        queue.add(root);
        while (!queue.isEmpty()) {
            final String id = queue.poll();
            if (!visited.contains(id)) { //checked visited in case there are cycles
                if (tree.containsKey(id)) {
                    queue.addAll(tree.get(id).getChildren());
                } else {
                    throw new UserException.BadInput("Could not find node " + id + " while traversing the tree");
                }
            } else {
                throw new UserException.BadInput("Tree contains a cycle. Attempted to visit node " + id + " twice during a breadth-first traversal.");
            }
            visited.add(id);
        }
        return visited;
    }

    /**
     * Get lowest common ancester of the set of given nodes.
     * Takes the intersection of node-to-root paths of all the nodes and finding the lowest one in the tree.
     */
    public String getLCA(final Collection<String> nodes) {
        Utils.validateArg(nodes.size() > 0, "Queried lowest common ancestor of a null set");
        final List<List<String>> paths = new ArrayList<>(nodes.size());
        for (final String node : nodes) {
            paths.add(getPathOf(node));
        }
        final List<String> firstPath = paths.remove(0);
        final Set<String> commonNodes = new HashSet<>(firstPath);
        for (final List<String> path : paths) {
            commonNodes.retainAll(path);
        }
        //Return first common node. Note paths are returned in order from lowest to highest (root at the end)
        for (final String node : firstPath) {
            if (commonNodes.contains(node)) return node;
        }
        //This should never happen if the tree structure has been checked
        throw new GATKException.ShouldNeverReachHereException("Could not find common ancester of node set.");
    }

    @SuppressWarnings("unchecked")
    public Collection<String> getChildrenOf(final String id) {
        Utils.validateArg(tree.containsKey(id), "Could not get children of node id " + id + " because it does not exist");
        return tree.get(id).getChildren();
    }

    public Set<String> getNodeIDs() {
        return tree.keySet();
    }

    public String getNameOf(final String id) {
        Utils.validateArg(tree.containsKey(id), "Could not get name of node id " + id + " because it does not exist");
        return tree.get(id).getName();
    }

    public String getParentOf(final String id) {
        Utils.validateArg(tree.containsKey(id), "Could not get parent of node id " + id + " because it does not exist");
        return tree.get(id).getParent();
    }

    public long getLengthOf(final String id) {
        Utils.validateArg(tree.containsKey(id), "Could not get length of node id " + id + " because it does not exist");
        return tree.get(id).getLength();
    }

    public String getRankOf(final String id) {
        Utils.validateArg(tree.containsKey(id), "Could not get rank of node id " + id + " because it does not exist");
        return tree.get(id).getRank();
    }

    public boolean hasNode(final String id) {
        return tree.containsKey(id);
    }

    /**
     * Removes all nodes not in the given set
     */
    public void retainNodes(final Set<String> idsToKeep) {
        Utils.validateArg(idsToKeep.contains(root), "Cannot remove root");
        final HashMap<String, PSTreeNode> newTree = new HashMap<>(idsToKeep.size());
        for (final String id : tree.keySet()) {
            if (idsToKeep.contains(id)) {
                final PSTreeNode info = tree.get(id);
                final PSTreeNode newInfo = info.copy();
                for (final String child : info.getChildren()) {
                    if (!idsToKeep.contains(child)) {
                        newInfo.removeChild(child);
                    }
                }
                if (!idsToKeep.contains(info.getParent())) {
                    newInfo.setParent(null);
                }
                newTree.put(id, newInfo);
            }
        }
        tree = newTree;
    }

    public void setNameOf(final String id, final String name) {
        Utils.validateArg(tree.containsKey(id), "Could not set name of node id " + id + " because it does not exist");
        Utils.validateArg(name != null, "Cannot set node name to null");
        Utils.validateArg(!id.equals(root), "Cannot set name of root");
        tree.get(id).setName(name);
    }

    public void setRankOf(final String id, final String rank) {
        Utils.validateArg(tree.containsKey(id), "Could not set rank of node id " + id + " because it does not exist");
        Utils.validateArg(rank != null, "Cannot set node rank to null");
        Utils.validateArg(!id.equals(root), "Cannot set rank of root");
        tree.get(id).setRank(rank);
    }

    public void setLengthOf(final String id, final long length) {
        Utils.validateArg(tree.containsKey(id), "Could not set rank of node id " + id + " because it does not exist");
        Utils.validateArg(!id.equals(root), "Cannot set rank of root");
        tree.get(id).setLength(length);
    }

    /**
     * Returns path of node id's from the input id to the root.
     */
    public List<String> getPathOf(String id) {
        final List<String> path = new ArrayList<>();
        final Set<String> visitedNodes = new HashSet<>(tree.size());
        while (id != null) {
            if (!visitedNodes.contains(id)) {
                visitedNodes.add(id);
                if (tree.containsKey(id)) {
                    path.add(id);
                    id = tree.get(id).getParent();
                } else {
                    throw new UserException.BadInput("Parent node " + id + " not found in tree while getting path");
                }
            } else {
                throw new UserException.BadInput("The tree contains a cycle at node " + id);
            }
        }
        return path;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final PSTree psTree = (PSTree) o;

        return root.equals(psTree.root) && tree.equals(psTree.tree);
    }

    @Override
    public int hashCode() {
        int result = root.hashCode();
        result = 31 * result + tree.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return getAbbreviatedNodeListString(tree.keySet());
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PSTree> {
        @Override
        public void write(final Kryo kryo, final Output output, final PSTree psTree) {
            psTree.serialize(kryo, output);
        }

        @Override
        public PSTree read(final Kryo kryo, final Input input,
                           final Class<PSTree> klass) {
            return new PSTree(kryo, input);
        }
    }

}
