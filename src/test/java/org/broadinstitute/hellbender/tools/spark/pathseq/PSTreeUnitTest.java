package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class PSTreeUnitTest extends CommandLineProgramTest {

    @Test
    public void testSerialize() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");

        try {
            final File tempFile = createTempFile("test", ".dat");
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Output output = new Output(new FileOutputStream(tempFile));
            kryo.writeObject(output, tree);
            output.close();

            final Input input = new Input(new FileInputStream(tempFile));
            final PSTree treeIn = kryo.readObject(input, PSTree.class);
            Assert.assertEquals(treeIn, tree);
        } catch (FileNotFoundException e) {
            throw new IOException("Error with Kryo IO", e);
        }
    }

    @Test
    public void testAddNode() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        try {
            tree.addNode(1, "new_root", 1, 0, "root");
            Assert.fail("Was able to add root node");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.addNode(PSTree.NULL_NODE, "n3", 3, 0, "none");
            Assert.fail("Was able to add node with null id");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.addNode(3, null, 3, 0, "none");
            Assert.fail("Was able to add node with null name");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.addNode(3, "n3", PSTree.NULL_NODE, 0, "none");
            Assert.fail("Was able to add node with null parent");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.addNode(3, "n3", 3, 0, null);
            Assert.fail("Was able to add node with null rank");
        } catch (IllegalArgumentException e) {
        }
    }

    @Test
    public void testCheckStructure() throws Exception {
        PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(4, "n4", 2, 0, "none");
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.checkStructure();

        tree = new PSTree(1);
        tree.addNode(2, "n2", 2, 0, "none");
        try {
            tree.checkStructure();
            Assert.fail("Tree validated when node self-referencing");
        } catch (UserException.BadInput e) {
        }

        tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");
        tree.addNode(4, "n4", 3, 0, "none");
        try {
            tree.checkStructure();
            Assert.fail("Tree validated when node parent-child points were incorrect");
        } catch (UserException.BadInput e) {
        }
    }

    @Test
    public void testGetChildrenOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");
        final Integer[] expectedChildren = {2, 3};
        Assert.assertEquals(tree.getChildrenOf(1).toArray(), expectedChildren);
    }

    @Test
    public void testGetNodeIDs() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        final Set<Integer> expectedIDs = new HashSet<>();
        expectedIDs.add(1);
        expectedIDs.add(2);
        Assert.assertEquals(tree.getNodeIDs(), expectedIDs);
    }

    @Test
    public void testGetNameOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "rank_of_2");
        Assert.assertEquals(tree.getNameOf(1), "root", "Name of root was not 'root'");
        Assert.assertEquals(tree.getNameOf(2), "n2");
    }

    @Test
    public void testGetParentOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        Assert.assertEquals(tree.getParentOf(1), PSTree.NULL_NODE, "Parent of root was not null");
        Assert.assertEquals(tree.getParentOf(2), 1);
    }

    @Test
    public void testGetRankOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "rank_of_2");
        Assert.assertEquals(tree.getRankOf(1), "root", "Rank of root was not 'root'");
        Assert.assertEquals(tree.getRankOf(2), "rank_of_2");
    }

    @Test
    public void testGetLengthOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 10, "rank_of_2");
        Assert.assertEquals(tree.getLengthOf(1), 0, "Length of root was not 0");
        Assert.assertEquals(tree.getLengthOf(2), 10);
    }

    @Test
    public void testHasNode() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        Assert.assertTrue(tree.hasNode(1), "Tree did not contain root");
        Assert.assertTrue(tree.hasNode(2));
        Assert.assertFalse(tree.hasNode(3));
    }

    @Test
    public void testRetainNodes() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");

        final Set<Integer> nodesToKeep = new HashSet<>();
        nodesToKeep.add(1);
        nodesToKeep.add(2);
        nodesToKeep.add(4);
        tree.retainNodes(nodesToKeep);
        Assert.assertEquals(tree.getNodeIDs(), nodesToKeep);
    }

    @Test
    public void testSetNameOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.setNameOf(2, "node_2");
        Assert.assertEquals(tree.getNameOf(2), "node_2");
        try {
            tree.setNameOf(1, "not_root");
            Assert.fail("Was able to set root name");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.setNameOf(PSTree.NULL_NODE, "");
            Assert.fail("Was able to set node name with id null");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.setNameOf(2, null);
            Assert.fail("Was able to set node name to null");
        } catch (IllegalArgumentException e) {
        }
    }

    @Test
    public void testSetRankOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.setRankOf(2, "rank_2");
        Assert.assertEquals(tree.getRankOf(2), "rank_2");
        try {
            tree.setRankOf(1, "not_root");
            Assert.fail("Was able to set root rank");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.setRankOf(PSTree.NULL_NODE, "");
            Assert.fail("Was able to set node rank with id null");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.setRankOf(2, null);
            Assert.fail("Was able to set node rank to null");
        } catch (IllegalArgumentException e) {
        }
    }

    @Test
    public void testSetLengthOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.setLengthOf(2, 5);
        Assert.assertEquals(tree.getLengthOf(2), 5);
        try {
            tree.setLengthOf(1, 5);
            Assert.fail("Was able to set root length");
        } catch (IllegalArgumentException e) {
        }
        try {
            tree.setLengthOf(PSTree.NULL_NODE, 5);
            Assert.fail("Was able to set node length with id null");
        } catch (IllegalArgumentException e) {
        }
    }

    @Test
    public void testGetPathOf() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");
        final Integer[] expectedPath = {4, 2, 1};
        final List<Integer> path = tree.getPathOf(4);
        Assert.assertEquals(path.toArray(), expectedPath);
        final List<Integer> path_null = tree.getPathOf(PSTree.NULL_NODE);
        Assert.assertEquals(path_null, new ArrayList<String>());

        try {
            tree.getPathOf(5);
            Assert.fail("Did not throw UserException when asking for path of non-existent node");
        } catch (UserException e) {
        }
    }

    @Test
    public void testGetLCA() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");
        tree.addNode(5, "n5", 3, 0, "none");
        tree.addNode(6, "n6", 3, 0, "none");
        tree.addNode(7, "n7", 5, 0, "none");

        final ArrayList<Integer> nodes = new ArrayList<>();
        try {
            tree.getLCA(nodes);
            Assert.fail("Did not throw exception when asking for LCA of empty set");
        } catch (Exception e) {
        }

        int lca;

        nodes.clear();
        nodes.add(1);
        lca = tree.getLCA(nodes);
        Assert.assertEquals(1, lca);

        nodes.clear();
        nodes.add(4);
        nodes.add(5);
        lca = tree.getLCA(nodes);
        Assert.assertEquals(1, lca);

        nodes.clear();
        nodes.add(6);
        nodes.add(7);
        lca = tree.getLCA(nodes);
        Assert.assertEquals(3, lca);

        nodes.clear();
        nodes.add(5);
        nodes.add(7);
        lca = tree.getLCA(nodes);
        Assert.assertEquals(5, lca);

        tree.addNode(8, "n8", 9, 0, "none");
        nodes.clear();
        nodes.add(8);
        nodes.add(2);
        try {
            tree.getLCA(nodes);
            Assert.fail("Did not throw exception when asking for LCA a disconnected node");
        } catch (Exception e) {
        }
    }

    @Test
    public void testRemoveUnreachableNodes() throws Exception {
        PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        Assert.assertEquals(tree.removeUnreachableNodes(), new HashSet<String>());

        tree.addNode(4, "n4", 8, 0, "none");
        tree.addNode(5, "n5", 6, 0, "none");
        final HashSet<Integer> trueUnreachable = new HashSet<>(2);
        trueUnreachable.add(4);
        trueUnreachable.add(5);
        trueUnreachable.add(6);
        trueUnreachable.add(8);
        Assert.assertEquals(tree.removeUnreachableNodes(), trueUnreachable);
    }

    @Test
    public void testToString() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.checkStructure();
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        tree.addNode(4, "n4", 2, 0, "none");
        Assert.assertTrue(!tree.toString().isEmpty());
        tree.toString();
    }

    @Test
    public void testHashCode() throws Exception {
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "n2", 1, 0, "none");
        tree.addNode(3, "n3", 1, 0, "none");
        final PSTree tree2 = new PSTree(1);
        tree2.addNode(2, "n2", 1, 0, "none");
        tree2.addNode(3, "n3", 1, 0, "none");
        final PSTree tree3 = new PSTree(1);
        tree3.addNode(2, "n2", 1, 0, "none");
        tree3.addNode(3, "not_n3", 1, 0, "none");
        Assert.assertEquals(tree.hashCode(), tree2.hashCode());
        Assert.assertNotEquals(tree.hashCode(), tree3.hashCode());
    }

}