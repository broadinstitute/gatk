package org.broadinstitute.hellbender.utils.diffengine;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public final class DiffNodeUnitTest extends BaseTest {
    // Data is:
    // MY_ROOT
    //   fields: A=A, B=B
    //   nodes: C, D
    //   C: fields: E=E, nodes: none
    //   D: fields: F=F, G=G, nodes: none
    static DiffNode MY_ROOT = DiffNode.rooted("MY_ROOT");
    static DiffValue Value_A = new DiffValue("A", MY_ROOT, "A");
    static DiffValue Value_B = new DiffValue("B", MY_ROOT, "B");
    static DiffNode NODE_C = DiffNode.empty("C", MY_ROOT);
    static DiffNode NODE_D = DiffNode.empty("D", MY_ROOT);
    static DiffValue Value_E = new DiffValue("E", NODE_C, "E");
    static DiffValue Value_F = new DiffValue("F", NODE_D, "F");
    static DiffValue Value_G = new DiffValue("G", NODE_D, "G");

    static {
        MY_ROOT.add(Value_A);
        MY_ROOT.add(Value_B);
        MY_ROOT.add(NODE_C);
        MY_ROOT.add(NODE_D);
        NODE_C.add(Value_E);
        NODE_D.add(Value_F);
        NODE_D.add(Value_G);
    }


    // --------------------------------------------------------------------------------
    //
    // Element testing routines
    //
    // --------------------------------------------------------------------------------

    private class ElementTest extends TestDataProvider {
        public DiffElement elt;
        public String name;
        public String fullName;
        public DiffElement parent;

        private ElementTest(DiffValue elt, DiffValue parent, String name, String fullName) {
            this(elt.getBinding(), parent.getBinding(), name, fullName);
        }

        private ElementTest(DiffElement elt, DiffElement parent, String name, String fullName) {
            super(ElementTest.class);
            this.elt = elt;
            this.name = name;
            this.fullName = fullName;
            this.parent = parent;
        }

        public String toString() {
            return String.format("ElementTest elt=%s name=%s fullName=%s parent=%s",
                    elt.toOneLineString(), name, fullName, parent.getName());
        }
    }

    @DataProvider(name = "elementdata")
    public Object[][] createElementData() {
        new ElementTest(MY_ROOT.getBinding(), DiffElement.ROOT, "MY_ROOT", "MY_ROOT");
        new ElementTest(NODE_C, MY_ROOT, "C", "MY_ROOT.C");
        new ElementTest(NODE_D, MY_ROOT, "D", "MY_ROOT.D");
        new ElementTest(Value_A, MY_ROOT, "A", "MY_ROOT.A");
        new ElementTest(Value_B, MY_ROOT, "B", "MY_ROOT.B");
        new ElementTest(Value_E, NODE_C, "E", "MY_ROOT.C.E");
        new ElementTest(Value_F, NODE_D, "F", "MY_ROOT.D.F");
        new ElementTest(Value_G, NODE_D, "G", "MY_ROOT.D.G");
        return TestDataProvider.getTests(ElementTest.class);
    }

    @Test(enabled = true, dataProvider = "elementdata")
    public void testElementMethods(ElementTest test) {
        Assert.assertNotNull(test.elt.getName());
        Assert.assertNotNull(test.elt.getParent());
        Assert.assertEquals(test.elt.getName(), test.name);
        Assert.assertEquals(test.elt.getParent(), test.parent);
        Assert.assertEquals(test.elt.fullyQualifiedName(), test.fullName);
    }

    // --------------------------------------------------------------------------------
    //
    // DiffValue testing routines
    //
    // --------------------------------------------------------------------------------

    private class LeafTest extends TestDataProvider {
        public DiffValue diffvalue;
        public Object value;

        private LeafTest(DiffValue diffvalue, Object value) {
            super(LeafTest.class);
            this.diffvalue = diffvalue;
            this.value = value;
        }

        public String toString() {
            return String.format("LeafTest diffvalue=%s value=%s", diffvalue.toOneLineString(), value);
        }
    }

    @DataProvider(name = "leafdata")
    public Object[][] createLeafData() {
        new LeafTest(Value_A, "A");
        new LeafTest(Value_B, "B");
        new LeafTest(Value_E, "E");
        new LeafTest(Value_F, "F");
        new LeafTest(Value_G, "G");
        return TestDataProvider.getTests(LeafTest.class);
    }

    @Test(enabled = true, dataProvider = "leafdata")
    public void testLeafMethods(LeafTest test) {
        Assert.assertNotNull(test.diffvalue.getValue());
        Assert.assertEquals(test.diffvalue.getValue(), test.value);
    }

    // --------------------------------------------------------------------------------
    //
    // Node testing routines
    //
    // --------------------------------------------------------------------------------

    private class NodeTest extends TestDataProvider {
        public DiffNode node;
        public Set<String> fields;
        public Set<String> subnodes;
        public Set<String> allNames;

        private NodeTest(DiffNode node, List<String> fields, List<String> subnodes) {
            super(NodeTest.class);
            this.node = node;
            this.fields = new HashSet<>(fields);
            this.subnodes = new HashSet<>(subnodes);
            this.allNames = new HashSet<>(fields);
            allNames.addAll(subnodes);
        }

        public String toString() {
            return String.format("NodeTest node=%s fields=%s subnodes=%s",
                    node.toOneLineString(), fields, subnodes);
        }
    }

    @DataProvider(name = "nodedata")
    public Object[][] createData1() {
        new NodeTest(MY_ROOT, Arrays.asList("A", "B"), Arrays.asList("C", "D"));
        new NodeTest(NODE_C, Arrays.asList("E"), Collections.<String>emptyList());
        new NodeTest(NODE_D, Arrays.asList("F", "G"), Collections.<String>emptyList());
        return TestDataProvider.getTests(NodeTest.class);
    }

    @Test(enabled = true, dataProvider = "nodedata")
    public void testNodeAccessors(NodeTest test) {
        Assert.assertNotNull(test.node.getElements());

        for ( String name : test.allNames ) {
            DiffElement elt = test.node.getElement(name);
            Assert.assertNotNull(elt, "Failed to find field " + elt + " in " + test.node);
            Assert.assertEquals(elt.getName(), name);
            Assert.assertEquals(elt.getValue().isAtomic(), test.fields.contains(name), "Failed atomic/compound expectation: " + test.node);
        }
    }

    // NOTE: add routines are being implicitly tested by the creation of the data structures

    @Test(enabled = true, dataProvider = "nodedata")
    public void testCounts(NodeTest test) {
        Assert.assertEquals(test.node.getElements().size(), test.allNames.size());
        Assert.assertEquals(test.node.getElementNames(), test.allNames);
    }

    // --------------------------------------------------------------------------------
    //
    // fromString testing routines
    //
    // --------------------------------------------------------------------------------

    private class FromStringTest extends TestDataProvider {
        public String string;
        public DiffElement expected;

        private FromStringTest(String string, DiffElement expected) {
            super(FromStringTest.class);
            this.string = string;
            this.expected = expected;
        }

        public String toString() {
            return String.format("FromStringTest string=%s expected=%s", string, expected.toOneLineString());
        }
    }

    @DataProvider(name = "fromstringdata")
    public Object[][] createFromData() {
        new FromStringTest("A=A", Value_A.getBinding());
        new FromStringTest("B=B", Value_B.getBinding());
        new FromStringTest("C=(E=E)", NODE_C.getBinding());
        new FromStringTest("D=(F=F G=G)", NODE_D.getBinding());
        return TestDataProvider.getTests(FromStringTest.class);
    }

    @Test(enabled = true, dataProvider = "fromstringdata")
    public void parseFromString(FromStringTest test) {
        logger.warn("Testing from string: " + test.string);
        DiffElement elt = DiffElement.fromString(test.string);
        Assert.assertEquals(elt.toOneLineString(), test.expected.toOneLineString());
    }
}