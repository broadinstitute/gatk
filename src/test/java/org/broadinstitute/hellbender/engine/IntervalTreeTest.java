package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.IntervalTree;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class IntervalTreeTest {

    static class TestElement {
        public final int start;
        public final int end;

        TestElement(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }

    private <E> List<E> mergeLists(final List<E> a, final List<E> b) {
        if (a instanceof ArrayList) {
            a.addAll(b);
            return a;
        } else if (b instanceof ArrayList) {
            b.addAll(a);
            return b;
        } else {
            final List<E> c = new ArrayList<>();
            c.addAll(a);
            c.addAll(b);
            return c;
        }
    }

    private int compare(final TestElement a, final TestElement b) {
        int cmp = a.start - b.start;
        if (cmp != 0) {
            return cmp;
        } else {
            return a.end - b.end;
        }
    }

    @Test(dataProvider = "randomElementLists")
    public void testList(final List<TestElement> testElements) {
        final IntervalTree<List<TestElement>> subject = new IntervalTree<>();
        for (final TestElement te : testElements) {
            subject.merge(te.start, te.end, Collections.singletonList(te), this::mergeLists);
        }
        final List<TestElement> sorted = testElements.stream().sorted(this::compare).collect(Collectors.toList());
        final Iterator<IntervalTree.Node<List<TestElement>>> it = subject.iterator();
        for (int i = 0; i < testElements.size();) {
            Assert.assertTrue(it.hasNext());
            final IntervalTree.Node<List<TestElement>> node = it.next();
            final List<TestElement> payload = node.getValue();
            for (int j = 0; j < payload.size(); j++) {
                TestElement element = sorted.get(i++);
                TestElement nodeElement = payload.get(j);
                Assert.assertEquals(nodeElement.start, node.getStart());
                Assert.assertEquals(nodeElement.end, node.getEnd());
                Assert.assertEquals(nodeElement.start, element.start);
                Assert.assertEquals(nodeElement.end, element.end);
            }
        }
    }


    @Test(dataProvider = "randomElementLists")
    public void testRemove(final List<TestElement> testElements) {
        final IntervalTree<List<TestElement>> subject = new IntervalTree<>();
        for (final TestElement te : testElements) {
            subject.merge(te.start, te.end, Collections.singletonList(te), this::mergeLists);
        }

        final Random rdn = new Random((13 + testElements.size() * 31) * 31 + (testElements.size() > 0? testElements.get(0).start : 0));
        final Set<TestElement> toRemove = testElements.stream().filter(x -> rdn.nextDouble() > 0.5).collect(Collectors.toSet());
        final List<TestElement> removeList = new ArrayList<>(toRemove);
        Collections.shuffle(removeList, rdn);
        for (final TestElement te : removeList) {
            final IntervalTree.Node<List<TestElement>> node = subject.find(te.start , te.end);
            if (node.getValue().size() == 1) {
                Assert.assertSame(node.getValue().get(0), te);
                subject.remove(te.start, te.end);
            } else {
                Assert.assertTrue(node.getValue().remove(te));
            }
        }
        final List<TestElement> sorted = testElements.stream().sorted(this::compare).filter(e -> !toRemove.contains(e)).collect(Collectors.toList());
        final Iterator<IntervalTree.Node<List<TestElement>>> it = subject.iterator();
        for (int i = 0; i < testElements.size() - toRemove.size();) {
            Assert.assertTrue(it.hasNext());
            final IntervalTree.Node<List<TestElement>> node = it.next();
            final List<TestElement> payload = node.getValue();
            for (int j = 0; j < payload.size(); j++) {
                TestElement element = sorted.get(i++);
                TestElement nodeElement = payload.get(j);
                Assert.assertEquals(nodeElement.start, node.getStart());
                Assert.assertEquals(nodeElement.end, node.getEnd());
                Assert.assertEquals(nodeElement.start, element.start);
                Assert.assertEquals(nodeElement.end, element.end);
            }
        }
        Assert.assertFalse(it.hasNext());
    }

    @Test(dataProvider = "randomElementLists")
    public void testListNew(final List<TestElement> testElements) {
        final IntervalTree2<TestElement> subject = new IntervalTree2<>();
        for (final TestElement te : testElements) {
            subject.insert(te.start, te.end, te);
        }
        final List<TestElement> sorted = testElements.stream().sorted(this::compare).collect(Collectors.toList());
        final IntervalTreeIterator2<TestElement> it = subject.iterator();
        for (int i = 0; i < testElements.size();) {
            Assert.assertTrue(it.hasNext());
            final IntervalTree2<TestElement>.Node node = it.next();
            final List<TestElement> payload = node.elements();
            for (int j = 0; j < payload.size(); j++) {
                if (payload.size() > 1) {
                }
                TestElement element = sorted.get(i++);
                TestElement nodeElement = payload.get(j);
                Assert.assertEquals(nodeElement.start, node.start());
                Assert.assertEquals(nodeElement.end, node.end());
                Assert.assertEquals(nodeElement.start, element.start, "i=" + i + "j=" + j);
                Assert.assertEquals(nodeElement.end, element.end);
            }
        }
    }

    @Test(dataProvider = "randomElementLists")
    public void testRemoveNew(final List<TestElement> testElements) {
        final IntervalTree2<TestElement> subject = new IntervalTree2<>();
        for (final TestElement te : testElements) {
            subject.insert(te.start, te.end, te);
        }
        double inbalance = subject.inbalance();
        if (inbalance >= 2.5) {
            System.err.println("whatta");
        }
        Assert.assertTrue(inbalance <= 1.0, "" + inbalance);
        Assert.assertTrue(inbalance >= 0.0, "" + inbalance);

        final Random rdn = new Random((13 + testElements.size() * 31) * 31 + (testElements.size() > 0? testElements.get(0).start : 0));
        final Set<TestElement> toRemove = testElements.stream().filter(x -> rdn.nextDouble() > 0.5).collect(Collectors.toSet());
        final List<TestElement> removeList = new ArrayList<>(toRemove);
        Collections.shuffle(removeList, rdn);
        IntervalTreeIterator2<TestElement> it = subject.iterator();
        for (final TestElement te : removeList) {
            final IntervalTree2<TestElement>.Node node = it.seek(te.start, te.end).next();
            if (node.elements().size() == 1) {
                Assert.assertSame(node.elements().get(0), te);
                it.remove();
            } else {
                Assert.assertTrue(node.elements().remove(te));
            }
        }
        final List<TestElement> sorted = testElements.stream().sorted(this::compare).filter(e -> !toRemove.contains(e)).collect(Collectors.toList());
        it = subject.iterator();
        for (int i = 0; i < sorted.size();) {
            Assert.assertTrue(it.hasNext());
            final IntervalTree2<TestElement>.Node node = it.next();
            final List<TestElement> payload = node.elements();
            for (int j = 0; j < payload.size(); j++) {
                TestElement element = sorted.get(i++);
                TestElement nodeElement = payload.get(j);
                Assert.assertEquals(nodeElement.start, node.start());
                Assert.assertEquals(nodeElement.end, node.end());
                Assert.assertEquals(nodeElement.start, element.start, "i=" + i + "j=" + j);
                Assert.assertEquals(nodeElement.end, element.end);
            }
        }
        Assert.assertFalse(it.hasNext());
        inbalance = subject.inbalance();
        if (inbalance >= 2.5) {
            System.err.println("whatta");
        }
        Assert.assertTrue(inbalance <= 1.0, "" + inbalance);
        Assert.assertTrue(inbalance >= 0.0, "" + inbalance);
    }

    @DataProvider
    private Object[][] randomElementLists() {
        try {
            final List<Object[]> result = new ArrayList<>(100);
            final Random rdn = new Random(13);
            for (int i = 0; i < 100; i++) {
                result.add(new Object[]{randomElements(rdn, 100, 50, 100, 50, i)});
            }
            return result.toArray(new Object[1000][]);
        } catch (final Throwable th) {
            th.printStackTrace(System.err);
            throw th;
        }
    }

    private List<TestElement> randomElements(final Random rdn, int startMean, int sizeMean, double startDev, double sizeDev, int count) {

        final List<TestElement> result = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            final int start = (int) (rdn.nextDouble() * startDev + startMean);
            final int size = Math.max(1, (int) (rdn.nextDouble() * sizeDev + sizeMean));
            result.add(new TestElement(start, start + size - 1));
        }
        return result;
    }
}
