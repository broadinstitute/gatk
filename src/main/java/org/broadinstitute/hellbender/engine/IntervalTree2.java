package org.broadinstitute.hellbender.engine;

import it.unimi.dsi.fastutil.ints.IntArrayFIFOQueue;

import java.util.*;

/**
 * @param <E> type of the element.
 */
public class IntervalTree2<E> {



    @FunctionalInterface
    private interface IntervalComparator {
        int compare(int s1, int e1, int s2, int e2);
    }

    /**
     * Comparator used to sort intervals within the tree.
     */
    private final IntervalComparator comparator;

    /**
     * Holds a reference to the root of the tree or {@link #NIL}
     * if there is no root (i.e. the tree is empty).
     */
    private Node root;

    /**
     * Sentinel or NIL value use to indicate the lack of a node.
     */
    private final Node NIL = new Node();

    /**
     * Tree modification time-stamp use to detect concurrent modifications
     */
    private volatile long stamp;

    /**
     * Number of intervals in the tree.
     */
    private int numberOfIntervals;

    /**
     * Number of elements in the tree.
     */
    private int numberOfElements;

    public class Node
    {

        private static final boolean BLACK = false;
        private static final boolean RED = true;

        private boolean color;
        private int start, end, max, min;
        private Node parent;
        private Node left;
        private Node right;
        private List<E> elements;

        private Node() {
            parent = this;
            start = end = max = min = -1;
            left = this;
            right = this;
            this.elements = Collections.emptyList();
            color = BLACK;
        }

        public int start() {
            return start;
        }

        public int end() {
            return end;
        }

        public List<E> elements() {
            return elements;
        }

        private Node(final Node parent, final int start, final int end, final E elem) {
            this.start = min = start;
            this.end = max = end;
            this.elements = Collections.singletonList(elem);
            this.parent = parent;
            this.left = this.right = NIL;
            this.color = RED;
        }

        private Node grandParent() {
            return parent.parent;
        }

        private Node firstDescendant() {
            Node result = this;
            while (result.left != NIL) {
                result = result.left;
            }
            return result;
        }

        private Node lastDescendant() {
            Node result = this;
            while (result.right != NIL) {
                result = result.right;
            }
            return result;
        }

        private Node nextNode() {
            if (right != NIL) {
                return right.firstDescendant();
            } else {
                Node result = this;
                while (result != NIL) {
                    if (result.parent.left == result) {
                        return result.parent;
                    } else {
                        result = result.parent;
                    }
                }
                return result;
            }
        }

        private Node previousNode() {
            if (left != NIL) {
                return left.lastDescendant();
            } else {
                Node result = this;
                while (result != NIL) {
                    if (result.parent.right == result) {
                        return result.parent;
                    } else {
                        result = result.parent;
                    }
                }
                return result;
            }
        }

        private Node uncle() {
            return parent.left == this ? parent.right : parent.left;
        }

        private void add(final E elem) {
            if (elements.size() == 1) {
                final List<E> newElements = new ArrayList<>(10);
                newElements.add(elements.get(0));
                elements = newElements;
            }
            elements.add(elem);
        }

        private void recalculateMinMaxFromOffspring() {
            min = start;
            max = end;
            if (left != NIL) {
                min = Math.min(min, left.min);
                max = Math.max(max, left.max);
            }
            if (right != NIL) {
                min = Math.min(min, right.min);
                max = Math.max(max, right.max);
            }
        }

        private boolean replaceOffspring(final Node a, final Node b) {
            if (left == a) {
                left = b;
                left.parent = this;
                return true;
            } else if (right == a) {
                right = b;
                right.parent = this;
                return true;
            } else {
                return false;
            }
        }

        private void adoptRight(final Node newChild) {
            right = newChild;
            if (newChild != NIL) {
                newChild.parent = this;
            }
        }

        private void adoptLeft(final Node newChild) {
            left = newChild;
            if (newChild != NIL) {
                newChild.parent = this;
            }
        }

        private void makeBlack() {
            color = BLACK;
        }

        private void makeRed() {
            color = RED;
        }

        public boolean isBlack() {
            return color == BLACK;
        }

        public boolean isRed() {
            return color == RED;
        }

        public boolean isRoot() {
            return this == root;
        }

        public Node sibling() {
            if (parent.isRoot()) {
                return NIL;
            } else {
                return parent.left == this ? parent.right : parent.left;
            }
        }

        public boolean isLeftChild() {
            return parent.left == this;
        }

        public boolean isRightChild() {
            return parent.right == this;
        }

        private boolean hasTwoOffspring() {
            return left != NIL && right != NIL;
        }

        public String toString() {
            if (this == NIL) {
                return "NIL " + (color == BLACK ? "BLACK" : "RED");
            } else {
                return String.format("[%d, %d] <<%d, %d>>", start, end, min, max)  + (color == BLACK ? "BLACK" : "RED");
            }
        }
    }

    public IntervalTree2() {
        this((s1, e1, s2, e2) -> {
           int cmp = s1 - s2;
           if (cmp != 0) {
               return cmp;
           } else {
               return e1 - e2;
           }
        });
    }

    public IntervalTree2(final IntervalComparator comparator) {
        this.comparator = comparator;
        stamp = Long.MIN_VALUE;
        numberOfIntervals = 0;
        root = NIL;
    }

    public int numberOfIntervals() {
        return numberOfIntervals;
    }

    public Node insert(final int start, final int end, final E elem) {
        ++stamp;
        final Node result = binaryTreeInsert(start, end, elem);
        rebalanceTree(result);
        numberOfIntervals++;
        return result;
    }

    public List<E> overlappingElements(final int start, final int end) {
        final Deque<Node> candidateNodes = new ArrayDeque<>(numberOfIntervals());
        candidateNodes.add(root);
        final ArrayList<E> result = new ArrayList<>(numberOfIntervals());
        while (!candidateNodes.isEmpty()) {
            final Node next = candidateNodes.remove();
            if (next != NIL) {
                if (next.min <= end && next.max >= start) {
                    candidateNodes.add(next.left);
                    candidateNodes.add(next.right);
                    if (next.start <= end && next.end >= start) {
                        result.addAll(next.elements);
                    }
                }
            }
        }
        return result;
    }

    public IntervalTreeIterator2<E> iterator() {
        final IntervalTree2<E> tree = this;

        return new IntervalTreeIterator2<E>() {

            private long stamp = tree.stamp;
            private Node next = tree.root.firstDescendant();
            private Node prev = NIL;
            private Node last = prev;

            @Override
            public Node next() {
                last = prev = peekNext();
                next = next.nextNode();
                return prev; // == next before calling this method.
            }

            public boolean hasNext() {
                return next != NIL;
            }

            @Override
            public Node peekNext() {
                if (stamp == tree.stamp && next != NIL) {
                    return next;
                } else if (next == NIL) {
                    throw new NoSuchElementException();
                } else {
                    throw new ConcurrentModificationException();
                }
            }

            @Override
            public Node previous() {
                last = next = peekPrevious();
                prev = prev.previousNode();
                return next; // == prev before calling thi method.
            }

            public boolean hasPrevious() {
                return prev != NIL;
            }

            public Node peekPrevious() {
                if (stamp == tree.stamp && prev != NIL) {
                    return prev;
                } else if (prev == NIL) {
                    throw new NoSuchElementException();
                } else {
                    throw new ConcurrentModificationException();
                }
            }


            @Override
            @SuppressWarnings("unchecked")
            public IntervalTreeIterator2<E> remove() {
                if (last == NIL) {
                    throw new IllegalStateException("must iterate over some element before removing");
                }
                final Node surrogate;
                if (last.hasTwoOffspring()) {
                    if (prev == last) {
                        surrogate = next;
                        next = last;
                        prev = prev.previousNode();
                    } else {
                        surrogate = prev;
                        prev = last;
                        next = next.nextNode();
                    }
                    last.start = surrogate.start;
                    last.end = surrogate.end;
                    last.elements = surrogate.elements;
                    last.recalculateMinMaxFromOffspring();
                    last = NIL;
                } else {
                    surrogate = last;
                    if (last == prev) {
                        prev = prev.previousNode();
                    } else {
                        next = next.nextNode();
                    }
                    last = NIL;
                }
                final Node successor = surrogate.left == NIL ?
                        surrogate.right : surrogate.left;
                if (surrogate == root) {
                    root = successor;
                    successor.parent = NIL;
                } else {
                    surrogate.parent.replaceOffspring(surrogate, successor);
                    surrogate.parent.recalculateMinMaxFromOffspring();
                    if (surrogate.isBlack()) {
                        repairAfterDelete(successor);
                    }
                }
                return this;
            }

            private void repairAfterDelete(Node node) {
                {
                    if (node.isRed()) {
                        node.makeBlack();
                        return;
                    }
                    while (!node.isRoot()) {
                        final Node parent = node.parent;
                        final boolean nodeIsLeft = node == parent.left;
                        final Node sibling = nodeIsLeft ? parent.right : parent.left;
                        if (sibling.isRed()) {
                            parent.makeRed();
                            sibling.makeBlack();
                            rotate(parent, nodeIsLeft);
                        }
                        final boolean siblingFamilyIsBlack = sibling.isBlack()
                                && sibling.left.isBlack()
                                && sibling.right.isBlack();
                        if (parent.isBlack()) {
                            if (siblingFamilyIsBlack) {
                                sibling.makeRed();
                                node = parent;
                                continue; // back to case 1 with node == parent.
                            }
                        } else {
                            if (siblingFamilyIsBlack) {
                                sibling.makeRed();
                                parent.isBlack();
                                break;
                            }
                        }
                        if (sibling.isBlack()) {
                            if (nodeIsLeft && sibling.right.isBlack() && sibling.left.isRed()) {
                                sibling.makeRed();
                                sibling.left.makeBlack();
                                rotateRight(sibling);
                            } else if (!nodeIsLeft && sibling.left.isBlack() && sibling.right.isRed()) {
                                sibling.makeRed();
                                sibling.right.makeBlack();
                                rotateLeft(sibling);
                            }
                        }
                        // case 6:
                        sibling.color = parent.color;
                        parent.makeBlack();
                        if (nodeIsLeft) {
                            sibling.right.makeBlack();
                            rotateLeft(parent);
                        } else {
                            sibling.left.makeBlack();
                            rotateRight(parent);
                        }
                    }
                }
            }


            @Override
            public IntervalTreeIterator2<E> seek(final int start, final int end) {
                // seek always resets last (i.e. remove is not allowed after a seek)
                last = NIL;
                // are we there already? quick check whether the next is in fact where the seek would go.
                if (next != NIL) {
                    if (comparator.compare(next.start, next.end, start, end) >= 0
                            && (prev == NIL || comparator.compare(prev.start, prev.end, start, end) < 0)) {
                        return this;
                    }
                }
                Node current;
                next = prev = NIL;
                current = root;

                while (current != NIL) {
                    final int cmp = comparator.compare(current.start, current.end, start, end);
                    if (cmp < 0) {
                        current = (prev = current).right;
                    } else if (cmp > 0) {
                        current = (next = current).left;
                    } else { // cmp == 0
                        // if this node has left descendants the
                        // prev is the last of them:
                        if ((next = current).left != NIL) {
                            prev = next.left.lastDescendant();
                        }
                        break;
                    }
                }
                return this;
            }
        };
    }

    private Node binaryTreeInsert(final int start, final int end, final E elem) {
        Node parent;
        Node current = parent = root;
        boolean wentLeft = false;
        while (current != NIL) {
            current.max = Math.max(current.max, end);
            current.min = Math.min(current.min, start);
            final int cmp = comparator.compare(current.start, current.end, start, end);
            if (cmp < 0) {
                wentLeft = false;
                parent = current;
                current = current.right;
            } else if (cmp > 0) {
                wentLeft = true;
                parent = current;
                current = current.left;
            } else { // { if (end == current.end) {
                current.add(elem);
                return current;
            }
        }
        current = new Node(parent, start, end, elem);
        if (parent != NIL) {
            if (wentLeft) {
                parent.left = current;
            } else {
                parent.right = current;
            }
        } else {
            root = current;
        }
        return current;
    }

    private void rebalanceTree(final Node start) {
        Node current = start;
        do {
            if (current.isRoot()) { // current is root.
                current.makeBlack();
                break;
            } else if (current.parent.isBlack()) {
                break;
            } else {
                final Node uncle = current.uncle();
                if (uncle.isRed()) {
                    current.parent.color = uncle.color = Node.BLACK;
                    // grand parent cannot be NULL:
                    final Node grandParent = current.grandParent();
                    grandParent.makeRed();
                    current = grandParent;
                } else {
                    repairTreeUsingRotations(current);
                }
            }
        } while (true);
    }

    private void repairTreeUsingRotations(final Node current) {
        // guaranteed that parent and grant-parent are not null nor NIL.
        final Node parent = current.parent;
        final Node grandParent = current.parent.parent;
        Node secondStepCurrent = parent;
        if (current == parent.right) {
            if (grandParent.left == parent) {
                rotateLeft(parent);
                secondStepCurrent = current;
            }
        } else if (current == parent.left) {
            if (grandParent.right == parent) {
                rotateRight(parent);
                secondStepCurrent = current;
            }
        }
        if (secondStepCurrent == grandParent.left) {
            rotateRight(grandParent);
        } else {
            rotateLeft(grandParent);
        }
        grandParent.makeRed();
        secondStepCurrent.makeBlack();
    }

    /**
     * Assumes that child is either the left or right offspring
     * for parent.
     * @param parent the parent node
     * @param child the offspring to rotate on.
     */
    private void rotateOn(final Node parent, final Node child) {
        if (parent.left == child) {
            rotateLeft(parent);
        } else {
            rotateRight(parent);
        }
    }

    private void rotateLeft(final Node parent) {
        final Node child = parent.right;
        final Node grandChild = child.left;
        parent.adoptRight(grandChild);
        child.adoptLeft(parent);
        if (parent.parent != NIL) {
            parent.parent.replaceOffspring(parent, child);
        } else {
            child.parent = NIL;
        }
        // child's new min-max is parent's old min-max:
        child.max = parent.max;
        child.min = parent.min;
        // we recalculate the parent new min-max.
        parent.recalculateMinMaxFromOffspring();
    }

    private void rotate(final Node node, final boolean left) {
        if (left) {
            rotateLeft(node);
        } else {
            rotateRight(node);
        }
    }

    private void rotateRight(final Node parent) {
        final Node child = parent.left;
        final Node grandChild = child.right;
        parent.adoptLeft(grandChild);
        child.adoptRight(parent);
        if (parent.parent != NIL) {
            parent.parent.replaceOffspring(parent, child);
        } else {
            child.parent = NIL;
        }
        // child's new min-max is parent's old min-max:
        child.max = parent.max;
        child.min = parent.min;
        // we recalculate the parent new min-max.
        parent.recalculateMinMaxFromOffspring();
    }

    /**
     * Compute the inbalance of the tree. This is a value from 0 to 100 where 0 corresponds to a fully balanced
     * tree and 1 with the worst case scenario for a RB-tree.
     * @return 0 to 1.
     */
    public double inbalance() {
        final Deque<Node> toProcess = new ArrayDeque<>(10);
        final IntArrayFIFOQueue heights = new IntArrayFIFOQueue(10);
        toProcess.push(root);
        heights.enqueue(0);
        int min = Integer.MAX_VALUE;
        int max = 0;
        while (!toProcess.isEmpty()) {
            final Node top = toProcess.pop();
            final int height = heights.dequeueInt();
            if (top == NIL) {
                min = Math.min(min, height);
                max = Math.max(max, height);
            } else {
                toProcess.add(top.left);
                toProcess.add(top.right);
                heights.enqueue(height + 1);
                heights.enqueue(height + 1);
            }
        }
        return max == 0 ? 0 : (max - min + 1.0) / (min + 1.0);
    }
}
