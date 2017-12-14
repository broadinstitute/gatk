package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.Collection;
import java.util.HashSet;
import java.util.Objects;

/**
 * Node class for PSTree
 */
@DefaultSerializer(PSTreeNode.Serializer.class)
public class PSTreeNode {

    private String name = null;
    private String rank = null;
    private int parent = 0;
    private long length = 0;
    private Collection<Integer> children;

    public PSTreeNode() {
        children = new HashSet<>();
    }

    private PSTreeNode(final Kryo kryo, final Input input) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        name = input.readString();
        rank = input.readString();
        parent = input.readInt();
        length = input.readLong();
        final int numChildren = input.readInt();
        children = new HashSet<>(numChildren);
        for (int i = 0; i < numChildren; i++) {
            children.add(Integer.valueOf(input.readString()));
        }

        kryo.setReferences(oldReferences);
    }

    private void serialize(final Kryo kryo, final Output output) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        output.writeString(name);
        output.writeString(rank);
        output.writeInt(parent);
        output.writeLong(length);
        output.writeInt(children.size());
        for (final int child : children) {
            output.writeString(String.valueOf(child));
        }

        kryo.setReferences(oldReferences);
    }

    public String getName() {
        return name;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public String getRank() {
        return rank;
    }

    public void setRank(final String rank) {
        this.rank = rank;
    }

    public int getParent() {
        return parent;
    }

    public void setParent(final int parent) {
        this.parent = parent;
    }

    public long getLength() {
        return length;
    }

    public void setLength(final long length) {
        this.length = length;
    }

    public Collection<Integer> getChildren() {
        return children;
    }

    public void addChild(final int child) {
        this.children.add(child);
    }

    public boolean removeChild(final int child) {
        if (this.children.contains(child)) {
            this.children.remove(child);
            return true;
        }
        return false;
    }

    @Override
    public String toString() {
        return name + "," + rank + "," + parent + "," + length + "," + children;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final PSTreeNode that = (PSTreeNode) o;

        if (length != that.length) return false;
        if (name != null ? !name.equals(that.name) : that.name != null) return false;
        if (rank != null ? !rank.equals(that.rank) : that.rank != null) return false;
        if (parent != that.parent) return false;
        return children.equals(that.children);
    }

    @Override
    public int hashCode() {
        int result = Objects.hashCode(name);
        result = 31 * result + Objects.hashCode(rank);
        result = 31 * result + parent;
        result = 31 * result + (int) (length ^ (length >>> 32));
        result = 31 * result + children.hashCode();
        return result;
    }

    public PSTreeNode copy() {
        final PSTreeNode n = new PSTreeNode();
        n.name = this.name;
        n.rank = this.rank;
        n.parent = this.parent;
        n.length = this.length;
        n.children = new HashSet<>(this.children);
        return n;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PSTreeNode> {
        @Override
        public void write(final Kryo kryo, final Output output, final PSTreeNode psTreeNode) {
            psTreeNode.serialize(kryo, output);
        }

        @Override
        public PSTreeNode read(final Kryo kryo, final Input input, final Class<PSTreeNode> klass) {
            return new PSTreeNode(kryo, input);
        }
    }
}