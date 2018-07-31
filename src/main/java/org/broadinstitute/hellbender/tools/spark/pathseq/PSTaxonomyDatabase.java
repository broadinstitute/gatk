package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

import java.util.HashMap;
import java.util.Map;

/**
 * Helper class for holding taxonomy data used by ClassifyReads
 */
@DefaultSerializer(PSTaxonomyDatabase.Serializer.class)
public class PSTaxonomyDatabase {
    public final PSTree tree;
    public final Map<String, Integer> accessionToTaxId; //Reference contig name to taxonomic ID

    public PSTaxonomyDatabase(final PSTree tree, final Map<String, Integer> map) {
        this.tree = tree;
        this.accessionToTaxId = map;
    }

    private PSTaxonomyDatabase(final Kryo kryo, final Input input) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        tree = kryo.readObject(input, PSTree.class);
        final int setSize = input.readInt();
        accessionToTaxId = new HashMap<>(setSize);
        for (int i = 0; i < setSize; i++) {
            final String key = input.readString();
            final String value = input.readString();
            accessionToTaxId.put(key, Integer.valueOf(value));
        }

        kryo.setReferences(oldReferences);
    }

    protected void serialize(final Kryo kryo, final Output output) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        kryo.writeObject(output, tree);
        output.writeInt(accessionToTaxId.size());
        for (final String key : accessionToTaxId.keySet()) {
            output.writeString(key);
            output.writeString(String.valueOf(accessionToTaxId.get(key)));
        }

        kryo.setReferences(oldReferences);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PSTaxonomyDatabase> {
        @Override
        public void write(final Kryo kryo, final Output output, final PSTaxonomyDatabase taxonomyDatabase) {
            taxonomyDatabase.serialize(kryo, output);
        }

        @Override
        public PSTaxonomyDatabase read(final Kryo kryo, final Input input,
                           final Class<PSTaxonomyDatabase> klass) {
            return new PSTaxonomyDatabase(kryo, input);
        }
    }
}
