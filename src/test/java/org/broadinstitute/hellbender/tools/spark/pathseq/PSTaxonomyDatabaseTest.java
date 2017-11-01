package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.HashMap;
import java.util.Map;

public class PSTaxonomyDatabaseTest extends GATKBaseTest {

    @Test
    public void testSerializeDeserialize() {

        final Map<String,Integer> accessionToTaxMap = new HashMap<>();
        accessionToTaxMap.put("A",1);
        accessionToTaxMap.put("B",1);
        accessionToTaxMap.put("C",2);
        accessionToTaxMap.put("D",3);

        final PSTree tree = new PSTree(1);
        tree.addNode(2, "node2", 1, 0, "genus");
        tree.addNode(3, "node3", 1, 0, "genus");
        tree.addNode(4, "node4", 2, 200, "species");
        tree.addNode(5, "node5", 2, 300, "species");
        tree.addNode(6, "node6", 3, 100, "species");

        final PSTaxonomyDatabase taxonomyDatabase = new PSTaxonomyDatabase(tree, accessionToTaxMap);

        final Kryo kryo = new Kryo();
        final Output output = new Output(new ByteArrayOutputStream());
        kryo.writeObject(output, taxonomyDatabase);

        final Input input = new Input(new ByteArrayInputStream(output.getBuffer()));
        final PSTaxonomyDatabase taxonomyDatabaseTest = kryo.readObject(input, PSTaxonomyDatabase.class);

        Assert.assertEquals(taxonomyDatabaseTest.tree, taxonomyDatabase.tree);
        Assert.assertEquals(taxonomyDatabaseTest.accessionToTaxId, taxonomyDatabase.accessionToTaxId);
    }

}