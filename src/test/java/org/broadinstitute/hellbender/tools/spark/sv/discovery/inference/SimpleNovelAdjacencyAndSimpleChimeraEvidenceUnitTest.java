package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;

public class SimpleNovelAdjacencyAndSimpleChimeraEvidenceUnitTest extends AssemblyBasedSVDiscoveryBaseTest {


    @DataProvider(name = "forKryoSerializationAndHashCode")
    private Object[][] forKryoSerializationAndHashCode() {
        final List<Object[]> data = new ArrayList<>();

        for (final Tuple2<? extends AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera, ? extends AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera>
                pair : getAllTestDataInPairs()) {
            final NovelAdjacencyAndAltHaplotype biPathBubble = pair._1.expectedNovelAdjacencyAndAltSeq;
            final SimpleChimera forwardRep = new SimpleChimera(pair._1.firstAlignment, pair._1.secondAlignment, Collections.emptyList(),
                    pair._1.evidenceAssemblyContigName, NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME,
                    pair._1.getAppropriateDictionary());
            final SimpleChimera reverseRep = new SimpleChimera(pair._2.firstAlignment, pair._2.secondAlignment, Collections.emptyList(),
                    pair._2.evidenceAssemblyContigName, NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME,
                    pair._2.getAppropriateDictionary());
            final List<SimpleChimera> evidence = Arrays.asList(forwardRep, reverseRep);
            data.add(new Object[]{new SimpleNovelAdjacencyAndChimericAlignmentEvidence(biPathBubble, evidence)});
        }
        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forKryoSerializationAndHashCode")
    public void testKryoSerializerAndHashCode(final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence) throws IOException {
        try (final ByteArrayOutputStream bos = new ByteArrayOutputStream()) {

            final Output out = new Output(bos);
            final Kryo kryo = new Kryo();
            kryo.writeClassAndObject(out, simpleNovelAdjacencyAndChimericAlignmentEvidence);
            out.flush();

            try ( final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray()) ) {
                final Input in = new Input(bis);
                @SuppressWarnings("unchecked")
                final SimpleNovelAdjacencyAndChimericAlignmentEvidence roundTrip = (SimpleNovelAdjacencyAndChimericAlignmentEvidence) kryo.readClassAndObject(in);
                Assert.assertEquals(roundTrip, simpleNovelAdjacencyAndChimericAlignmentEvidence);
                Assert.assertEquals(roundTrip.hashCode(), simpleNovelAdjacencyAndChimericAlignmentEvidence.hashCode());
            }
        }
    }
}
