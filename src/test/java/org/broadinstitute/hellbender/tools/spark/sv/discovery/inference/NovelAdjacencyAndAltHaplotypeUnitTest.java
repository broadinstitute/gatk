package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class NovelAdjacencyAndAltHaplotypeUnitTest extends AssemblyBasedSVDiscoveryBaseTest {

    @DataProvider(name = "forKryoSerializationAndHashCode")
    private Object[][] forKryoSerializationAndHashCode() {
        final List<Object[]> data = new ArrayList<>();
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera testData : getAllTestData()) {
            data.add(new Object[]{testData.expectedNovelAdjacencyAndAltSeq});
        }
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forKryoSerializationAndHashCode")
    public void testKryoSerializerAndHashCode(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) throws IOException {
        try (final ByteArrayOutputStream bos = new ByteArrayOutputStream()) {

            final Output out = new Output(bos);
            final Kryo kryo = new Kryo();
            kryo.writeClassAndObject(out, novelAdjacencyAndAltHaplotype);
            out.flush();

            try ( final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray()) ) {
                final Input in = new Input(bis);
                @SuppressWarnings("unchecked")
                final NovelAdjacencyAndAltHaplotype roundTrip = (NovelAdjacencyAndAltHaplotype) kryo.readClassAndObject(in);
                Assert.assertEquals(roundTrip, novelAdjacencyAndAltHaplotype);
                Assert.assertEquals(roundTrip.hashCode(), novelAdjacencyAndAltHaplotype.hashCode());
            }
        }
    }

    @DataProvider(name = "forToSimpleOrBNDTypes")
    private Object[][] forToSimpleOrBNDTypes() {
        final List<Object[]> data = new ArrayList<>(20);

        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera assemblyBasedSVDiscoveryTestDataForSimpleChimera : getAllTestData()) {
            data.add(new Object[]{assemblyBasedSVDiscoveryTestDataForSimpleChimera.expectedNovelAdjacencyAndAltSeq,
                                  assemblyBasedSVDiscoveryTestDataForSimpleChimera.expectedSvTypes,
                                  assemblyBasedSVDiscoveryTestDataForSimpleChimera.getAppropriateRef(),
                                  assemblyBasedSVDiscoveryTestDataForSimpleChimera.getAppropriateDictionary()}
            );
        }

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forToSimpleOrBNDTypes")
    public void testToSimpleOrBNDTypes(final NovelAdjacencyAndAltHaplotype breakpoints,
                                       final List<SvType> expectedInferredTypes,
                                       final ReferenceMultiSparkSource ref,
                                       final SAMSequenceDictionary dict) {
        final List<SvType> actual = breakpoints.toSimpleOrBNDTypes(ref, dict);
        Assert.assertEquals(actual, expectedInferredTypes);
        for (int i = 0; i < actual.size(); ++i) {
            Assert.assertEquals(actual.get(i).hashCode(), expectedInferredTypes.get(i).hashCode());
        }
    }
}