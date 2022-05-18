package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.generic.GenericRecordBuilder;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

class ExtractCohortEngineTest extends CommandLineProgramTest {
    // Refactoring Opportunity -- we dynamically read the schema from the BQ tables when reading with
    // the Read/Write API, so we never declare it explicitly. However, this is because we use a GenericRecord
    // for Avro instead of a custom class.. if we did that, we wouldn't need this and might even have better validation
    // when reading
    private Schema vetSchema = SchemaBuilder.record("vet")
            .namespace("org.broadinstitute.dsp")
            .fields()
            .optionalLong("sample_id")
            .optionalLong("location")
            .optionalString("ref")
            .optionalString("alt")
            .optionalString("AS_RAW_MQ")
            .optionalString("AS_RAW_MQRankSum")
            .optionalString("QUALapprox")
            .optionalString("AS_QUALapprox")
            .optionalString("AS_RAW_ReadPosRankSum")
            .optionalString("AS_SB_TABLE")
            .optionalString("AS_VarDP")
            .optionalString("call_GT")
            .optionalString("call_AD")
            .optionalLong("call_GQ")
            .optionalString("call_PGT")
            .optionalString("call_PID")
            .optionalString("call_PL")
            .endRecord();

    // Refactoring Opportunity -- this enum is statically stateful...
    static {
        ChromosomeEnum.setRefVersion("hg38");
    }

    private static final ReferenceDataSource reference = new ReferenceFileSource(IOUtils.getPath(b38_reference_20_21));

    private ExtractCohortEngine getMinimalEngine(SimpleInterval interval, Map<Long, String> sampleIdToName) {
        VCFHeader vcfHeader = new VCFHeader();
        VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Collections.emptyList(), null, Collections.emptyList(), false, false);

        // The ExtractCohortEngine is fairly stateful, so even though we are going to unit test
        // individual methods we need to construct the engine so that those methods can pull instance state
        //
        // Refactoring Opportunity: this is a bit ugly since we have to specify things we know we won't need for our test,
        // however, the general philosophy in GATK seems to be to create this monolithic engine statefully and then use it.
        return new ExtractCohortEngine(null,
                null,
                vcfHeader,
                annotationEngine,
                reference,
                sampleIdToName,
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                Arrays.asList(interval),
                SchemaUtils.encodeLocation(interval.getContig(), interval.getStart()),
                SchemaUtils.encodeLocation(interval.getContig(), interval.getEnd()),
                null,
                null,
                1000000,
                false,
                null,
                null,
                null,
                null,
                false,
                false,
                ExtractCohort.VQSLODFilteringType.NONE,
                false,
                GQStateEnum.SIXTY,
                false);
    }

    @Test
    public void testSpanningDeletion() throws Exception {
        Map<Long, String> sampleIdToName = new HashMap<>();
        sampleIdToName.put(1L, "sample_1");
        sampleIdToName.put(2L, "sample_2");
        sampleIdToName.put(3L, "sample_3");

        ExtractCohortEngine engine = getMinimalEngine(new SimpleInterval("chr20", 100000, 100010), sampleIdToName);

        // we have 3 samples (sample1, sample2, sample3) and two locations of interest  CTTGAGCATG
        // chr20:100000 -- alt is a 10bp deletion, sample 1 is 0/1 sample 2 is 1/1 and sample 3 is 0/0
        // chr20:100005 -- alt is a SNP, sample 1 should be */0, sample 2 is */* and sample 3 is 0/1
        // chr20:100010 -- alt is a SNP, sample 1 should be 0/0, sample 2 is 0/0 and sample 3 is 0/1

        // TODO another test case would be a 1/2 sample with different length deletions and see it go from */* to */0 or something

        long location1 = SchemaUtils.encodeLocation("chr20", 100000);
        GenericRecord r11 = createVetRecord(location1, 1, "CTTGAGCATG", "C", "0/1", "30");
        GenericRecord r12 = createVetRecord(location1, 2, "CTTGAGCATG", "C", "1/1", "40");
        GenericRecord r23 = createVetRecord(location1 + 5, 3, "G", "T", "0/1", "50");
        GenericRecord r33 = createVetRecord(location1 + 10, 3, "C", "T", "0/1", "50");

        List<VariantContext> results = new ArrayList<>();

        List<GenericRecord> sortedRefRanges = new ArrayList<>();
        List<GenericRecord> sortedVet = Arrays.asList(r11, r12, r23, r33);

        engine.createVariantsFromSortedRanges(
                new TreeSet<>(sampleIdToName.keySet()),
                sortedVet,
                sortedRefRanges,
                new HashMap<>(),
                new HashMap<>(),
                new HashMap<>(),
                true,
                (vc) -> results.add(vc)
        );

        // we should have three "rows"
        Assert.assertEquals(results.size(), 3);

        VariantContext vc0 = results.get(0);
        Assert.assertEquals(vc0.getGenotype("sample_1").getGenotypeString(), "CTTGAGCATG/C");
        Assert.assertEquals(vc0.getGenotype("sample_1").getGQ(), 30);
        Assert.assertEquals(vc0.getGenotype("sample_2").getGenotypeString(), "C/C");
        Assert.assertEquals(vc0.getGenotype("sample_2").getGQ(), 40);
        Assert.assertEquals(vc0.getGenotype("sample_3").getGenotypeString(), "CTTGAGCATG/CTTGAGCATG");
        Assert.assertEquals(vc0.getGenotype("sample_3").getGQ(), 60);

        VariantContext vc1 = results.get(1);
        Assert.assertEquals(vc1.getGenotype("sample_1").getGenotypeString(), "G/*");
        Assert.assertEquals(vc1.getGenotype("sample_1").getGQ(), 30);
        Assert.assertEquals(vc1.getGenotype("sample_2").getGenotypeString(), "*/*");
        Assert.assertEquals(vc1.getGenotype("sample_2").getGQ(), 40);
        Assert.assertEquals(vc1.getGenotype("sample_3").getGenotypeString(), "G/T");
        Assert.assertEquals(vc1.getGenotype("sample_3").getGQ(), 50);

        VariantContext vc3 = results.get(2);
        Assert.assertEquals(vc3.getGenotype("sample_1").getGenotypeString(), "C/C");
        Assert.assertEquals(vc3.getGenotype("sample_1").getGQ(), 60);
        Assert.assertEquals(vc3.getGenotype("sample_2").getGenotypeString(), "C/C");
        Assert.assertEquals(vc3.getGenotype("sample_2").getGQ(), 60);
        Assert.assertEquals(vc3.getGenotype("sample_3").getGenotypeString(), "C/T");
        Assert.assertEquals(vc3.getGenotype("sample_3").getGQ(), 50);

        // for easier debugging visually
        // dumpToVcf(new File("out.vcf"), sampleIdToName.values(), results);
    }

    private void dumpToVcf(File outputVcf, Collection<String> sampleNames, List<VariantContext> variants) {
      VariantContextWriter w = GATKVariantContextUtils.createVCFWriter(outputVcf.toPath(), reference.getSequenceDictionary(), false);
      VCFHeader header = ExtractCohort.generateVcfHeader(new HashSet<String>(sampleNames), reference.getSequenceDictionary(), new HashSet<VCFHeaderLine>());
      w.writeHeader(header);
      for (VariantContext vc : variants) {
        w.add(vc);
      }
      w.close();
    }

    private GenericData.Record createVetRecord(long location, long sample_id, String ref, String alt, String gt, String gq) {
        GenericRecordBuilder builder = new GenericRecordBuilder(vetSchema);
        builder
                .set("sample_id", sample_id)
                .set("location", location)
                .set("ref", ref)
                .set("alt", alt)
                .set("call_GT", gt)
                .set("call_GQ", gq);
        return builder.build();
    }

}
