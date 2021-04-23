package org.broadinstitute.hellbender.tools.htsgetreader;

import java.io.IOException;
import java.net.URI;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class HtsgetResponseUnitTest extends GATKBaseTest {
    @Test
    public void testDeserialization() throws JsonParseException, JsonMappingException, IOException {
        final ObjectMapper mapper = new ObjectMapper();
        mapper.enable(DeserializationFeature.UNWRAP_ROOT_VALUE);
        final String respJson = "{\n   \"htsget\" : {\n      \"format\" : \"BAM\",\n      \"urls\" : [\n         {\n            \"url\" : \"data:application/vnd.ga4gh.bam;base64,QkFNAQ==\",\n            \"class\" : \"header\"\n         },\n         {\n            \"url\" : \"https://htsget.blocksrv.example/sample1234/header\",\n            \"class\" : \"header\"\n         },\n         {\n            \"url\" : \"https://htsget.blocksrv.example/sample1234/run1.bam\",\n            \"headers\" : {\n               \"Authorization\" : \"Bearer xxxx\",\n               \"Range\" : \"bytes=65536-1003750\"\n             },\n            \"class\" : \"body\"\n         },\n         {\n            \"url\" : \"https://htsget.blocksrv.example/sample1234/run1.bam\",\n            \"headers\" : {\n               \"Authorization\" : \"Bearer xxxx\",\n               \"Range\" : \"bytes=2744831-9375732\"\n            },\n            \"class\" : \"body\"\n         }\n      ]\n   }\n}";
        final HtsgetResponse resp = mapper.readValue(respJson, HtsgetResponse.class);

        Assert.assertEquals(resp.getFormat(), HtsgetFormat.BAM);

        Assert.assertEquals(resp.getBlocks().get(0).getUri(), URI.create("data:application/vnd.ga4gh.bam;base64,QkFNAQ=="));
        Assert.assertEquals(resp.getBlocks().get(0).getDataClass(), HtsgetClass.header);

        Assert.assertEquals(resp.getBlocks().get(2).getHeaders().get("Authorization"), "Bearer xxxx");
        Assert.assertEquals(resp.getBlocks().get(2).getHeaders().get("Range"), "bytes=65536-1003750");
    }
}
