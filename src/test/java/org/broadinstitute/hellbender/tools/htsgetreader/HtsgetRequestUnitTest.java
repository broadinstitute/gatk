package org.broadinstitute.hellbender.tools.htsgetreader;

import java.net.URI;
import java.net.URISyntaxException;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class HtsgetRequestUnitTest extends GATKBaseTest {

    private static final String endpoint = "https://example.com";

    @Test
    public void testOnlyId() throws URISyntaxException {
        final HtsgetRequestBuilder req = new HtsgetRequestBuilder(new URI(endpoint), "1");
        Assert.assertEquals(req.toURI().toString(), "https://example.com/1");
    }

    @Test
    public void testBasicFields() throws URISyntaxException {
        final HtsgetRequestBuilder req = new HtsgetRequestBuilder(new URI(endpoint), "1")
            .withFormat(HtsgetFormat.BAM)
            .withDataClass(HtsgetClass.body)
            .withInterval(new SimpleInterval("chr1:1-16"));
        Assert.assertEquals(req.toURI().toString(), "https://example.com/1?format=BAM&class=body&referenceName=chr1&start=0&end=16");
    }

    @Test
    public void testCompositeFields() throws URISyntaxException {
        final HtsgetRequestBuilder req = new HtsgetRequestBuilder(new URI(endpoint), "1")
            .withField(HtsgetRequestField.QNAME)
            .withField(HtsgetRequestField.FLAG)
            .withTag("tag1")
            .withTag("tag2")
            .withNotag("tag3")
            .withNotag("tag4");
        final String query = req.toURI().toString();
        Assert.assertTrue(query.contains("fields=QNAME%2CFLAG") || query.contains("fields=FLAG%2CQNAME"));
        Assert.assertTrue(query.contains("tags=tag1%2Ctag2") || query.contains("tags=tag2%2Ctag1"));
        Assert.assertTrue(query.contains("notags=tag3%2Ctag4") || query.contains("notags=tag4%2Ctag3"));
    }

    @DataProvider(name = "invalidParams")
    public Object[][] invalidParams() throws URISyntaxException {
        return new Object[][]{
            // class=header while interval also specified
            {new HtsgetRequestBuilder(new URI(endpoint), "1")
                .withDataClass(HtsgetClass.header)
                .withInterval(new SimpleInterval("chr1:1-16"))},
            // class=header while field also specified
            {new HtsgetRequestBuilder(new URI(endpoint), "1")
                .withDataClass(HtsgetClass.header)
                .withField(HtsgetRequestField.QNAME)},
            // class=header while tag also specified
            {new HtsgetRequestBuilder(new URI(endpoint), "1")
                .withDataClass(HtsgetClass.header)
                .withTag("NH")},
            // class=header while notag also specified
            {new HtsgetRequestBuilder(new URI(endpoint), "1")
                .withDataClass(HtsgetClass.header)
                .withNotag("NH")},
            // tags and notags overlap
            {new HtsgetRequestBuilder(new URI(endpoint), "1")
                .withDataClass(HtsgetClass.header)
                .withTag("NH")
                .withNotag("NH")},
            // .bam file requested while format is variant
            {new HtsgetRequestBuilder(new URI(endpoint), "example.bam")
                .withFormat(HtsgetFormat.VCF)
            },
            // .vcf file requested while format is read
            {new HtsgetRequestBuilder(new URI(endpoint), "example.vcf")
                .withFormat(HtsgetFormat.BAM)
            }
        };
    }

    // Expect a validation failure for invalid combinations of query parameters
    @Test(dataProvider = "invalidParams", expectedExceptions = UserException.class)
    public void testValidationFailure(final HtsgetRequestBuilder query) {
        query.toURI();
    }
}
