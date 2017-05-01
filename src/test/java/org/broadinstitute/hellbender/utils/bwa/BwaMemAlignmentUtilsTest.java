package org.broadinstitute.hellbender.utils.bwa;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class BwaMemAlignmentUtilsTest {
    @Test
    public void testSATags() {
        final List<String> refNames = new ArrayList<>(5);
        refNames.add("chr1");
        refNames.add("chr4");
        refNames.add("chr12");
        refNames.add("chr20");
        refNames.add("chrUn_KI270442v1");
        final List<BwaMemAlignment> alignments = new ArrayList<>(7);
        alignments.add(new BwaMemAlignment(0,4,71180,71338,404,568,40,20,78,48,"404S94M3I10M5I54M66S","11G22T14T15T0T2A13G7A3T1T8C21C29",null,-1,0,0));
        alignments.add(new BwaMemAlignment(2064,1,49146596,49146678,517,599,0,2,72,71,"517H82M37H","60C10T10",null,-1,0,0));
        alignments.add(new BwaMemAlignment(2048,4,70496,70594,98,192,18,8,58,48,"98H98M440H","10A10G21T14G15T1G3T7A9",null,-1,0,0));
        alignments.add(new BwaMemAlignment(2048,0,790150,790188,283,321,8,0,38,35,"283H38M315H","38",null,-1,0,0));
        alignments.add(new BwaMemAlignment(2048,3,31231128,31231162,345,379,0,0,34,34,"345H34M257H","34",null,-1,0,0));
        alignments.add(new BwaMemAlignment(2048,2,97988540,97988587,250,297,0,3,32,30,"250H47M339H","29G2A4G9",null,-1,0,0));
        alignments.add(new BwaMemAlignment(2064,1,49141579,49141624,583,628,0,3,30,30,"583H45M8H","19C0A4A19",null,-1,0,0));
        final String[] expectedTags = new String[alignments.size()];
        expectedTags[0] = "chr4,49146597,-,517S82M37S,0,2;chrUn_KI270442v1,70497,+,98S98M440S,18,8;chr1,790151,+,283S38M315S,8,0;chr20,31231129,+,345S34M257S,0,0;chr12,97988541,+,250S47M339S,0,3;chr4,49141580,-,583S45M8S,0,3;";
        expectedTags[1] = "chrUn_KI270442v1,71181,+,404S94M3I10M5I54M66S,40,20;chrUn_KI270442v1,70497,+,98S98M440S,18,8;chr1,790151,+,283S38M315S,8,0;chr20,31231129,+,345S34M257S,0,0;chr12,97988541,+,250S47M339S,0,3;chr4,49141580,-,583S45M8S,0,3;";
        expectedTags[2] = "chrUn_KI270442v1,71181,+,404S94M3I10M5I54M66S,40,20;chr4,49146597,-,517S82M37S,0,2;chr1,790151,+,283S38M315S,8,0;chr20,31231129,+,345S34M257S,0,0;chr12,97988541,+,250S47M339S,0,3;chr4,49141580,-,583S45M8S,0,3;";
        expectedTags[3] = "chrUn_KI270442v1,71181,+,404S94M3I10M5I54M66S,40,20;chr4,49146597,-,517S82M37S,0,2;chrUn_KI270442v1,70497,+,98S98M440S,18,8;chr20,31231129,+,345S34M257S,0,0;chr12,97988541,+,250S47M339S,0,3;chr4,49141580,-,583S45M8S,0,3;";
        expectedTags[4] = "chrUn_KI270442v1,71181,+,404S94M3I10M5I54M66S,40,20;chr4,49146597,-,517S82M37S,0,2;chrUn_KI270442v1,70497,+,98S98M440S,18,8;chr1,790151,+,283S38M315S,8,0;chr12,97988541,+,250S47M339S,0,3;chr4,49141580,-,583S45M8S,0,3;";
        expectedTags[5] = "chrUn_KI270442v1,71181,+,404S94M3I10M5I54M66S,40,20;chr4,49146597,-,517S82M37S,0,2;chrUn_KI270442v1,70497,+,98S98M440S,18,8;chr1,790151,+,283S38M315S,8,0;chr20,31231129,+,345S34M257S,0,0;chr4,49141580,-,583S45M8S,0,3;";
        expectedTags[6] = "chrUn_KI270442v1,71181,+,404S94M3I10M5I54M66S,40,20;chr4,49146597,-,517S82M37S,0,2;chrUn_KI270442v1,70497,+,98S98M440S,18,8;chr1,790151,+,283S38M315S,8,0;chr20,31231129,+,345S34M257S,0,0;chr12,97988541,+,250S47M339S,0,3;";
        final Map<BwaMemAlignment,String> actualTagMap = BwaMemAlignmentUtils.createSATags(alignments,refNames);
        for ( int idx = 0; idx != expectedTags.length; ++idx ) {
            Assert.assertEquals(actualTagMap.get(alignments.get(idx)), expectedTags[idx]);
        }
    }
}
