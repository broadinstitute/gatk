package org.broadinstitute.hellbender.utils.haplotype;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class FlowBasedHaplotypeUnitTest extends GATKBaseTest {

    @DataProvider(name="haplotypeGenerator")
    public Object[][] haplotypeTestDataGenerator(){
        String [] haplotypeSeqs = {"ATCGCAGGGAATTGTCCCCATGAAACTAAG",
                                "TGGGCTACCCCGTATATTTCGATTGCATTA",
                                "CCGCCTATTCGCTCTATCGCATCAAATCAA",
                                "GACGGCCTAGCTGCTCGTAGGCATCCTATA",
                                "ACTAACCGCATTTAACGCTCACGCATAAAG",
                                "TGAGTTTTCCACGACGTATTTCAGCTAAGA",
                                "AACTTTCACGTCACACAAGATTCCAGGTAC",
                                "GCACCGTCGTTTCGCCAATAAAATCGACTA",
                                "ATTTCCGCCCCTTAGACATTTGTATAACAT",
                                "GGACTTCGAAATTTTACAGATCATCGCTAC"};
        int[][] expectedFlow =  { {0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 3, 0, 2, 0, 0, 2, 0, 0, 1, 1, 0, 4, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 3, 1, 0, 1, 2, 0, 1},
                {1, 0, 0, 3, 0, 0, 1, 0, 1, 1, 4, 1, 1, 1, 0, 0, 1, 1, 0, 0, 3, 0, 1, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 1},
                {0, 0, 2, 1, 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 3, 0, 0, 1, 0, 1, 0, 0, 2},
                {0, 0, 0, 1, 0, 1, 1, 2, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2, 0, 1, 1, 0, 0, 1, 1},
                {0, 1, 1, 0, 1, 2, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 3, 2, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 3, 0, 1},
                {1, 0, 0, 1, 0, 1, 0, 1, 4, 0, 2, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 3, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 2, 0, 1, 0, 1},
                {0, 2, 1, 0, 3, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 0, 2, 0, 2, 0, 0, 1, 0, 2, 1, 1, 1},
                {0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 1, 3, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 4, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1},
                {0, 1, 0, 0, 3, 0, 2, 1, 0, 0, 4, 0, 2, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 3, 0, 0, 1, 1, 1, 0, 0, 1, 2, 1, 0, 0, 1, 0, 0, 1},
                {0, 0, 0, 2, 0, 1, 1, 0, 2, 0, 1, 1, 0, 3, 0, 0, 4, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1}};


        int [] trimLeft = {0,1,2,3,4,5,6,7,8,9};


        final List<Object[]> tests = new LinkedList<>();
        for (int i = 0; i < haplotypeSeqs.length; i++){
            tests.add( new Object[]{ haplotypeSeqs[i], expectedFlow[i], trimLeft[i], trimLeft[i]});

        }

        return tests.toArray(new Object[][]{});

    }

    @Test(dataProvider = "haplotypeGenerator")
    public void testFlowBasedHaplotype(String inputSeq, int[] expectedKey, int trimLeftIgnore, int trimRightIgnore){
        Haplotype hap = new Haplotype(inputSeq.getBytes());
        FlowBasedHaplotype fbh = new FlowBasedHaplotype(hap, "TACG");
        Assert.assertEquals(fbh.getKey(), expectedKey);


    }
    @Test(dataProvider = "haplotypeGenerator")
    public void testFindLeftClipping(String inputSeq, int[] expectedKey, int trimLeft, int trimRightIgnore) {
        Haplotype hap = new Haplotype(inputSeq.getBytes());
        Haplotype hapTrimmed = new Haplotype(inputSeq.substring(trimLeft).getBytes());
        FlowBasedHaplotype fbh = new FlowBasedHaplotype(hap, "TACG");
        int[] leftClip = fbh.findLeftClipping(trimLeft);
        int[] trimmedKey= Arrays.copyOfRange(fbh.getKey(), leftClip[0], fbh.getKeyLength());
        trimmedKey[0]-=leftClip[1];
        FlowBasedHaplotype fbhTrimmed = new FlowBasedHaplotype(hapTrimmed, "TACG");
        if (trimLeft>0) {
            expectedKey = Arrays.copyOfRange(fbhTrimmed.getKey(), findFirstNonzero(fbhTrimmed.getKey()), fbhTrimmed.getKeyLength());
        } else {
            expectedKey = fbhTrimmed.getKey();
        }
        Assert.assertEquals(trimmedKey, expectedKey);

    }


    @Test(dataProvider = "haplotypeGenerator")
    public void testFindRightClipping(String inputSeq, int[] expectedKey, int trimLeftIgnore, int trimRight) {
        Haplotype hap = new Haplotype(inputSeq.getBytes());
        Haplotype hapTrimmed = new Haplotype(inputSeq.substring(0, inputSeq.length() - trimRight).getBytes());
        FlowBasedHaplotype fbh = new FlowBasedHaplotype(hap, "TACG");
        int[] rightClip = fbh.findRightClipping(trimRight);
        int[] trimmedKey= Arrays.copyOfRange(fbh.getKey(), 0, fbh.getKeyLength()-rightClip[0]);
        trimmedKey[trimmedKey.length-1]-=rightClip[1];
        FlowBasedHaplotype fbhTrimmed = new FlowBasedHaplotype(hapTrimmed, "TACG");
        if (trimRight>0) {
            expectedKey = Arrays.copyOfRange(fbhTrimmed.getKey(), 0, findLastNonzero(fbhTrimmed.getKey())+1);
        } else {
            expectedKey = fbhTrimmed.getKey();
        }
        Assert.assertEquals(trimmedKey, expectedKey);

    }

    private int findFirstNonzero(int [] key){
        for (int i = 0; i < key.length; i++){
            if (key[i]>0){
                return i;
            }
        }
        return key.length;
    }

    private int findLastNonzero(int [] key){
        for (int i = key.length-1; i >= 0; i--){
            if (key[i]>0){
                return i;
            }
        }
        return -1;
    }
}