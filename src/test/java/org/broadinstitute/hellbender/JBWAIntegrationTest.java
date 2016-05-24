package org.broadinstitute.hellbender;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.File;

public class JBWAIntegrationTest extends BaseTest {

    private void skipJBWATestOnUnsupportedPlatforms() {
        if ( ! NativeUtils.runningOnLinux() && ! NativeUtils.runningOnMac() ) {
            throw new SkipException("jbwa not available on this platform");
        }
        if ( NativeUtils.runningOnPPCArchitecture() ) {
            throw new SkipException("jbwa not available for this architecture");
        }
    }

    @Test
    public void testJBWAIsLoadable(){
        skipJBWATestOnUnsupportedPlatforms();

        BWANativeLibrary.load();
    }


    @Test
    public void testJBWAAlignSingleRead() throws Exception {
        skipJBWATestOnUnsupportedPlatforms();

        BWANativeLibrary.load();

        BwaIndex index= new BwaIndex(new File(b37_reference_20_21));
        final BwaMem bwaMem = new BwaMem(index);
        //real read taken from src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam
        final String name="20FUKAAXX100202:3:46:9213:168594";
        final byte[] seqs= "GTTTTGTTTACTACAGCTTTGTAGTAAATTTTGAACTCTAAAGTGTTAGTTCTCTAACTTTGTTTGTTTTTCAAGAGTGTTTTGACTCTTCTTACTGCATC".getBytes(); ;
        final byte[] quals= "DGFDGFDHFFFFGFEFHEGFFFGGHEHFHGHHGEGGGGGFGFHGHGHEHGGGFGAEFGDACAHHDHGCGFGGGFGDHGHFFHDDCGGDGEE".getBytes();
        final ShortRead read = new ShortRead(name, seqs, quals);
        final AlnRgn[] align = bwaMem.align(read);
        Assert.assertEquals(align.length, 1);
        Assert.assertEquals(align[0].getChrom(), "20");
        Assert.assertEquals(align[0].getCigar(), "101M");
        Assert.assertEquals(align[0].getMQual(), 60);
        Assert.assertEquals(align[0].getPos(), 9999997-1); // note difference from the bam file (9999997 in bam, 9999996 here)
        Assert.assertEquals(align[0].getNm(), 0);
        bwaMem.dispose();
    }
}
