package org.broadinstitute.hellbender.tools.walkers.conversion;

public class ConversionTestUtils {

    public static String getTestDataDir(){
        return "src/test/resources/org/broadinstitute/hellbender/tools/walkers/conversion/testFiles";
    }

    public static String getDecoysGtf(){
        return getTestDataDir()+"/decoys.gtf";
    }

    public static String getGencodeGtf(){
        return getTestDataDir()+"/gencode.v38.chr_patch_hapl_scaff.annotation.gtf";
    }

    public static String getMapk1Gtf(){
        return getTestDataDir()+"/mapk1.gtf";
    }

    public static String getReferenceDict(){
        return getTestDataDir()+"/reference.dict";
    }
}
