package org.broadinstitute.hellbender.tools.walkers.conversion;

public class ConversionTestUtils {

    public static String getTestDataDir(){
        return "src/test/resources/org/broadinstitute/hellbender/tools/walkers/conversion/";
    }

    public static String getDecoySampleGtf(){
        return getTestDataDir()+"/decoySample.gtf";
    }

    public static String getMapk1Gtf(){
        return getTestDataDir()+"/mapk1.gtf";
    }

    public static String getReferenceDict(){
        return getTestDataDir()+"/reference.dict";
    }

    public static String getDecoySamplesGeneBed(){
        return getTestDataDir() + "/decoySampleGene.bed";
    }

    public static String getDecoySamplesTranscriptBed(){
        return getTestDataDir() + "/decoySampleTranscript.bed";
    }

    public static String getMapk1GeneBed(){
        return getTestDataDir() + "/mapk1Gene.bed";
    }

    public static String getMapk1TranscriptBed(){
        return getTestDataDir() + "/mapk1Transcript.bed";
    }
}
