package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.GATKBaseTest;

public class ConversionTestUtils {

    public static String getTestDataDir(){
        return GATKBaseTest.toolsTestDir + "/walkers/gtfToBed/";
    }

    public static String getDecoySampleGtf(){
        return getTestDataDir()+"/decoySample.gtf";
    }

    public static String getMapk1Gtf(){
        return GATKBaseTest.largeFileTestDir + "/gtfToBed/mapk1.gtf";
    }

    public static String getManyTranscriptsGtf(){
        return getTestDataDir() + "/manyTranscripts.gtf";
    }

    public static String getReferenceDict(){
        return getTestDataDir()+"/hg38GencodeReference.dict";
    }

    public static String getChr14GeneGtf(){
        return  GATKBaseTest.largeFileTestDir + "/gtfToBed/ENSG00000215398.12.gtf";
    }

    public static String getChr14GeneBed(){
        return getTestDataDir() + "/ENSG00000215398.12.bed";
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

    public static String getManyTranscriptsBed(){
        return getTestDataDir() + "/manyTranscripts.bed";
    }

    public static String getNotBasicBed(){
        return getTestDataDir() + "/notBasic.bed";
    }

    public static String getMouseGtf(){
        return getTestDataDir() + "/mouseSample.gtf";
    }

    public static String getMouseDict(){
        return getTestDataDir() + "/mouseReference.dict";
    }

    public static String getMouseBed(){
        return getTestDataDir() + "/mouse.bed";
    }

    public static String getChr14Fasta() {
        return GATKBaseTest.largeFileTestDir + "/GRCm38_primary_assembly_genome/chr14.GRCm38.primary_assembly.genome.fa.gz";
    }
}
