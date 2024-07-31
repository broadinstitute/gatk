package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Converts Gencode Gtf files to Bed file format with each row of bed file being either a gene or a transcript. ",
        oneLineSummary = "Gtf to Bed",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class GtfToBed extends FeatureWalker<GencodeGtfFeature> {
    public static final String SORT_BY_TRANSCRIPT_LONG_NAME = "dont-mix-contigs";

    @Argument(fullName = "Gencode Gtf input file",
            shortName = "G", doc = "Gencode GTF file")
    public GATKPath inputFile;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output BED file")
    public GATKPath outputFile;

    @Argument(doc = "Make each row of BED file sorted by transcript (false is sorted by gene)", shortName = "T", fullName = SORT_BY_TRANSCRIPT_LONG_NAME, optional = true)
    public boolean sortByTranscript = false;

    @Argument(doc = "path to sequence dictionary", fullName = "Dictionary", shortName = "D")
    public GATKPath dictionaryPath;

    private final Map<String, GtfInfo> idToInfo = new HashMap<>();

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(GencodeGtfFeature.class);
    }

    //runs per line of gtf file
    @Override
    public void apply(GencodeGtfFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        //list of all features of the gene
        List<GencodeGtfFeature> geneFeatures = feature.getAllFeatures();
        //process each gtf feature in the list of gene features
        for (GencodeGtfFeature gtfFeature : geneFeatures) {
            processFeature(gtfFeature);
        }
    }

    private void processFeature(GencodeGtfFeature gtfFeature) {
        //the basic tag is in optional fields
        List<GencodeGtfFeature.OptionalField<?>> optionalFields = getOptionalFields(gtfFeature);

        //if the gtf feature is a Gene
        if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.GENE) {
            processGeneFeature(gtfFeature);
        }

        //if the gtf feature is a transcript and has the basic tag
        for (GencodeGtfFeature.OptionalField<?> field : optionalFields) {
            if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.TRANSCRIPT && "basic".equals(field.getValue())) {
                processTranscriptFeature(gtfFeature);
            }
        }
    }

    //gets the tag out of the list of optional fields
    private List<GencodeGtfFeature.OptionalField<?>> getOptionalFields(GencodeGtfFeature gtfFeature) {
        List<GencodeGtfFeature.OptionalField<?>> optionalFields = null;
        try {
            optionalFields = gtfFeature.getOptionalField("tag");
        } catch (Exception e) {
            System.err.println("Could not retrieve optional fields: " + e.getMessage());
        }
        return optionalFields;
    }


    private void processGeneFeature(GencodeGtfFeature gtfFeature) {
        //gene start
        int geneStart = gtfFeature.getStart();
        //gene end
        int geneEnd = gtfFeature.getEnd();
        Interval interval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
        //put the interval, type as gene, and the name of gene
        GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
        //store in hashmap with key as geneId
        idToInfo.put(gtfFeature.getGeneId(), gtfInfo);
    }

    private void processTranscriptFeature(GencodeGtfFeature gtfFeature) {
        Interval interval = new Interval(gtfFeature.getContig(), gtfFeature.getStart(), gtfFeature.getEnd());
        //put the interval, type as transcript, and the name of the gene it's in
        GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.TRANSCRIPT, gtfFeature.getGeneName());
        //store in hashmap with key as transcriptId
        idToInfo.put(gtfFeature.getTranscriptId(), gtfInfo);
        //update start/end of corresponding gene if needed
        updateGeneStartIfNeeded(gtfFeature);
        updateGeneEndIfNeeded(gtfFeature);
    }

    private void updateGeneStartIfNeeded(GencodeGtfFeature gtfFeature) {
        //get the start value of the gene
        int geneStart = idToInfo.get(gtfFeature.getGeneId()).getStart();
        //if the transcript start is less than the gene start
        if (gtfFeature.getStart() < geneStart) {
            //set the gene start to be the transcript start
            geneStart = gtfFeature.getStart();
            updateGeneInterval(gtfFeature, geneStart, idToInfo.get(gtfFeature.getGeneId()).getEnd());
        }
    }

    private void updateGeneEndIfNeeded(GencodeGtfFeature gtfFeature) {
        //get the end value of the gene
        int geneEnd = idToInfo.get(gtfFeature.getGeneId()).getEnd();
        //if the transcript start is greater than the gene start
        if (gtfFeature.getEnd() > geneEnd) {
            //set the gene end to be the transcript end
            geneEnd = gtfFeature.getEnd();
            updateGeneInterval(gtfFeature, idToInfo.get(gtfFeature.getGeneId()).getStart(), geneEnd);
        }
    }

    //updates an interval of the gene if it needs to be changed
    private void updateGeneInterval(GencodeGtfFeature gtfFeature, int geneStart, int geneEnd) {
        Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
        GtfInfo gtfGeneInfo = new GtfInfo(geneInterval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
        idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
    }

   //immediately after it has gone through each line of gtf (apply method)
    @Override
    public Object onTraversalSuccess(){
        //get the user input dictionary
        SAMSequenceDictionary sequenceDictionary = getSequenceDictionary(String.valueOf(dictionaryPath));

        //create linked hash map to store sorted values of idToInfo
        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = getSortedIdToInfo(sequenceDictionary);

        //if user wants to sort by transcript only use transcripts else only use genes
        GtfInfo.Type selectedType = sortByTranscript ? GtfInfo.Type.TRANSCRIPT : GtfInfo.Type.GENE;
        writeToBed(selectedType, karyotypeIdToInfo);

        return null;
    }

    private SAMSequenceDictionary getSequenceDictionary(String dictionaryPath) {
        SamReader reader = SamReaderFactory.makeDefault().open(new File(dictionaryPath)); //TODO: figure out what makeDefault() does
        return reader.getFileHeader().getSequenceDictionary();
    }

    private LinkedHashMap<String, GtfInfo> getSortedIdToInfo(SAMSequenceDictionary sequenceDictionary) {
        //create a list that has the keys and values of idToInfo
        List<Map.Entry<String, GtfInfo>> entries = new ArrayList<>(idToInfo.entrySet());
        //sort the list using CompareGtfInfo
        entries.sort(new CompareGtfInfo(sequenceDictionary));

        //put each (sorted) entry in the list into a linked hashmap
        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = new LinkedHashMap<>();
        for (Map.Entry<String, GtfInfo> entry : entries) {
            karyotypeIdToInfo.put(entry.getKey(), entry.getValue());
        }
        return karyotypeIdToInfo;
    }

    //writes to bed file
    private void writeToBed(GtfInfo.Type type, Map<String, GtfInfo> sortedMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(String.valueOf(outputFile)))) {
            for (Map.Entry<String, GtfInfo> entry : sortedMap.entrySet()) {
                if (entry.getValue().getType() == type) {
                    String line = formatBedLine(entry, type);
                    writer.write(line + System.lineSeparator());
                }
            }
        } catch (IOException e) {
            throw new GATKException("Error writing to BED file", e);
        }
    }

    //formats each line of the bed file depending on whether user has selected gene or transcript
    private String formatBedLine(Map.Entry<String, GtfInfo> entry, GtfInfo.Type type) {
        GtfInfo info = entry.getValue();
        String line = info.getInterval().getContig() + "\t" +
                info.getInterval().getStart() + "\t" +
                info.getInterval().getEnd() + "\t" +
                info.getGeneName();

        if (type == GtfInfo.Type.TRANSCRIPT) {
            line += "," + entry.getKey();
        }
        return line;
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFile;
    }
}

