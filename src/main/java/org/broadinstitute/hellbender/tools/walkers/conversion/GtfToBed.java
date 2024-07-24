package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Converts Gtf files to Bed file format", //TODO: make better summaries
        oneLineSummary = "Gtf to Bed",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class GtfToBed extends FeatureWalker<GencodeGtfFeature> {
    public static final String SORT_BY_TRANSCRIPT_LONG_NAME = "dont-mix-contigs";

    @Argument(fullName = "Gtf",
            shortName = "G", doc = "Gencode GTF file")
    public GATKPath inputFile;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output BED file")
    public GATKPath outputFile;

    @Argument(doc = "Make each row of BED file sorted by transcript (false is sorted by gene)", fullName = SORT_BY_TRANSCRIPT_LONG_NAME, optional = true)
    public boolean sortByTranscript = false;

    @Argument(doc = "path sequence dictionary", fullName = "Dictionary", shortName = "D")
    public GATKPath dictionaryPath;

    private final Map<String, GtfInfo> idToInfo = new HashMap<>();

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(GencodeGtfFeature.class);
    }

    @Override
    public void apply(GencodeGtfFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        List<GencodeGtfFeature> geneFeatures = feature.getAllFeatures();

        for (GencodeGtfFeature gtfFeature : geneFeatures) {
            processFeature(gtfFeature);
        }
    }

    private void processFeature(GencodeGtfFeature gtfFeature) {
        List<GencodeGtfFeature.OptionalField<?>> optionalFields = getOptionalFields(gtfFeature);

        if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.GENE) {
            processGeneFeature(gtfFeature);
        }

        for (GencodeGtfFeature.OptionalField<?> field : optionalFields) {
            if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.TRANSCRIPT && "basic".equals(field.getValue())) {
                processTranscriptFeature(gtfFeature);
            }
        }
    }

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
        int geneStart = gtfFeature.getStart();
        int geneEnd = gtfFeature.getEnd();
        Interval interval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
        GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
        idToInfo.put(gtfFeature.getGeneId(), gtfInfo);
    }

    private void processTranscriptFeature(GencodeGtfFeature gtfFeature) {
        Interval interval = new Interval(gtfFeature.getContig(), gtfFeature.getStart(), gtfFeature.getEnd());
        GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.TRANSCRIPT, gtfFeature.getGeneName());
        idToInfo.put(gtfFeature.getTranscriptId(), gtfInfo);
        updateGeneStartIfNeeded(gtfFeature);
        updateGeneEndIfNeeded(gtfFeature);
    }

    private void updateGeneStartIfNeeded(GencodeGtfFeature gtfFeature) {
        int geneStart = idToInfo.get(gtfFeature.getGeneId()).getStart();
        if (gtfFeature.getStart() < geneStart) {
            geneStart = gtfFeature.getStart();
            updateGeneInterval(gtfFeature, geneStart, idToInfo.get(gtfFeature.getGeneId()).getEnd());
        }
    }

    private void updateGeneEndIfNeeded(GencodeGtfFeature gtfFeature) {
        int geneEnd = idToInfo.get(gtfFeature.getGeneId()).getEnd();
        if (gtfFeature.getEnd() > geneEnd) {
            geneEnd = gtfFeature.getEnd();
            updateGeneInterval(gtfFeature, idToInfo.get(gtfFeature.getGeneId()).getStart(), geneEnd);
        }
    }

    private void updateGeneInterval(GencodeGtfFeature gtfFeature, int geneStart, int geneEnd) {
        Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
        GtfInfo gtfGeneInfo = new GtfInfo(geneInterval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
        idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
    }



    @Override
    public Object onTraversalSuccess() { //after going through every line of the gtf
        SamReader reader = SamReaderFactory.makeDefault().open(new File(String.valueOf(dictionaryPath))); //TODO: what is makeDefault() doing
        SAMSequenceDictionary sequenceDictionary = reader.getFileHeader().getSequenceDictionary(); //get the sequence dictionary

        List<Map.Entry<String, GtfInfo>> entries = new ArrayList<>(idToInfo.entrySet()); //convert the hashmap (populate in apply()) to a list
        entries.sort(new CompareGtfInfo(sequenceDictionary)); //sort the list by contig and then start

        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = new LinkedHashMap<>(); //create a new linked hashmap to preserve order
        for (Map.Entry<String, GtfInfo> entry : entries) {
            karyotypeIdToInfo.put(entry.getKey(), entry.getValue()); // put the sorted ordering of the entries the linked hashmap
        }
        if (sortByTranscript) { //if the user has selected the transcript option
            writeToBed(GtfInfo.Type.TRANSCRIPT, karyotypeIdToInfo);
        } else {
            writeToBed(GtfInfo.Type.GENE, karyotypeIdToInfo); //if user has selected gene option (default)
        }
        return null;
    }


    private void writeToBed(GtfInfo.Type type, Map<String, GtfInfo> sortedMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(String.valueOf(outputFile)))) { //write to bed file
            for (Map.Entry<String, GtfInfo> entry : sortedMap.entrySet()) { //for each entry in the sorted map
                if (entry.getValue().getType() == type) {
                    String line = entry.getValue().getInterval().getContig() + "\t" + //chr + tab
                            entry.getValue().getInterval().getStart() + "\t" + //start + tab
                            entry.getValue().getInterval().getEnd() + "\t" + //end + tab
                            entry.getValue().getGeneName(); //gene name

                    if (type == GtfInfo.Type.TRANSCRIPT) { //if the user selected transcript option
                        line = line + "," + entry.getKey(); //add a comma and the transcript id to the line
                    }
                    writer.write(line + System.lineSeparator()); //go to next line in file
                }
            }
        } catch (IOException e) {
            e.printStackTrace(); //TODO: figure out better exception
        }
    }


    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFile;
    }
}

