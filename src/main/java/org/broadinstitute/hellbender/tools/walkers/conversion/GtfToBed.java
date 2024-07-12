package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.util.Interval;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(
        summary = "Converts Gtf files to Bed file format", //TODO: make better summaries
        oneLineSummary = "Gtf to Bed",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class GtfToBed extends FeatureWalker <GencodeGtfFeature>{

    @Argument(fullName = "Gtf",
            shortName = "G", doc = "Gencode GTF file")
    public GATKPath inputFile;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output BED file")
    public String bedFilePath;

    private final Map<String, GtfInfo> idToInfo = new HashMap<>();

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(GencodeGtfFeature.class);
    }

    //Apply runs PER line
    @Override
    public void apply(GencodeGtfFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        int geneStart = 0;
        int geneEnd = 0;
        if(feature.getOptionalField("tag").contains("basic")) { //TODO: check if the key and value is actually correct for basic tag
            //Todo: I don't need to open the GTF file right?
            if (feature.getFeatureType() == GencodeGtfFeature.FeatureType.GENE) {
                geneStart = feature.getStart();
                geneEnd = feature.getEnd();
                Interval interval = new Interval(feature.getContig(), geneStart, geneEnd);
                GtfInfo gtfInfo = new GtfInfo(GtfInfo.Type.GENE, feature.getGeneName(), interval);
                idToInfo.put(feature.getGeneId(), gtfInfo);
            }
            if (feature.getFeatureType() == GencodeGtfFeature.FeatureType.TRANSCRIPT) {
                Interval interval = new Interval(feature.getContig(), feature.getStart(), feature.getEnd());
                GtfInfo gtfInfo = new GtfInfo(GtfInfo.Type.TRANSCRIPT, feature.getGeneName(), interval);
                idToInfo.put(feature.getTranscriptId(), gtfInfo);
                if (feature.getStart() < idToInfo.get(feature.getGeneId()).getStart()) { //Todo: Validate the assumption that the hashmap will always have the gene before the transcript is processed
                    geneStart = feature.getStart(); // doing this again instead store it in a variable
                    Interval geneInterval = new Interval(feature.getContig(), geneStart, geneEnd); //interval class is immutable, so I need to create a new interval
                    GtfInfo gtfGeneInfo = new GtfInfo(GtfInfo.Type.GENE, feature.getGeneName(), geneInterval);
                    idToInfo.put(feature.getGeneId(), gtfGeneInfo);
                }
                //Todo: Validate the assumption that the hashmap will always have the gene before the transcript is processed
                if (feature.getEnd() > idToInfo.get(feature.getGeneId()).getEnd()) {
                    geneEnd = feature.getEnd();
                    Interval geneInterval = new Interval(feature.getContig(), geneStart, geneEnd); //interval class is immutable, so I need to create a new interval
                    GtfInfo gtfGeneInfo = new GtfInfo(GtfInfo.Type.GENE, feature.getGeneName(), geneInterval);
                    idToInfo.put(feature.getGeneId(), gtfGeneInfo);
                }

                feature.getOptionalField("basic");

            }
        }
    }

    @Override
    public Object onTraversalSuccess(){
        /** TODO: make the user options and add sorting also move out of apply method*/
        //if user selects option 1:
        writeToBed(GtfInfo.Type.GENE);

        //if user selects option2:
        writeToBed(GtfInfo.Type.TRANSCRIPT);
        return null;
    }

    public void writeToBed(GtfInfo.Type type){
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(bedFilePath))) {
            for (Map.Entry<String, GtfInfo> entry : idToInfo.entrySet()) {
                if (entry.getValue().getType() == type) {
                    String line = entry.getValue().getInterval().getContig() + "\t" + //chr
                            entry.getValue().getInterval().getStart() + "\t" + //start
                            entry.getValue().getInterval().getEnd() + "\t" + //end
                            entry.getValue().getGeneName(); //gene name
                    if(type == GtfInfo.Type.TRANSCRIPT){
                        line = line + "," + entry.getKey();
                    }
                    writer.write(line);
                }
            }
        } catch (IOException e){
            e.printStackTrace(); //TODO: figure out better exception // Is there logging done in other programs. maybe write to a log file
        }
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFile;
    }
}

