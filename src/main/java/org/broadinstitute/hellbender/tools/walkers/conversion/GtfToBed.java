package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.util.Interval;
import htsjdk.tribble.Feature;
import org.apache.parquet.filter2.predicate.Operators;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfCodec;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeatureBaseData;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Converts Gtf files to Bed file format", //TODO: make better summaries
        oneLineSummary = "Gtf to Bed",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class GtfToBed extends FeatureWalker<GencodeGtfFeature> {

    @Argument(fullName = "Gtf",
            shortName = "G", doc = "Gencode GTF file")
    public GATKPath inputFile;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output BED file")
    public GATKPath outputFile;

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
        List<GencodeGtfFeature> geneFeature = feature.getAllFeatures();
        for (GencodeGtfFeature gtfFeature : geneFeature) {
            List<GencodeGtfFeature.OptionalField<?>> optionalFields = gtfFeature.getOptionalField("tag");
                    if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.GENE) {
                        System.out.println("Type is: " + gtfFeature.getFeatureType() + " gene id = " + gtfFeature.getGeneId());
                        geneStart = gtfFeature.getStart();
                        geneEnd = gtfFeature.getEnd();
                        Interval interval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
                        GtfInfo gtfInfo = new GtfInfo(GtfInfo.Type.GENE, gtfFeature.getGeneName(), interval);
                        idToInfo.put(gtfFeature.getGeneId(), gtfInfo);
                    }
           for (GencodeGtfFeature.OptionalField<?> field : optionalFields) {

                    if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.TRANSCRIPT && field.getValue().equals("basic")) {
                        System.out.println("type is: " + gtfFeature.getFeatureType() + " transcript id = " + gtfFeature.getTranscriptId());
                        Interval interval = new Interval(gtfFeature.getContig(), gtfFeature.getStart(), gtfFeature.getEnd());
                        GtfInfo gtfInfo = new GtfInfo(GtfInfo.Type.TRANSCRIPT, gtfFeature.getGeneName(), interval);
                        idToInfo.put(gtfFeature.getTranscriptId(), gtfInfo);

                        if (gtfFeature.getStart() < idToInfo.get(gtfFeature.getGeneId()).getStart()) { //Todo: Validate the assumption that the hashmap will always have the gene before the transcript is processed
                            geneStart = gtfFeature.getStart(); // doing this again instead store it in a variable
                            Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd); //interval class is immutable, so I need to create a new interval
                            GtfInfo gtfGeneInfo = new GtfInfo(GtfInfo.Type.GENE, gtfFeature.getGeneName(), geneInterval);
                            idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
                        }

                        if (gtfFeature.getEnd() > idToInfo.get(gtfFeature.getGeneId()).getEnd()) {
                            geneEnd = gtfFeature.getEnd();
                            Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd); //interval class is immutable, so I need to create a new interval
                            GtfInfo gtfGeneInfo = new GtfInfo(GtfInfo.Type.GENE, gtfFeature.getGeneName(), geneInterval);
                            idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
                        }
                    }
           }
        }
    }


    @Override
    public Object onTraversalSuccess() {
        List<Map.Entry<String, GtfInfo>> list = new ArrayList<>(idToInfo.entrySet());
        list.sort(Map.Entry.comparingByValue(Comparator.comparing(GtfInfo::getInterval)));
        Map<String, GtfInfo> sortedIdtoInfo = new LinkedHashMap<>();
        for (Map.Entry<String, GtfInfo> entry : list) {
            sortedIdtoInfo.put(entry.getKey(), entry.getValue());
        }

        /** TODO: make the user options and add sorting also move out of apply method*/
        //if user selects option 1:
        //writeToBed(GtfInfo.Type.GENE, sortedIdtoInfo);

        //if user selects option2:
        //writeToBed(GtfInfo.Type.TRANSCRIPT, sortedIdtoInfo);
        return null;
    }

    public void writeToBed(GtfInfo.Type type, Map<String, GtfInfo> sortedMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(String.valueOf(outputFile)))) { //TODO: figure out what should be the output file path
            for (Map.Entry<String, GtfInfo> entry : sortedMap.entrySet()) {
                if (entry.getValue().getType() == type) {
                    String line = "contig = " + entry.getValue().getInterval().getContig() + "\t" + //chr
                            " start = " + entry.getValue().getInterval().getStart() + "\t" + //start
                            " end = " + entry.getValue().getInterval().getEnd() + "\t" + //end
                            " gene name = " + entry.getValue().getGeneName(); //gene name
                    if (type == GtfInfo.Type.TRANSCRIPT) {
                        line = line + "," + "transcript id = " + entry.getKey();
                    }
                    writer.write(line + System.lineSeparator());
                }
            }
        } catch (IOException e) {
            e.printStackTrace(); //TODO: figure out better exception // Is there logging done in other programs. maybe write to a log file
        }
        System.out.println("writing to bed");
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFile;
    }
}

