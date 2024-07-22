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

    //File dictFile = new File("/Users/shahsana/TestGatk/reference.dict");


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
            List<GencodeGtfFeature.OptionalField<?>> optionalFields = null;
            try {
                optionalFields = gtfFeature.getOptionalField("tag");
            } catch (Exception e) {
                System.err.println("Could not retrieve optional fields" + e.getMessage());
                //continue; //TODO: do i need continue? or return? or nothing?
            }
            try {
                if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.GENE) {
                    System.out.println("Type is: " + gtfFeature.getFeatureType() + " gene id = " + gtfFeature.getGeneId());
                    geneStart = gtfFeature.getStart();
                    geneEnd = gtfFeature.getEnd();
                    Interval interval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
                    GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
                    idToInfo.put(gtfFeature.getGeneId(), gtfInfo);
                }
                for (GencodeGtfFeature.OptionalField<?> field : optionalFields) {
                    if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.TRANSCRIPT && field.getValue().equals("basic")) {
                        System.out.println("type is: " + gtfFeature.getFeatureType() + " transcript id = " + gtfFeature.getTranscriptId());
                        Interval interval = new Interval(gtfFeature.getContig(), gtfFeature.getStart(), gtfFeature.getEnd());
                        GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.TRANSCRIPT, gtfFeature.getGeneName());
                        idToInfo.put(gtfFeature.getTranscriptId(), gtfInfo);

                        if (gtfFeature.getStart() < idToInfo.get(gtfFeature.getGeneId()).getStart()) { //Todo: Validate the assumption that the hashmap will always have the gene before the transcript is processed
                            geneStart = gtfFeature.getStart();
                            Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
                            GtfInfo gtfGeneInfo = new GtfInfo(geneInterval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
                            idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
                        }

                        if (gtfFeature.getEnd() > idToInfo.get(gtfFeature.getGeneId()).getEnd()) {
                            geneEnd = gtfFeature.getEnd();
                            Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
                            GtfInfo gtfGeneInfo = new GtfInfo(geneInterval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
                            idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
                        }
                    }
                }

            } catch (Exception e) {
                System.err.println("Error processing feature: " + e.getMessage());
            }
        }
    }


    @Override
    public Object onTraversalSuccess() {

        SamReader reader = SamReaderFactory.makeDefault().open(new File(String.valueOf(dictionaryPath)));
        SAMSequenceDictionary sequenceDictionary = reader.getFileHeader().getSequenceDictionary();

        List<Map.Entry<String, GtfInfo>> entries = new ArrayList<>(idToInfo.entrySet());
        entries.sort(new CompareGtfInfo(sequenceDictionary));

        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = new LinkedHashMap<>();
        for (Map.Entry<String, GtfInfo> entry : entries) {
            karyotypeIdToInfo.put(entry.getKey(), entry.getValue());
        }
        if (sortByTranscript) {
            writeToBed(GtfInfo.Type.TRANSCRIPT, karyotypeIdToInfo);
        } else {
            writeToBed(GtfInfo.Type.GENE, karyotypeIdToInfo);
        }
        return null;
    }


    public void writeToBed(GtfInfo.Type type, Map<String, GtfInfo> sortedMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(String.valueOf(outputFile)))) { //TODO: figure out what should be the output file path
            for (Map.Entry<String, GtfInfo> entry : sortedMap.entrySet()) {
                if (entry.getValue().getType() == type) {
                    String line = entry.getValue().getInterval().getContig() + "\t" + //chr
                            entry.getValue().getInterval().getStart() + "\t" + //start
                            entry.getValue().getInterval().getEnd() + "\t" + //end
                            entry.getValue().getGeneName(); //gene name

                    if (type == GtfInfo.Type.TRANSCRIPT) {
                        line = line + "," + entry.getKey();
                    }
                    writer.write(line + System.lineSeparator());
                }
            }
        } catch (IOException e) {
            e.printStackTrace(); //TODO: figure out better exception
        }
        System.out.println("writing to bed");
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFile;
    }
}

