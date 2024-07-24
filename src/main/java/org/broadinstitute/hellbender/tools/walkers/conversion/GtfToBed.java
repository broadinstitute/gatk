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
    public Object onTraversalSuccess(){
        SAMSequenceDictionary sequenceDictionary = getSequenceDictionary(String.valueOf(dictionaryPath));
        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = getSortedIdToInfo(sequenceDictionary);

        GtfInfo.Type selectedType = sortByTranscript ? GtfInfo.Type.TRANSCRIPT : GtfInfo.Type.GENE;
        writeToBed(selectedType, karyotypeIdToInfo);

        return null;
    }

    private SAMSequenceDictionary getSequenceDictionary(String dictionaryPath) {
        SamReader reader = SamReaderFactory.makeDefault().open(new File(dictionaryPath)); //TODO: figure out what makeDefault() does
        return reader.getFileHeader().getSequenceDictionary();
    }

    private LinkedHashMap<String, GtfInfo> getSortedIdToInfo(SAMSequenceDictionary sequenceDictionary) {
        List<Map.Entry<String, GtfInfo>> entries = new ArrayList<>(idToInfo.entrySet());
        entries.sort(new CompareGtfInfo(sequenceDictionary));

        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = new LinkedHashMap<>();
        for (Map.Entry<String, GtfInfo> entry : entries) {
            karyotypeIdToInfo.put(entry.getKey(), entry.getValue());
        }
        return karyotypeIdToInfo;
    }

    private void writeToBed(GtfInfo.Type type, Map<String, GtfInfo> sortedMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(String.valueOf(outputFile)))) {
            for (Map.Entry<String, GtfInfo> entry : sortedMap.entrySet()) {
                if (entry.getValue().getType() == type) {
                    String line = formatBedLine(entry, type);
                    writer.write(line + System.lineSeparator());
                }
            }
        } catch (IOException e) {
            System.err.println("Error writing to BED file: " + e.getMessage());
        }
    }

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

