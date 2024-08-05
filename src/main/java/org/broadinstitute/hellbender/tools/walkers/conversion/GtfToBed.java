package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.*; //TODO: make these individual imports
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

/**
 * <p>Convert GTF files to BED format with options for gene and transcript level processing.
 * This tool allows for the extraction of gene and transcript information from GTF files and
 * outputs the data in BED format.</p>
 *
 *
 * <p>The conversion process includes sorting entries
 * by karyotype, providing flexibility in the selection of either gene or transcript level
 * data, and an option to only use basic tags. It ensures that the BED output is sorted and formatted correctly for subsequent use.</p>
 *
 * <h3>Usage examples</h3>
 * <p>Example commands to run GtfToBed for typical scenarios:</p>
 *
 * <h4>(i) Convert GTF to BED with gene level data</h4>
 * <p>This mode extracts and converts gene data from the input GTF file to BED format:</p>
 *
 * <pre>
 *     java -jar GtfToBed.jar \
 *     -G input.gtf \
 *     -SD dictionary.dict \
 *     -T False \
 *     -B False \
 *     -O output.bed \
 * </pre>
 *
 * <h4>(ii) Convert GTF to BED with transcript level data</h4>
 * <p>This mode extracts and converts transcript data from the input GTF file to BED format:</p>
 *
 * <pre>
 *     java -jar GtfToBed.jar \
 *     -G input.gtf \
 *     -SD dictionary.dict \
 *     -T True \
 *     -B False \
 *     -O output.bed \
 * </pre>
 *
 * <h4>(iii) Convert GTF to BED with transcript level data filtering for only those with the basic tag</h4>
 *  * <p>This mode extracts and converts basic transcript data from the input GTF file to BED format:</p>
 *  *
 *  * <pre>
 *     java -jar GtfToBed.jar \
 *     -G input.gtf \
 *     -SD dictionary.dict \
 *     -T True \
 *     -B True \
 *     -O output.bed \
 *  * </pre>
 */

@CommandLineProgramProperties(
        summary = "Converts GTF files to Bed file format with each row of bed file being either a gene or a transcript.",
        oneLineSummary = "GTF to BED",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)

//TODO: make args all caps
public class GtfToBed extends FeatureWalker<GencodeGtfFeature> {
    public static final String SORT_BY_TRANSCRIPT_LONG_NAME = "SORT_BY_TRANSCRIPT";
    public static final String SEQUENCE_DICTIONARY_LONG_NAME = "SEQUENCE_DICTIONARY";
    public static final String SORT_BY_BASIC_LONG_NAME = "SORT_BY_BASIC";



    @Argument(fullName = "GTF",
            shortName = "G", doc = "Path to Gencode GTF file")
    public GATKPath inputFile;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output BED file")
    public GATKPath outputFile;

    @Argument(fullName = SORT_BY_TRANSCRIPT_LONG_NAME, shortName = "T", doc = "Make each row of BED file sorted by transcript", optional = true)
    public boolean sortByTranscript = false;

    @Argument( fullName = SEQUENCE_DICTIONARY_LONG_NAME, shortName = "SD", doc = "Path to sequence dictionary")
    public GATKPath sequenceDictionary;

    @Argument( fullName = SORT_BY_BASIC_LONG_NAME, shortName = "B", doc = "Only use basic transcripts")
    public boolean sortByBasic = false;

    private final Map<String, GtfInfo> idToInfo = new HashMap<>();

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(GencodeGtfFeature.class);
    }

    // runs per line of gtf file
    @Override
    public void apply(GencodeGtfFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // list of all features of the gene
        List<GencodeGtfFeature> geneFeatures = feature.getAllFeatures();
        // process each gtf feature in the list of gene features
        for (GencodeGtfFeature gtfFeature : geneFeatures) {
            processFeature(gtfFeature);
        }
    }

    // gets both genes and transcripts from GTF and puts it in a hashmap
    private void processFeature(GencodeGtfFeature gtfFeature) {
        // the basic tag is in optional fields
        List<GencodeGtfFeature.OptionalField<?>> optionalFields = getOptionalFields(gtfFeature);

        // if the gtf feature is a Gene
        if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.GENE) {
            processGeneFeature(gtfFeature);
        }

        //  check if the gtf feature is a transcript. If user only wants basic transcripts check that it has the basic tag
        if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.TRANSCRIPT) {
            if (sortByBasic) {
                for (GencodeGtfFeature.OptionalField<?> field : optionalFields) {
                    if ("basic".equals(field.getValue())) {
                        processTranscriptFeature(gtfFeature);
                    }
                }
            } else {
                processTranscriptFeature(gtfFeature);
            }
        }
    }

    // gets the tag out of the list of optional fields
    private List<GencodeGtfFeature.OptionalField<?>> getOptionalFields(GencodeGtfFeature gtfFeature) {
        List<GencodeGtfFeature.OptionalField<?>> optionalFields = null;
        try {
            optionalFields = gtfFeature.getOptionalField("tag");
        } catch (Exception e) {
            System.err.println("Could not retrieve optional fields: " + e.getMessage());
        }
        return optionalFields;
    }

    // stores the gene ID and Interval info in hashmap
    private void processGeneFeature(GencodeGtfFeature gtfFeature) {
        final int geneStart = gtfFeature.getStart();
        final int geneEnd = gtfFeature.getEnd();
        final Interval interval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);

        // put the interval, type as gene, and the name of gene
        final GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.GENE, gtfFeature.getGeneName());

        // store in hashmap with key as geneId
        idToInfo.put(gtfFeature.getGeneId(), gtfInfo);
    }

    // stores the transcript ID and Interval info in hashmap
    private void processTranscriptFeature(GencodeGtfFeature gtfFeature) {
        //get interval and put the interval, type as transcript, and the name of the gene it's in
        final Interval interval = new Interval(gtfFeature.getContig(), gtfFeature.getStart(), gtfFeature.getEnd());
        final GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.TRANSCRIPT, gtfFeature.getGeneName());

        //store in hashmap with key as transcriptId
        idToInfo.put(gtfFeature.getTranscriptId(), gtfInfo);

        //update start/end of corresponding gene if needed
        updateGeneStart(gtfFeature);
        updateGeneEnd(gtfFeature);
    }

    // update the gene interval start position based on the transcript
    private void updateGeneStart(GencodeGtfFeature gtfFeature) {
        // get the start value of the gene
        int geneStart = idToInfo.get(gtfFeature.getGeneId()).getStart();

        // if the transcript start is less than the gene start
        if (gtfFeature.getStart() < geneStart) {
            // set the gene start to be the transcript start
            geneStart = gtfFeature.getStart();
            updateGeneInterval(gtfFeature, geneStart, idToInfo.get(gtfFeature.getGeneId()).getEnd());
        }
    }

    // update the gene interval end position based on the transcript
    private void updateGeneEnd(GencodeGtfFeature gtfFeature) {
        // get the end value of the gene
        int geneEnd = idToInfo.get(gtfFeature.getGeneId()).getEnd();

        // if the transcript start is greater than the gene start
        if (gtfFeature.getEnd() > geneEnd) {
            // set the gene end to be the transcript end
            geneEnd = gtfFeature.getEnd();
            updateGeneInterval(gtfFeature, idToInfo.get(gtfFeature.getGeneId()).getStart(), geneEnd);
        }
    }

    // updates an interval of the gene if it needs to be changed
    private void updateGeneInterval(GencodeGtfFeature gtfFeature, int geneStart, int geneEnd) {
        Interval geneInterval = new Interval(gtfFeature.getContig(), geneStart, geneEnd);
        GtfInfo gtfGeneInfo = new GtfInfo(geneInterval, GtfInfo.Type.GENE, gtfFeature.getGeneName());
        idToInfo.put(gtfFeature.getGeneId(), gtfGeneInfo);
    }

    // runs immediately after it has gone through each line of gtf (apply method)
    @Override
    public Object onTraversalSuccess() {
        // get the user input dictionary
        SAMSequenceDictionary sequenceDictionary = getSequenceDictionary(String.valueOf(this.sequenceDictionary));

        // create linked hash map to store sorted values of idToInfo
        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = getSortedMap(sequenceDictionary);

        // if user wants to sort by transcript only use transcripts else only use genes
        GtfInfo.Type selectedType = sortByTranscript ? GtfInfo.Type.TRANSCRIPT : GtfInfo.Type.GENE;
        writeToBed(selectedType, karyotypeIdToInfo);

        return null;
    }

    // reads the dictionary
    private SAMSequenceDictionary getSequenceDictionary(String dictionaryPath) {
        SamReader reader = SamReaderFactory.makeDefault().open(new File(dictionaryPath)); //TODO: figure out what makeDefault() does

        return reader.getFileHeader().getSequenceDictionary();
    }

    // sorts the map containing the features based on contig and start position
    private LinkedHashMap<String, GtfInfo> getSortedMap(SAMSequenceDictionary sequenceDictionary) {
        // create a list that has the keys and values of idToInfo and sort the list using CompareGtfInfo
        List<Map.Entry<String, GtfInfo>> entries = new ArrayList<>(idToInfo.entrySet());
        entries.sort(new CompareGtfInfo(sequenceDictionary));

        // put each (sorted) entry in the list into a linked hashmap
        LinkedHashMap<String, GtfInfo> karyotypeIdToInfo = new LinkedHashMap<>();
        for (Map.Entry<String, GtfInfo> entry : entries) {
            karyotypeIdToInfo.put(entry.getKey(), entry.getValue());
        }

        return karyotypeIdToInfo;
    }

    // writes to bed file
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

    // formats each line of the bed file depending on whether user has selected gene or transcript
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

