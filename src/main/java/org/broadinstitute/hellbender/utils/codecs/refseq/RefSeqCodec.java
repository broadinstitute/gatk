package org.broadinstitute.hellbender.utils.codecs.refseq;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;

/**
 * Allows for reading in RefSeq information
 * TODO this header needs to be rewritten
 * <p>
 * Parses a sorted UCSC RefSeq file (see below) into relevant features: the gene name, the unique gene name (if multiple transcrips get separate entries), exons, gene start/stop, coding start/stop,
 * strandedness of transcription. 
 * </p>
 *
 * <p>
 * Instructions for generating a RefSeq file for use with the RefSeq codec can be found on the documentation guide here
 * <a href="https://gatkforums.broadinstitute.org/gatk/discussion/1329/where-can-i-get-a-gene-list-in-refseq-format">https://gatkforums.broadinstitute.org/gatk/discussion/1329/where-can-i-get-a-gene-list-in-refseq-format</a>
 * </p>
 *
 * You will need to consult individual walkers for the argument name
 *
 * <h2>File format example</h2>
 * If you want to define your own file for use, the format is (tab delimited):
 * bin, name, chrom, strand, transcription start, transcription end, coding start, coding end, num exons, exon starts, exon ends, id, alt. name, coding start status (complete/incomplete), coding end status (complete,incomplete)
 * and exon frames, for example:
 * <pre>
 * 76 NM_001011874 1 - 3204562 3661579 3206102 3661429 3 3204562,3411782,3660632, 3207049,3411982,3661579, 0 Xkr4 cmpl cmpl 1,2,0,
 * </pre>
 * for more information see <a href="http://skip.ucsc.edu/cgi-bin/hgTables?hgsid=5651&hgta_doSchemaDb=mm8&hgta_doSchemaTable=refGene">here</a>
 * <p>
 *     
 * </p>
 *
 * @since 2010
 */
public class RefSeqCodec extends AsciiFeatureCodec<RefSeqFeature>{

    // codec file extension
    protected static final String FILE_EXT = "refseq";
    public static final String COMMENT_LINE_CHARACTER = "#";
    public static final String LINE_DELIMITER = "\t";
    public static final int MINIMUM_LINE_FIELD_COUNT = 16;
    private boolean zero_coding_length_user_warned = false;

    // Indexes for interpreting lines:
    public static final int TRANSCRIPT_ID_INDEX = 1;
    public static final int CONTIG_INDEX = 2;
    public static final int STRAND_INDEX = 3;
    public static final int INTERVAL_LEFT_BOUND_INDEX = 4;
    public static final int INTERVAL_RIGHT_BOUND_INDEX = 5;
    public static final int CODING_START_INDEX = 6;
    public static final int CODING_STOP_INDEX = 7;
    public static final int EXON_STARTS_INDEX = 9;
    public static final int EXON_STOPS_INDEX = 10;
    public static final int GENE_NAME_INDEX = 12;
    public static final int EXON_FRAMES_INDEX = 15;


    public RefSeqCodec() {
        super(RefSeqFeature.class);
    }

    @Override
    public Feature decodeLoc(final LineIterator lineIterator) {
        final String line = lineIterator.next();
        if (line.startsWith(COMMENT_LINE_CHARACTER)){
            return null;
        }
        final String fields[] = line.split(LINE_DELIMITER);
        if (fields.length < MINIMUM_LINE_FIELD_COUNT){
            throw new TribbleException("RefSeq (decodeLoc) : Unable to parse line -> " + line + ", we expected at least 16 columns, we saw " + fields.length);
        }
        final String contig_name = fields[CONTIG_INDEX];
        try {
            return new RefSeqFeature(new SimpleInterval(contig_name, Integer.parseInt(fields[INTERVAL_LEFT_BOUND_INDEX])+1, Integer.parseInt(fields[INTERVAL_RIGHT_BOUND_INDEX])));
        //TODO maybe except for malformed simple intervals? Genome locs had that
        } catch ( NumberFormatException e ) {
            throw new UserException.MalformedFile("Could not parse location from line: " + line);
        }
    }

    /** Fills this object from a text line in RefSeq (UCSC) text dump file */
    @Override
    public RefSeqFeature decode(String line) {
        if (line.startsWith(COMMENT_LINE_CHARACTER)) {
            return null;
        }
        String fields[] = line.split(LINE_DELIMITER);

        // we reference postion 15 in the split array below, make sure we have at least that many columns
        if (fields.length < MINIMUM_LINE_FIELD_COUNT) {
            throw new TribbleException("RefSeq (decode) : Unable to parse line -> " + line + ", we expected at least 16 columns, we saw " + fields.length);
        }
        String contig_name = fields[CONTIG_INDEX];
        RefSeqFeature feature = new RefSeqFeature(new SimpleInterval(contig_name, Integer.parseInt(fields[INTERVAL_LEFT_BOUND_INDEX])+1, Integer.parseInt(fields[INTERVAL_RIGHT_BOUND_INDEX])));

        feature.setTranscript_id(fields[TRANSCRIPT_ID_INDEX]);
        if ( fields[STRAND_INDEX].length()==1 && fields[STRAND_INDEX].charAt(0)=='+') {
            feature.setStrand(1);

        } else if ( fields[STRAND_INDEX].length()==1 && fields[STRAND_INDEX].charAt(0)=='-') {
            feature.setStrand(-1);

        } else {
            throw new UserException.MalformedFile("Expected strand symbol (+/-), found: "+fields[STRAND_INDEX] + " for line=" + line);
        }

        int coding_start = Integer.parseInt(fields[CODING_START_INDEX])+1;
        int coding_stop = Integer.parseInt(fields[CODING_STOP_INDEX]);

        if ( coding_start > coding_stop ) {
            if ( ! zero_coding_length_user_warned ) {
                Utils.warnUser("RefSeq file contains transcripts with zero coding length. "+
                        "Such transcripts will be ignored (this warning is printed only once)");
                zero_coding_length_user_warned = true;
            }
            return null;
        }

        feature.setTranscript_interval(new SimpleInterval(contig_name, Integer.parseInt(fields[INTERVAL_LEFT_BOUND_INDEX])+1, Integer.parseInt(fields[INTERVAL_RIGHT_BOUND_INDEX])));
        feature.setTranscript_coding_interval(new SimpleInterval(contig_name, coding_start, coding_stop));
        feature.setGene_name(fields[GENE_NAME_INDEX]);
        String[] exon_starts = fields[EXON_STARTS_INDEX].split(",");
        String[] exon_stops = fields[EXON_STOPS_INDEX].split(",");
        String[] eframes = fields[EXON_FRAMES_INDEX].split(",");

        if ( exon_starts.length != exon_stops.length ) {
            throw new UserException.MalformedFile("Data format error: numbers of exon start and stop positions differ for line=" + line);
        }

        if ( exon_starts.length != eframes.length ) {
            throw new UserException.MalformedFile("Data format error: numbers of exons and exon frameshifts differ for line=" + line);
        }

        ArrayList<SimpleInterval> exons = new ArrayList<>(exon_starts.length);
        ArrayList<Integer> exon_frames = new ArrayList<Integer>(eframes.length);

        for ( int i = 0 ; i < exon_starts.length  ; i++ ) {
            // NOTE, we add 1 here to account for represention issues in the exon counts
            exons.add(new SimpleInterval(contig_name, Integer.parseInt(exon_starts[i])+1, Integer.parseInt(exon_stops[i]) ) );
            exon_frames.add(Integer.decode(eframes[i]));
        }

        feature.setExons(exons);
        feature.setExon_frames(exon_frames);
        return feature;
    }

    /**
     * Can the file be decoded?
     * @param path path the file to test for parsability with this codec
     * @return true if the path has the correct file extension, false otherwise
     */
    @Override
    public boolean canDecode(final String path) { return path.endsWith("." + FILE_EXT); }

    @Override
    public Object readActualHeader(LineIterator lineIterator) {
        // No header for this format
        return null;
    }
}
