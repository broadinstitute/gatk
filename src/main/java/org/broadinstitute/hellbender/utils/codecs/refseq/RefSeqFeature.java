package org.broadinstitute.hellbender.utils.codecs.refseq;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;
import java.util.List;

/**
 * The ref seq feature. See how to generate these here: https://gatkforums.broadinstitute.org/gatk/discussion/1329/where-can-i-get-a-gene-list-in-refseq-format
 */
//TODO if there is cause to use this class for any purpose other than DepthOfCoverage then we should port remaining functionality from GATK3
public class RefSeqFeature implements RefSeqTranscript, Feature {

    private String transcript_id;
    private int strand;
    private SimpleInterval transcript_interval;
    private SimpleInterval transcript_coding_interval;
    private List<SimpleInterval> exons;
    private String gene_name;
    private List<Integer> exon_frames;
    private String name;

    public RefSeqFeature(SimpleInterval genomeLoc) {
        this.transcript_interval = genomeLoc;
    }

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId() { return transcript_id; }

    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand() { return strand; }

    @Override
    public SimpleInterval getLocation() {
        return transcript_interval;
    }

    /** Returns genomic interval of the coding sequence (does not include UTRs, but still includes introns, since it's a single interval on the DNA) */
    public SimpleInterval getCodingLocation() { return transcript_coding_interval; }

    /** Name of the gene this transcript corresponds to (NOT gene id such as Entrez etc) */
    public String getGeneName() { return gene_name; }

    /** Number of exons in this transcript */
    public int getNumExons() { return exons.size(); }

    /** Genomic location of the n-th exon; throws an exception if n is out of bounds */
    public SimpleInterval getExonLocation(int n) {
        if ( n >= exons.size() || n < 0 ) {
            throw new GATKException("Index out-of-bounds. RefSeqTranscript has " + exons.size() +" exons; requested: "+n);
        }
        return exons.get(n);
    }

    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<SimpleInterval> getExons() { return Collections.unmodifiableList(exons); }

    /** Returns a uniquified name of the Gene and TranscriptID*/
    public String getTranscriptUniqueGeneName() {
        return String.format("%s(%s)",getGeneName(),getTranscriptId());
    }

    /**
     * Returns a count of the total number of reference bases spanned by gene summary. Will total the length of the exons or
     * if absent, the lengthOnReference for the gene itself.
     *
     * NOTE: This currently makes the assumption that genes do not ever have overlapping exons. I do not know if this is a fair
     *       assumption given extant RefSeqGeneList files.
     */
    public int getTotalExonLength() {
        if (exons.isEmpty()) {
            return getLengthOnReference();
        }
        return exons.stream().mapToInt(Locatable::getLengthOnReference).sum();
    }

    /**
     * Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript.
     *
     * NOTE: this is is checking that the locatable is entirely contained within at least one exon.
     * */
    @Override
    public boolean contains(Locatable that) {
        if (exons.isEmpty()) {
            return getLocation().contains(that);
        }
        for ( SimpleInterval exon : exons ) {
            if ( IntervalUtils.overlaps(exon, that) ) {
                return true;
            }
        }
        return false;
    }

    public String toString() {
            StringBuilder b = new StringBuilder("000\t"); // first field is unused but required in th ecurrent format; just set to something
            b.append(transcript_id);   // #1
            b.append('\t');
            b.append(getLocation().getContig()); // #2
            b.append('\t');
            b.append( (strand==1?'+':'-') ); // #3
            b.append('\t');
            b.append( (getLocation().getStart() - 1) ); // #4
            b.append('\t');
            b.append( getLocation().getEnd());  // #5
            b.append('\t');
            b.append( (transcript_coding_interval.getStart() - 1) ); // #6
            b.append('\t');
            b.append( transcript_coding_interval.getEnd());  // #7
            b.append('\t');
            b.append(exons.size()); // #8
            b.append('\t');
            for ( SimpleInterval loc : exons ) { b.append( (loc.getStart()-1) ); b.append(','); } // #9
            b.append('\t');
            for ( SimpleInterval loc : exons ) { b.append( loc.getEnd() ); b.append(','); } // #10
            b.append("\t0\t"); // # 11 - unused?
            b.append(gene_name); // # 12
            b.append("\tcmpl\tcmpl\t"); // #13, #14 - unused?
            for ( Integer f : exon_frames ) { b.append( f ); b.append(','); } // #15


            return b.toString();
        }

    public void setTranscript_id(String transcript_id) {
        this.transcript_id = transcript_id;
    }

    public void setStrand(int strand) {
        this.strand = strand;
    }

    public void setTranscript_interval(SimpleInterval transcript_interval) {
        this.transcript_interval = transcript_interval;
    }

    public void setTranscript_coding_interval(SimpleInterval transcript_coding_interval) {
        this.transcript_coding_interval = transcript_coding_interval;
    }

    public void setExons(List<SimpleInterval> exons) {
        this.exons = exons;
    }

    public void setGene_name(String gene_name) {
        this.gene_name = gene_name;
    }

    public void setExon_frames(List<Integer> exon_frames) {
        this.exon_frames = exon_frames;
    }

    public void setName(String name) {
        this.name = name;
    }

    @Override
    public String getContig() {
        return transcript_interval.getContig();
    }

    @Override
    public int getStart() {
        return transcript_interval.getStart();
    }

    @Override
    public int getEnd() {
        return transcript_interval.getEnd();
    }
}
