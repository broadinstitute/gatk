package org.broadinstitute.hellbender.utils.codecs.refseq;

import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.List;

public interface Transcript {

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId();
    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand();
    /** Returns transcript's full genomic interval (includes all exons with UTRs) */
    public GenomeLoc getLocation();
    /** Returns genomic interval of the coding sequence (does not include
     * UTRs, but still includes introns, since it's a single interval on the DNA)
     */
    public GenomeLoc getCodingLocation();
    /** Name of the gene this transcript corresponds to (typically NOT gene id such as Entrez etc,
     * but the implementation can decide otherwise)
     */
    public String getGeneName();
    /** Number of exons in this transcript */
    public int getNumExons();
    /** Genomic location of the n-th exon; expected to throw an exception (runtime) if n is out of bounds */
    public GenomeLoc getExonLocation(int n);

    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<GenomeLoc> getExons();

    /** Returns true if the specified interval 'that' overlaps with the full genomic interval of this transcript */
    public boolean overlapsP (GenomeLoc that);

    /** Returns true if the specified interval 'that' overlaps with the coding genomic interval of this transcript.
      * NOTE: since "coding interval" is still a single genomic interval, it will not contain UTRs of the outermost exons,
      * but it will still contain introns and/or exons internal to this genomic locus that are not spliced into this transcript.
      * @see #overlapsExonP
      */
    public boolean overlapsCodingP (GenomeLoc that);

    /** Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript */
    public boolean overlapsExonP (GenomeLoc that);


}
