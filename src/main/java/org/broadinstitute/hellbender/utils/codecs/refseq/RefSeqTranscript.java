package org.broadinstitute.hellbender.utils.codecs.refseq;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 22, 2009
 * Time: 5:22:30 PM
 * To change this template use File | Settings | File Templates.
 */
// TODO if we want to expand transcript functionality beyond what is used for DepthOfCoverage there is more functionality that we want to port from GATK3
public interface RefSeqTranscript extends Locatable {

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId();
    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand();
    /** Returns transcript's full genomic interval (includes all exons with UTRs) */
    public SimpleInterval getLocation();
    /** Returns genomic interval of the coding sequence (does not include
     * UTRs, but still includes introns, since it's a single interval on the DNA)
     */
    public SimpleInterval getCodingLocation();
    /** Name of the gene this transcript corresponds to (typically NOT gene id such as Entrez etc,
     * but the implementation can decide otherwise)
     */
    public String getGeneName();
    /** Number of exons in this transcript */
    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<SimpleInterval> getExons();

}
