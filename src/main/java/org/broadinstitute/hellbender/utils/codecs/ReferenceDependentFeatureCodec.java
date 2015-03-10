package org.broadinstitute.hellbender.utils.codecs;

import org.broadinstitute.hellbender.utils.GenomeLocParser;

/**
 * An interface marking that a given Tribble feature/codec is actually dependent on context within the
 * reference, rather than having a dependency only on the contig, start, and stop of the given feature.
 * A HACK.  Tribble should contain all the information in needs to decode the unqualified position of
 * a feature.
 */
// TODO: All "reference-dependent" codecs are currently broken. Need to
// TODO: refactor them so that they don't require a GenomeLocParser
public interface ReferenceDependentFeatureCodec {
    /**
     * Sets the appropriate GenomeLocParser, providing additional context when decoding larger and more variable features.
     * @param genomeLocParser The parser to supply. 
     */
    public void setGenomeLocParser(GenomeLocParser genomeLocParser);
}
