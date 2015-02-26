/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
