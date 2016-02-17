/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.contexts;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

/**
 * The section of the reference that overlaps with the given
 * read / locus. 
 *
 * @author hanna
 * @version 0.1
 */
public class ReferenceContext {
    /**
     * Facilitates creation of new GenomeLocs.
     */
    final private GenomeLocParser genomeLocParser;

    /**
     * The locus.
     */
    final private GenomeLoc locus;

    /**
     * The window of reference information around the current locus.
     */
    final private GenomeLoc window;

    /**
     * The bases in the window around the current locus.  If null, then bases haven't been fetched yet.
     * Bases are always upper cased
     */
    private byte[] basesCache = null;

    /**
     * Lazy loader to fetch reference bases
     */
    final private ReferenceContextRefProvider basesProvider;

    /**
     * Interface to create byte[] contexts for lazy loading of the reference
     */
    public static interface ReferenceContextRefProvider {
        /**
         * You must provide a routine that gets the byte[] bases that would have been passed into the
         * ReferenceContext.  The RC will handling caching.  The value of this interface and routine is
         * that it is only called when the bytes are actually requested by the walker, not up front.  So
         * if the walker doesn't need the refBases for whatever reason, there's no overhead to
         * provide them.
         *
         * @return
         */
        @Ensures({"result != null"})
        public byte[] getBases();
    }

    private static class ForwardingProvider implements ReferenceContextRefProvider {
        byte[] bases;

        public ForwardingProvider( byte base ) {
            this(new byte[] { base });
        }

        public ForwardingProvider( byte[] bases ) {
            this.bases = bases;
        }

        public byte[] getBases() { return bases; }
    }

    /**
     * Contructor for a simple, windowless reference context.
     * @param locus locus of interest.
     * @param base reference base at that locus.
     */
    @Requires({
            "genomeLocParser != null",
            "locus != null",
            "locus.size() > 0"})
    public ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, byte base ) {
        this( genomeLocParser, locus, locus, new ForwardingProvider(base) );
    }

    @Requires({
            "genomeLocParser != null",
            "locus != null",
            "locus.size() > 0",
            "window != null",
            "window.size() > 0",
            "bases != null && bases.length > 0"})
    public ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, GenomeLoc window, byte[] bases ) {
        this( genomeLocParser, locus, window, new ForwardingProvider(bases) );
    }

    @Requires({
            "genomeLocParser != null",
            "locus != null",
            "locus.size() > 0",
            "window != null",
            "window.size() > 0",
            "basesProvider != null"})
    public ReferenceContext( GenomeLocParser genomeLocParser, GenomeLoc locus, GenomeLoc window, ReferenceContextRefProvider basesProvider ) {
        this.genomeLocParser = genomeLocParser;
        this.locus = locus;
        this.window = window;
        this.basesProvider = basesProvider;
    }

    /**
     * Utility function to load bases from the provider to the cache, if necessary
     */
    @Ensures({
            "basesCache != null",
            "old(basesCache) == null || old(basesCache) == basesCache"})
    private void fetchBasesFromProvider() {
        if ( basesCache == null ) {
            basesCache = basesProvider.getBases();

            // must be an assertion that only runs when the bases are fetch to run in a reasonable amount of time
            assert BaseUtils.isUpperCase(basesCache);
        }
    }

    /**
     * @return The genome loc parser associated with this reference context
     */
    @Ensures("result != null")
    public GenomeLocParser getGenomeLocParser() {
        return genomeLocParser;
    }

    /**
     * The locus currently being examined.
     * @return The current locus.
     */
    @Ensures("result != null")
    public GenomeLoc getLocus() {
        return locus;
    }

    @Ensures("result != null")
    public GenomeLoc getWindow() {
        return window;
    }

    /**
     * Get the base at the given locus.
     * @return The base at the given locus from the reference.
     */
    public byte getBase() {
        return getBases()[(locus.getStart() - window.getStart())];
    }

    /**
     * All the bases in the window currently being examined.
     * @return All bases available.  If the window is of size [0,0], the array will
     *         contain only the base at the given locus.
     */
    @Ensures({"result != null", "result.length > 0"})
    public byte[] getBases() {
        fetchBasesFromProvider();
        return basesCache;
    }

    /**
     * All the bases in the window from the current base forward to the end of the window.
     */
    @Ensures({"result != null", "result.length > 0"})
    public byte[] getForwardBases() {
        final byte[] bases = getBases();
        final int mid = locus.getStart() - window.getStart();
        // todo -- warning of performance problem, especially if this is called over and over
        return new String(bases).substring(mid).getBytes();
    }

    @Deprecated
    public char getBaseAsChar() {
        return (char)getBase();
    }

    /**
     * Get the base at the given locus.
     * @return The base at the given locus from the reference.
     */
    @Deprecated()
    public int getBaseIndex() {
        return BaseUtils.simpleBaseToBaseIndex(getBase());
    }
}
