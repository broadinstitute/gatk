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

package org.broadinstitute.hellbender.utils.codecs.hapmap;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;

import java.io.IOException;
import java.util.Arrays;

/**
 * A codec for the file types produced by the HapMap consortium
 *
 * <p>
 *     The format includes eleven standard fields, plus genotypes for each of the samples included
 *     in the file:
 *
 * <pre>
 *     Col1: refSNP rs# identifier at the time of release (NB might merge with another rs# in the future)
 *     Col2: SNP alleles according to dbSNP
 *     Col3: chromosome that SNP maps to
 *     Col4: chromosome position of SNP, in basepairs on reference sequence
 *     Col5: strand of reference sequence that SNP maps to
 *     Col6: version of reference sequence assembly
 *     Col7: HapMap genotype center that produced the genotypes
 *     Col8: LSID for HapMap protocol used for genotyping
 *     Col9: LSID for HapMap assay used for genotyping
 *     Col10: LSID for panel of individuals genotyped
 *     Col11: QC-code, currently 'QC+' for all entries (for future use)
 *     Col12 and on: observed genotypes of samples, one per column, sample identifiers in column headers (Coriell catalog numbers, example: NA10847). Duplicate samples have .dup suffix.
 * </pre>
 * </p>
 *
 * <p>
 *  See also: @See <a href="http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/">HapMap genotypes download</a>
 * </p>
 *
 * <h2>File format example</h2>
 * From <a href="http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/latest/forward/non-redundant/genotypes_chr1_ASW_r27_nr.b36_fwd.txt.gz">genotypes_chr1_ASW_r27_nr.b36_fwd.txt.gz</a>:
 * <pre>
 *     rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode NA19625 NA19700 NA19701 NA19702 NA19703 NA19704 NA19705 NA19708 NA19712 NA19711 NA19818 NA19819 NA19828 NA19835 NA19834 NA19836 NA19902 NA19901 NA19900 NA19904 NA19919 NA19908 NA19909 NA19914 NA19915 NA19916 NA19917 NA19918 NA19921 NA20129 NA19713 NA19982 NA19983 NA19714 NA19985 NA20128 NA20126 NA20127 NA20277 NA20276 NA20279 NA20282 NA20281 NA20284 NA20287 NA20288 NA20290 NA20289 NA20291 NA20292 NA20295 NA20294 NA20297 NA20300 NA20301 NA20302 NA20317 NA20319 NA20322 NA20333 NA20332 NA20335 NA20334 NA20337 NA20336 NA20340 NA20341 NA20343 NA20342 NA20344 NA20345 NA20346 NA20347 NA20348 NA20349 NA20350 NA20357 NA20356 NA20358 NA20359 NA20360 NA20363 NA20364
 *     rs9629043 C/T chr1 554636 + ncbi_b36 broad urn:LSID:affymetrix.hapmap.org:Protocol:GenomeWideSNP_6.0:3 urn:LSID:broad.hapmap.org:Assay:SNP_A-8575115:3 urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:3 QC+ CC CC CC CC CC CC CC CC CC CC CC CC NN CC CC CC CT CT CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CT CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC CC
 *     rs28446478 G/T chr1 576058 + ncbi_b36 sanger urn:LSID:illumina.hapmap.org:Protocol:Human_1M_BeadChip:3 urn:LSID:sanger.hapmap.org:Assay:H1Mrs28446478:3 urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:3 QC+ GT TT GT TT TT TT TT GT GT TT TT TT TT GT GT GT GT TT GT TT GT GT TT GT GT TT TT TT GT GT TT TT TT GT TT GT TT GT GT GT GT GT TT GT TT TT GT GT TT TT TT TT TT TT GT GT GT GT TT TT TT TT GT TT GT TT TT GT TT TT TT GT TT TT TT GT GT TT GT TT GT TT TT
 *     rs12565286 C/G chr1 711153 + ncbi_b36 broad urn:LSID:affymetrix.hapmap.org:Protocol:GenomeWideSNP_6.0:3 urn:LSID:broad.hapmap.org:Assay:SNP_A-8709646:3 urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:3 QC+ GG GG GG GG GG GG GG GG CG GG GG GG GG GG GG GG GG GG GG CG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG CG GG GG GG GG GG GG GG CG CG GG GG GG GG GG GG GG GG GG CG CG GG GG GG GG GG GG GG GG GG GG CG NN GG GG GG GG GG GG NN GG NN NN
 * </pre>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class RawHapMapCodec extends AsciiFeatureCodec<RawHapMapFeature> {
    // the minimum number of features in the HapMap file line
    private static final int minimumFeatureCount = 11;

    private String headerLine;

    public RawHapMapCodec() {
        super(RawHapMapFeature.class);
    }

    /**
     * decode the hapmap record
     * @param line the input line to decode
     * @return a HapMapFeature, with the given fields 
     */
    public RawHapMapFeature decode(String line) {
        String[] array = line.split("\\s+");

        // make sure the split was successful - that we got an appropriate number of fields
        if (array.length < minimumFeatureCount)
            throw new IllegalArgumentException("Unable to parse line " + line + ", the length of split features is less than the minimum of " + minimumFeatureCount);

        // create a new feature given the array
        return new RawHapMapFeature(array[0],
                array[1].split("/"),
                array[2],
                Long.valueOf(array[3]),
                Strand.toStrand(array[4]),
                array[5],
                array[6],
                array[7],
                array[8],
                array[9],
                array[10],
                Arrays.copyOfRange(array,11,array.length),
                headerLine);
    }

    @Override
    public Object readActualHeader(final LineIterator lineIterator) {
        this.headerLine = lineIterator.next();
        return headerLine;
    }

    @Override
    public FeatureCodecHeader readHeader(final LineIterator lineIterator) throws IOException {
        final String header = (String) readActualHeader(lineIterator);
        // TODO: This approach may cause issues with files formatted with \r\n-style line-endings.
        return new FeatureCodecHeader(header, header.length() + 1);
    }

    @Override
    public boolean canDecode(final String path) {
        return path.toLowerCase().endsWith(".hapmap");
    }
}
