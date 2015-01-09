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

package org.broadinstitute.hellbender.utils.variant;

/**
 * <i>Homo sapiens</i> genome constants.
 *
 * <p>NOTE: reference to these constants is an indication that your code is (human) species assumption dependant.</p>
 */
public class HomoSapiensConstants {

    /**
     * Standard heterozygous rate for SNP variation.
     */
    public static final double SNP_HETEROZYGOSITY = 1e-3;

    /**
     * Standard heterozygous rate for INDEL variation.
     */
    public static final double INDEL_HETEROZYGOSITY = 1.0/8000;

    /**
     * Standard ploidy for autosomal chromosomes.
     */
    public static final int DEFAULT_PLOIDY = 2;
}
