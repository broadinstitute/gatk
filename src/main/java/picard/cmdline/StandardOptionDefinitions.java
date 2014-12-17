/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.cmdline;

/**
 * A set of String constants in which the name of the constant (minus the _SHORT_NAME suffix)
 * is the standard long Option name, and the value of the constant is the standard shortName.
 */
public class StandardOptionDefinitions {
    public static final String INPUT_SHORT_NAME = "I";
    public static final String OUTPUT_SHORT_NAME = "O";
    public static final String REFERENCE_SHORT_NAME = "R";
    public static final String SAMPLE_ALIAS_SHORT_NAME = "ALIAS";
    public static final String LIBRARY_NAME_SHORT_NAME = "LIB";
    public static final String EXPECTED_INSERT_SIZE_SHORT_NAME = "INSERT";
    public static final String LANE_SHORT_NAME = "L";
    public static final String SEQUENCE_DICTIONARY_SHORT_NAME = "SD";
    public static final String METRICS_FILE_SHORT_NAME = "M";
    public static final String ASSUME_SORTED_SHORT_NAME = "AS";
    public static final String PF_READS_ONLY_SHORT_NAME = "PF";
    public static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "MQ";
    public static final String READ_GROUP_ID_SHORT_NAME = "RG";
    public static final String PROGRAM_RECORD_ID_SHORT_NAME = "PG";
    public static final String MINIMUM_LOD_SHORT_NAME = "LOD";
    public static final String SORT_ORDER_SHORT_NAME = "SO";
    public static final String USE_ORIGINAL_QUALITIES_SHORT_NAME = "OQ";
}
