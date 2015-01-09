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

package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.exceptions.UserException;

public enum GATKReportVersion {
    /**
     * Differences between other versions:
     * - Does not allow spaces in cells.
     * - Mostly fixed width but has a bug where the string width of floating point
     * values was not measured correctly leading to columns that aren't aligned
     */
    V0_1("v0.1"),

    /**
     * Differences between other versions:
     * - Spaces allowed in cells, for example in sample names with spaces in them ex: "C507/FG-CR 6".
     * - Fixed width fixed for floating point values
     */
    V0_2("v0.2"),

    /*
    * Differences between v0.x
    * - Added table and report headers
    * - Headers changed format, include the number of tables, rows, and metadata for gathering
    * - IS GATHERABLE
    */
    V1_0("v1.0"),

    /*
    * Differences between v1.0
    * - column numbers in header reflect the actual count of columns
    * - primary keys are never displayed
    */
    V1_1("v1.1");

    private final String versionString;

    private GATKReportVersion(String versionString) {
        this.versionString = versionString;
    }

    @Override
    public String toString() {
        return versionString;
    }

    public boolean equals(GATKReportVersion that) {
        return (versionString.equals(that.versionString));
    }

    /**
     * Returns the GATK Report Version from the file header.
     *
     * @param header Header from the file starting with ##:GATKReport.v[version]
     * @return The version as an enum.
     */
    public static GATKReportVersion fromHeader(String header) {
        if ( header == null )
            throw new UserException.BadInput("The GATK report has no version specified in the header");

        if (header.startsWith("##:GATKReport.v0.1 "))
            return GATKReportVersion.V0_1;

        if (header.startsWith("##:GATKReport.v0.2 "))
            return GATKReportVersion.V0_2;

        if (header.startsWith("#:GATKReport.v1.0"))
            return GATKReportVersion.V1_0;

        if (header.startsWith("#:GATKReport.v1.1"))
            return GATKReportVersion.V1_1;

        throw new UserException.BadInput("The GATK report has an unknown/unsupported version in the header: " + header);
    }
}
