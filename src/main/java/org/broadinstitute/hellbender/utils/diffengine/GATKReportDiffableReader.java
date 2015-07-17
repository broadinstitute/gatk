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

package org.broadinstitute.hellbender.utils.diffengine;

import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;


/**
 * Class implementing diffnode reader for GATKReports
 */

// TODO Version check to be added at the report level

public class GATKReportDiffableReader implements DiffableReader {
    @Override
    public String getName() {
        return "GATKReport";
    }

    @Override
    public DiffElement readFromFile(File file, int maxElementsToRead) {
        DiffNode root = DiffNode.rooted(file.getName());
        try {
            // one line reads the whole thing into memory
            GATKReport report = new GATKReport(file);

            for (GATKReportTable table : report.getTables()) {
                root.add(tableToNode(table, root));
            }

            return root.getBinding();
        } catch (Exception e) {
            return null;
        }
    }

    private DiffNode tableToNode(GATKReportTable table, DiffNode root) {
        DiffNode tableRoot = DiffNode.empty(table.getTableName(), root);

        tableRoot.add("Description", table.getTableDescription());
        tableRoot.add("NumberOfRows", table.getNumRows());

        for ( GATKReportColumn column : table.getColumnInfo() ) {
            DiffNode columnRoot = DiffNode.empty(column.getColumnName(), tableRoot);

            columnRoot.add("Width", column.getColumnFormat().getWidth());
            // NOTE: as the values are trimmed during parsing left/right alignment is not currently preserved
            columnRoot.add("Displayable", true);

            for ( int i = 0; i < table.getNumRows(); i++ ) {
                String name = column.getColumnName() + (i+1);
                columnRoot.add(name, table.get(i, column.getColumnName()).toString());
            }

            tableRoot.add(columnRoot);
        }

        return tableRoot;
    }

    @Override
    public boolean canRead(File file) {
        try {
            final String HEADER = GATKReport.GATKREPORT_HEADER_PREFIX;
            final char[] buff = new char[HEADER.length()];
            final FileReader FR = new FileReader(file);
            FR.read(buff, 0, HEADER.length());
            FR.close();
            String firstLine = new String(buff);
            return firstLine.startsWith(HEADER);
        } catch (IOException e) {
            return false;
        }
    }
}
