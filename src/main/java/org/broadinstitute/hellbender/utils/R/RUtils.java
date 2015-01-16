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

package org.broadinstitute.hellbender.utils.R;

import org.apache.commons.lang3.StringUtils;

import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;

public class RUtils {
    /**
     * Converts a collection of values to an R compatible list. A null list will return NA,
     * otherwise the values will be escaped with single quotes and combined with c().
     * @param list Collection of values
     * @return The R representation of the list
     */
    public static String toStringList(Collection<? extends CharSequence> list) {
        if (list == null)
            return "NA";
        if (list.size() == 0)
            return "c()";
        return "c('" + StringUtils.join(list, "','") + "')";
    }

    /**
     * Converts a collection of values to an R compatible list. A null list will return NA,
     * otherwise the values will be combined with c().
     * @param list Collection of values
     * @return The R representation of the list
     */
    public static String toNumberList(Collection<? extends Number> list) {
        return list == null ? "NA": "c(" + StringUtils.join(list, ",") + ")";
    }

    /**
     * Converts a collection of values to an R compatible list. A null list will return NA,
     * otherwise the date will be escaped with single quotes and combined with c().
     * @param list Collection of values
     * @return The R representation of the list
     */
    public static String toDateList(Collection<? extends Date> list) {
        return toDateList(list, "''yyyy-MM-dd''");
    }

    /**
     * Converts a collection of values to an R compatible list formatted by pattern.
     * @param list Collection of values
     * @param pattern format pattern string for each date
     * @return The R representation of the list
     */
    public static String toDateList(Collection<? extends Date> list, String pattern) {

        if (list == null)
            return "NA";
        SimpleDateFormat format = new SimpleDateFormat(pattern);
        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        boolean first = true;
        for (Date date : list) {
            if (!first) sb.append(",");
            sb.append(format.format(date));
            first = false;
        }
        sb.append(")");
        return sb.toString();
    }
}
