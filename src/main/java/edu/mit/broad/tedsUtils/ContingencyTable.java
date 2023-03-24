/*
 * $Id: ContingencyTable.java 51142 2007-11-05 17:19:49Z tsharpe $
 * WHITEHEAD INSTITUTE 
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2002 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * A 2-dimensional table of double values.
 * The ChiSquare implementation gets its data from this interface so that
 * you can adapt any source of data to use as its input.
 *
 * @author Ted Sharpe
 * @version $Revision: 11017 $
 */
public interface ContingencyTable
{
    int getNX();
    int getNY();
    double get(int x, int y);
}
