/*
 * $Id: NumberList.java 70611 2008-07-29 18:06:55Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * A series of numbers.
 *
 * @author tsharpe
 * @version $Revision: 70611 $
 */
public interface NumberList
{
    int size();
    int intVal( int idx );
    double dblVal( int idx );
    Iterator iterator();

    interface Iterator
    {
        boolean hasNext();
        int nextInt();
        double nextDbl();
    }
}
