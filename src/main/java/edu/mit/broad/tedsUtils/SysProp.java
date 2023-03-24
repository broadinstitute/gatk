/*
 * $Id: SysProp.java 38013 2007-04-11 21:58:13Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * Utility methods for getting system properties.
 *
 * @author tsharpe
 * @version $Revision: 38013 $
 */
public class SysProp
{
    public static String sVal( String propName, String dflt )
    {
        return System.getProperty(propName,dflt);
    }

    public static boolean bVal( String propName, boolean dflt )
    {
        return System.getProperty(propName) == null ? dflt : Boolean.getBoolean(propName);
    }

    public static int iVal( String propName, int dflt )
    {
        return Integer.getInteger(propName,dflt).intValue();
    }

    public static double dVal( String propName, double dflt )
    {
        String propVal = System.getProperty(propName);
        return propVal == null ? dflt : Double.parseDouble(propVal);
    }

    public static float fVal( String propName, float dflt )
    {
        String propVal = System.getProperty(propName);
        return propVal == null ? dflt : Float.parseFloat(propVal);
    }
}
