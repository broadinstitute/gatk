/*
 * $Id: QualityTrimmer.java 49350 2007-10-04 20:48:40Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * Some methods for calculating quality clips.<br>
 * We're using Java-standard, 0-based indexing.<br>
 * The left quality clip methods give you the index of the first "good" score.<br>
 * You should trim away from 0 to this value to be left with the high-quality stuff.<br>
 * The right quality clip methods give you the index of the first "bad" score.<br>
 * You should trim away from this index to the end to be left with the high-quality stuff.<br>
 * @author tsharpe
 * @version $Revision: 37214 $
 */
public class QualityTrimmer
{
    /**
     * Compute the left quality clip according to the Jaffe criterion.<br>
     * This method returns the index into the quals array of the first score that starts a
     * group of 12 consecutive scores that have values no less than 13.<br>
     * It returns the length of the array if there is no such window.  (In other words, the
     * whole enchilada should be trimmed away.)
     * @param quals The quality scores.
     * @return Start of the "good" scores.
     */
    public static int computeJaffeLeftClip( byte[] quals )
    {
        return computeLeftClip(quals,JAFFE_SCORE,JAFFE_WINDOW);
    }

    /**
     * Compute the right quality clip according to the time-honored, Pipeline/Squid tradition
     * for calculating the read_length quality statistic.<br>
     * This method returns the end of the rightmost window in which 20 consecutive quality scores
     * have an average value of at least 15.<br>
     * It returns 0 if there is no such window.<br>
     * @param quals The quality scores.
     * @return Start of the "bad" scores.
     */

    public static int computeReadLengthRightClip( byte[] quals )
    {
        return computeRightAvgClip(quals,READLENGTH_SCORE,READLENGTH_WINDOW);
    }

    /**
     * Compute the average quality in the time-honored, Pipeline/Squid tradition for
     * calculating the avg_quality statistic.
     * @param quals The quality scores.
     * @return The mean of the quality values for positions 100 up to, but not including 300.
     */
    public static int computeAvgQuality( byte[] quals )
    {
        return computeAverageQuality(quals,AVG_QUALITY_START,AVG_QUALITY_END);
    }

    /**
     * Do these quality scores earn a sequencing_pass?
     * @param quals The quality scores.
     * @return True, if the computeAvgQuality >= 20.
     */
    public static boolean isSeqPass( byte[] quals )
    {
        return computeAvgQuality(quals) >= SEQPASS_AVG_QUALITY;
    }

    /**
     * Counts the number of bases analyzable for SNPs in the time-honored way.<br>
     * This is the number of scores of quality 20 or better that are flanked on
     * each side by five scores no worse than 15.
     * @param quals The quality scores.
     * @return The ANA_SNP count.
     */
    public static int computeAnaSNPCount( byte[] quals )
    {
        if ( quals == null )
            throw new IllegalArgumentException("quals array may not be null");

        int result = 0;
        int nnn = quals.length;
        if ( nnn >= ANASNP_WINDOW )
        {
            int found = 0;
            for ( int iii = 0; iii < nnn; ++iii )
            {
                if ( quals[iii] < ANASNP_FLANKING_SCORE )
                {
                    found = 0;
                }
                else if ( ++found >= ANASNP_WINDOW )
                {
                    if ( quals[iii-ANASNP_FLANKING_WINDOW] >= ANASNP_CENTRAL_SCORE )
                        ++result;
                }
            }
        }
        return result;
    }

    /**
     * Left quality clip according to a minimum-score criterion.<br>
     * This method returns the index into the quals array of the first score that starts a
     * group of windowLen consecutive scores that have values no less than minQual.<br>
     * It returns the length of the array if there is no such window.
     * @param quals The quality scores.
     * @param minQual The minimum acceptable quality score.
     * @param windowLen The number of consecutive scores.
     * value.
     * @return Start of the "good" scores.
     */
    public static int computeLeftClip( byte[] quals, int minQual, int windowLen )
    {
        if ( quals == null )
            throw new IllegalArgumentException("quals array may not be null");

        if ( windowLen <= 0 )
            throw new IllegalArgumentException("window length must be positive");

        int nnn = quals.length;
        int result = nnn;
        int found = 0;
        for ( int iii = 0; iii < nnn; ++iii )
        {
            if ( quals[iii] < minQual )
            {
                found = 0;
            }
            else if ( ++found == windowLen )
            {
                result = iii + 1 - windowLen;
                break;
            }
        }

        return result;
    }

    /**
     * Left quality clip according to an average-score criterion.<br>
     * This method returns the index into the quals array of the first score that starts a
     * group of windowLen consecutive scores that have an average (mean) value no less
     * than avgQual.
     * It returns the length of the array if there is no such window.
     * @param quals The quality scores.
     * @param avgQual The minimum acceptable average quality score.
     * @param windowLen The number of consecutive scores.
     * @return Start of the "good" scores.
     */
    public static int computeLeftAvgClip( byte[] quals, int avgQual, int windowLen )
    {
        if ( quals == null )
            throw new IllegalArgumentException("quals array may not be null");

        if ( windowLen <= 0 )
            throw new IllegalArgumentException("window length must be positive");

        int nnn = quals.length;
        int result = nnn;

        if ( nnn >= windowLen )
        {
            int sum = 0;
            int idx2 = 0;
            while ( idx2 < windowLen )
            {
                sum += quals[idx2++];
            }

            int idx1 = 0;
            int minSum = avgQual * windowLen;
            while ( true )
            {
                if ( sum >= minSum )
                {
                    result = idx1;
                    break;
                }
                if ( idx2 >= nnn )
                {
                    break;
                }
                sum -= quals[idx1++];
                sum += quals[idx2++];
            }
        }

        return result;
    }

    /**
     * Compute right clip according to a minimum-score criterion.<br>
     * This method returns the end of the rightmost group of windowLen consecutive scores that
     * have values no less than minQual.<br>
     * It returns 0 if there is no such window.
     * @param quals The quality scores.
     * @param minQual The minimum acceptable quality score.
     * @param windowLen The number of consecutive scores.
     * @return Start of the "bad" scores.
     */
    public static int computeRightClip( byte[] quals, int minQual, int windowLen )
    {
        if ( quals == null )
            throw new IllegalArgumentException("quals array may not be null");

        if ( windowLen <= 0 )
            throw new IllegalArgumentException("window length must be positive");

        int result = 0;
        int found = 0;

        for ( int iii = quals.length - 1; iii >= 0; --iii )
        {
            if ( quals[iii] < minQual )
            {
                found = 0;
            }
            else if ( ++found == windowLen )
            {
                result = iii + windowLen;
                break;
            }
        }

        return result;
    }

    /**
     * Compute right clip according to a minimum-score criterion.<br>
     * This method returns the end of the rightmost group of windowLen consecutive scores that
     * have an average (mean) value no less than avgQual.<br>
     * It returns 0 if there is no such window.
     * @param quals The quality scores.
     * @param avgQual The minimum acceptable average quality score.
     * @param windowLen The number of consecutive scores.
     * @return Start of the "bad" scores.
     */
    public static int computeRightAvgClip( byte[] quals, int avgQual, int windowLen )
    {
        if ( quals == null )
            throw new IllegalArgumentException("quals array may not be null");

        if ( windowLen <= 0 )
            throw new IllegalArgumentException("window length must be positive");

        int result = 0;
        int idx2 = quals.length;

        if ( idx2 >= windowLen )
        {
            int idx1 = idx2 - windowLen;
            int sum = 0;
            while ( idx2 > idx1 )
            {
                sum += quals[--idx2];
            }

            int minSum = avgQual * windowLen;
            idx2 = quals.length;
            while ( true )
            {
                if ( sum >= minSum )
                {
                    result = idx2;
                    break;
                }
                if ( idx1 <= 0 )
                {
                    break;
                }
                sum -= quals[--idx2];
                sum += quals[--idx1];
            }
        }

        return result;
    }

    /**
     * Computes the mean of quality scores.
     * If the start and end are "out of bounds", they are silently corrected.
     * @param quals The quality scores.
     * @param start The starting index.
     * @param end The ending index.
     * @return The average.
     */
    public static int computeAverageQuality( byte[] quals, int start, int end )
    {
        if ( quals == null )
            throw new IllegalArgumentException("quals array may not be null");

        int nnn = quals.length;
        if ( start > nnn )
            start = nnn;
        if ( end > nnn )
            end = nnn;

        int result = 0;
        if ( start < end )
        {
            int totQ = 0;
            for ( int iii = start; iii < end; ++iii )
            {
                totQ += quals[iii];
            }
            result = totQ / (end - start);
        }
        return result;
    }

    private static final int JAFFE_SCORE = Integer.getInteger("QTJaffeScore",13).intValue();
    private static final int JAFFE_WINDOW = Integer.getInteger("QTJaffeWindow",12).intValue();
    private static final int READLENGTH_SCORE = Integer.getInteger("QTReadLengthScore",15).intValue();
    private static final int READLENGTH_WINDOW = Integer.getInteger("QTReadLengthWindow",20).intValue();
    private static final int AVG_QUALITY_START = Integer.getInteger("QTAvgQualityStart",100).intValue();
    private static final int AVG_QUALITY_END = Integer.getInteger("QTAvgQualityEnd",300).intValue();
    private static final int SEQPASS_AVG_QUALITY = Integer.getInteger("QTSeqPassAvgQuality",20).intValue();
    public static final int ANASNP_FLANKING_WINDOW = Integer.getInteger("QTAnaSNPFlankingWindow",5).intValue();
    public static final int ANASNP_WINDOW = Integer.getInteger("QTAnaSNPWindow",2*ANASNP_FLANKING_WINDOW + 1).intValue();
    public static final int ANASNP_FLANKING_SCORE = Integer.getInteger("QTAnaSNPFlankingScore",15).intValue();
    public static final int ANASNP_CENTRAL_SCORE = Integer.getInteger("QTAnaSNPCentralScore",20).intValue();
}
