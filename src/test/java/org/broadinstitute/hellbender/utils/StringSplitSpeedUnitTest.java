package org.broadinstitute.hellbender.utils;

import htsjdk.tribble.util.ParsingUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A class to test the speed of different string.split implementations.
 * Created by jonn on 11/1/17.
 */
public class StringSplitSpeedUnitTest extends GATKBaseTest {

    private static final String stringToSplit =
            "Arma virumque cano, Troiae qui primus ab oris " +
            "Italiam, fato profugus, Laviniaque venit " +
            "litora, multum ille et terris iactatus et alto " +
            "vi superum saevae memorem Iunonis ob iram; " +
            "multa quoque et bello passus, dum conderet urbem, " +
            "inferretque deos Latio, genus unde Latinum, " +
            "Albanique patres, atque altae moenia Romae.  " +
            "Musa, mihi causas memora, quo numine laeso, " +
            "quidve dolens, regina deum tot volvere casus " +
            "insignem pietate virum, tot adire labores " +
            "impulerit. Tantaene animis caelestibus irae?  " +
            "Urbs antiqua fuit, Tyrii tenuere coloni, " +
            "Karthago, Italiam contra Tiberinaque longe " +
            "ostia, dives opum studiisque asperrima belli; " +
            "quam Iuno fertur terris magis omnibus unam " +
            "posthabita coluisse Samo; hic illius arma, " +
            "hic currus fuit; hoc regnum dea gentibus esse, " +
            "si qua fata sinant, iam tum tenditque fovetque. " +
            "Progeniem sed enim Troiano a sanguine duci " +
            "audierat, Tyrias olim quae verteret arces; " +
            "hinc populum late regem belloque superbum " +
            "venturum excidio Libyae: sic volvere Parcas.  " +
            "Id metuens, veterisque memor Saturnia belli, " +
            "prima quod ad Troiam pro caris gesserat Argis- " +
            "necdum etiam causae irarum saevique dolores " +
            "exciderant animo: manet alta mente repostum " +
            "iudicium Paridis spretaeque iniuria formae, " +
            "et genus invisum, et rapti Ganymedis honores.  " +
            "His accensa super, iactatos aequore toto " +
            "Troas, reliquias Danaum atque immitis Achilli,  " +
            "arcebat longe Latio, multosque per annos " +
            "errabant, acti fatis, maria omnia circum.  " +
            "Tantae molis erat Romanam condere gentem!";

    private static final List<String> wordsToSplitOn = new ArrayList<>( Arrays.asList(stringToSplit.split(" ")) );

    private static final Set<String> singleCharStringsToSplitOn =
                        Arrays.stream(stringToSplit.split(""))
                                .map(s -> (s.equals("?")) ? "\\?" : s)
                                .collect(Collectors.toSet());

    private static final int numIterations = 100000;

    //=========================================================================

    void printTimingString(final String decorator,
                           final long time_ns,
                           final int numSplits) {

        final double time_ms =  time_ns / 1000000.0;

        final double timePerSplit_ns = ((double)time_ns) / ((double)numSplits) / ((double)numIterations);
        final double timePerSplit_ms = timePerSplit_ns / 1000000.0;

        System.out.println("\t" + decorator + " Total Time:\t" + time_ns + "ns\t" + time_ms +
                            "ms\tPer Split:\t" + timePerSplit_ns + "ns\t\t" + timePerSplit_ms + "ms");
    }

    //==================================================================================================================

    // Disabled so that we don't waste time.
    // Cached results here:
//    --------------------------------------------------------------------------------
//    Timing Results:
//    --------------------------------------------------------------------------------
//    Java Split String:
//          Words Total Time:	131867865203ns	131867.865203ms	Per Split:	6048.98464233945ns		0.00604898464233945ms
//	        "Chars" Total Time:	12917243085ns	12917.243085ms	Per Split:	3004.010019767442ns		0.003004010019767442ms
//    HTSJDK ParsingUtils::split
//	        "Chars" Total Time:	5882790859ns	5882.790859ms	Per Split:	1368.0908974418605ns	0.0013680908974418606ms
//    GATK Utils::split
//          Words Total Time:	38734463275ns	38734.463275ms	Per Split:	1776.8102419724771ns	0.0017768102419724772ms
//	        "Chars" Total Time:	7120052467ns	7120.052467ms	Per Split:	1655.826155116279ns		0.0016558261551162792ms
//================================================================================
    @Test(enabled = false)
    void compareTimingForSplitString() {

//    Java String.split,
//    HTSJDK's ParsingUtils::split
//    GATK3's Utils::split.

        //------------------------------------------------------------------------------------------------------------------
        // Baseline Java String.split:

        //-------------------------------------------------
        // First we do words:
        final long javaSplitWordStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : wordsToSplitOn ) {
                stringToSplit.split(word);
            }
        }
        final long javaSplitWordStopTime = System.nanoTime();
        final long javaSplitWordTotalTime_ns = javaSplitWordStopTime - javaSplitWordStartTime;

        //-------------------------------------------------
        // Now we do single char words:
        final long javaSplitSingleCharStringStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : singleCharStringsToSplitOn ) {
                stringToSplit.split(word);
            }
        }
        final long javaSplitSingleCharStringStopTime = System.nanoTime();
        final long javaSplitSingleCharStringTotalTime_ns = javaSplitSingleCharStringStopTime - javaSplitSingleCharStringStartTime;

        //------------------------------------------------------------------------------------------------------------------
        // HTSJDK's ParsingUtils::split:

        //-------------------------------------------------
        // First we do words:

        // NOT APPLICABLE FOR HTSJDK's ParsingUtils::split!

        //-------------------------------------------------
        // Now we do single char words:
        final long htsjdkSplitSingleCharStringStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : singleCharStringsToSplitOn ) {
                ParsingUtils.split( stringToSplit, word.charAt(0) );
            }
        }
        final long htsjdkSplitSingleCharStringStopTime = System.nanoTime();
        final long htsjdkSplitSingleCharTotalTime_ns = htsjdkSplitSingleCharStringStopTime - htsjdkSplitSingleCharStringStartTime;

        //------------------------------------------------------------------------------------------------------------------
        // GATK's Utils::split:

        //-------------------------------------------------
        // First we do words:
        final long gatkSplitWordStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : wordsToSplitOn ) {
                Utils.split( stringToSplit, word );
            }
        }
        final long gatkSplitWordStopTime = System.nanoTime();
        final long gatkSplitWordTotalTime_ns = gatkSplitWordStopTime - gatkSplitWordStartTime;

        //-------------------------------------------------
        // Now we do single char words:
        final long gatkSplitSingleCharStringStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : singleCharStringsToSplitOn ) {
                Utils.split( stringToSplit, word );
            }
        }
        final long gatkSplitSingleCharStringStopTime = System.nanoTime();
        final long gatkSplitSingleCharStringTotalTime_ns = gatkSplitSingleCharStringStopTime - gatkSplitSingleCharStringStartTime;

        //------------------------------------------------------------------------------------------------------------------

        // Print results:
        System.out.println("================================================================================");
        System.out.println("Timing Results:");
        System.out.println("--------------------------------------------------------------------------------");
        System.out.println("Java Split String:");
        printTimingString("Words", javaSplitWordTotalTime_ns, wordsToSplitOn.size());
        printTimingString("\"Chars\"", javaSplitSingleCharStringTotalTime_ns, singleCharStringsToSplitOn.size());
        System.out.println("HTSJDK ParsingUtils::split");
        printTimingString("\"Chars\"", htsjdkSplitSingleCharTotalTime_ns, singleCharStringsToSplitOn.size());
        System.out.println("GATK Utils::split");
        printTimingString("Words", gatkSplitWordTotalTime_ns, wordsToSplitOn.size());
        printTimingString("\"Chars\"", gatkSplitSingleCharStringTotalTime_ns, singleCharStringsToSplitOn.size());
        System.out.println("================================================================================");
        
    }

}
