package org.broadinstitute.hellbender.utils;

import htsjdk.tribble.util.ParsingUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A class to test the speed of different string.split implementations.
 * This test is disabled by default because it takes 3.5 minutes to complete and we already have the information it provides.
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

    private static final double MS_TO_NS = 1000000.0;
    private static final double S_TO_MS  = 1000.0;

    //=========================================================================

    private static void printTimingString(final String decorator,
                           final long time_ns,
                           final int numSplits) {

        final double time_ms =  time_ns / MS_TO_NS;

        final double timePerSplit_ns = ((double)time_ns) / ((double)numSplits) / ((double)numIterations);
        final double timePerSplit_ms = timePerSplit_ns / MS_TO_NS;

        System.out.println("\t" + decorator + " Total Time:\t" + time_ns + "ns\t" + time_ms +
                            "ms\tPer Split:\t" + timePerSplit_ns + "ns\t\t" + timePerSplit_ms + "ms");
    }

    private static void printTimingTable(final long   javaSplitWordTotalTime_ns,
                                         final long   javaSplitSingleCharStringTotalTime_ns,
                                         final long   htsjdkSplitSingleCharTotalTime_ns,
                                         final long   gatkSplitWordTotalTime_ns,
                                         final long   gatkSplitSingleCharStringTotalTime_ns,
                                         final long   overallElapsedTime_ns,
                                         final double overallElapsedTime_ms) {
        System.out.println("================================================================================");
        System.out.println("Timing Results:");
        System.out.println("--------------------------------------------------------------------------------");
        System.out.println("Overall Elapsed Time: " + overallElapsedTime_ns + "ns, " + overallElapsedTime_ms + "ms");
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

    private static void printTimingTableMarkdown(final long   javaSplitWordTotalTime_ns,
                                                 final long   javaSplitSingleCharStringTotalTime_ns,
                                                 final long   htsjdkSplitSingleCharTotalTime_ns,
                                                 final long   gatkSplitWordTotalTime_ns,
                                                 final long   gatkSplitSingleCharStringTotalTime_ns ) {

        System.out.println( "| Method | Benchmark | Total Time (ns) | Total Time (ms) | Time Per Split Operation (ns) | Time Per Split Operation (ms) |" );
        System.out.println( "| --- | --- | --- | --- | --- | --- |" );
        printTimingMarkdownLine( "| Java String::split | Split on Words | ", javaSplitWordTotalTime_ns, wordsToSplitOn.size() );
        printTimingMarkdownLine( "| Java String::split | Split on Chars | ", javaSplitSingleCharStringTotalTime_ns, singleCharStringsToSplitOn.size() );
        System.out.println( "| HTSJDK ParsingUtils::split | Split on Words | NA | NA | NA | NA |" );
        printTimingMarkdownLine( "| HTSJDK ParsingUtils::split | Split on Chars | ", htsjdkSplitSingleCharTotalTime_ns, singleCharStringsToSplitOn.size() );
        printTimingMarkdownLine( "| GATK Utils::split | Split on Words | ", gatkSplitWordTotalTime_ns, wordsToSplitOn.size() );
        printTimingMarkdownLine( "| GATK Utils::split | Split on Chars | ", gatkSplitSingleCharStringTotalTime_ns, singleCharStringsToSplitOn.size() );
    }

    private static void printTimingMarkdownLine(final String rowHeader,
                                                final long time_ns,
                                                final int numSplits) {

        final double time_ms =  time_ns / MS_TO_NS;

        final double timePerSplit_ns = ((double)time_ns) / ((double)numSplits) / ((double)numIterations);
        final double timePerSplit_ms = timePerSplit_ns / MS_TO_NS;

        System.out.println( rowHeader + time_ns + " | " + time_ms + " | " + timePerSplit_ns + " | " + timePerSplit_ms + " |" );
    }

    //==================================================================================================================

    // Disabled so that we don't waste time.
    // Cached results here:
//--------------------------------------------------------------------------------
//    Overall Elapsed Time: 194969742609ns, 194969.742609ms, 3.25min
//--------------------------------------------------------------------------------
//    Java Split String:
//    Words Total Time:	123581339505ns	123581.339505ms	Per Split:	5668.868784633028ns		0.0056688687846330275ms
//	"Chars" Total Time:	13324148242ns	13324.148242ms	Per Split:	3098.6391260465116ns		0.0030986391260465116ms
//    HTSJDK ParsingUtils::split
//	"Chars" Total Time:	6054845269ns	6054.845269ms	Per Split:	1408.1035509302328ns		0.0014081035509302328ms
//    GATK Utils::split
//    Words Total Time:	45913524687ns	45913.524687ms	Per Split:	2106.124985642202ns		0.002106124985642202ms
//	"Chars" Total Time:	6095292091ns	6095.292091ms	Per Split:	1417.509788604651ns		0.001417509788604651ms
//================================================================================
    @Test(enabled = false)
    void compareTimingForSplitString() {

        final long overallStartTime = System.nanoTime();

        //------------------------------------------------------------------------------------------------------------------
        // Baseline Java String.split:

        //-------------------------------------------------
        // First we do words:
        System.out.print("Testing Java String.split on words..."); System.out.flush();
        final long javaSplitWordStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : wordsToSplitOn ) {
                stringToSplit.split(word);
            }
        }
        final long javaSplitWordStopTime = System.nanoTime();
        final long javaSplitWordTotalTime_ns = javaSplitWordStopTime - javaSplitWordStartTime;
        System.out.println( " Done! (" + (javaSplitWordTotalTime_ns / MS_TO_NS / S_TO_MS) + "s)" );

        //-------------------------------------------------
        // Now we do single char words:
        System.out.print("Testing Java String.split on single characters..."); System.out.flush();
        final long javaSplitSingleCharStringStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : singleCharStringsToSplitOn ) {
                stringToSplit.split(word);
            }
        }
        final long javaSplitSingleCharStringStopTime = System.nanoTime();
        final long javaSplitSingleCharStringTotalTime_ns = javaSplitSingleCharStringStopTime - javaSplitSingleCharStringStartTime;
        System.out.println( " Done! (" + (javaSplitSingleCharStringTotalTime_ns / MS_TO_NS / S_TO_MS) + "s)" );


        //------------------------------------------------------------------------------------------------------------------
        // HTSJDK's ParsingUtils::split:

        //-------------------------------------------------
        // First we do words:

        // NOT APPLICABLE FOR HTSJDK's ParsingUtils::split!

        //-------------------------------------------------
        // Now we do single char words:
        System.out.print("Testing HTSJDK's ParsingUtils::split on single characters..."); System.out.flush();
        final long htsjdkSplitSingleCharStringStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : singleCharStringsToSplitOn ) {
                ParsingUtils.split( stringToSplit, word.charAt(0) );
            }
        }
        final long htsjdkSplitSingleCharStringStopTime = System.nanoTime();
        final long htsjdkSplitSingleCharTotalTime_ns = htsjdkSplitSingleCharStringStopTime - htsjdkSplitSingleCharStringStartTime;
        System.out.println( " Done! (" + (htsjdkSplitSingleCharTotalTime_ns / MS_TO_NS / S_TO_MS) + "s)" );

        //------------------------------------------------------------------------------------------------------------------
        // GATK's Utils::split:

        //-------------------------------------------------
        // First we do words:
        System.out.print("Testing GATK's Utils::split: on words..."); System.out.flush();
        final long gatkSplitWordStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : wordsToSplitOn ) {
                Utils.split( stringToSplit, word );
            }
        }
        final long gatkSplitWordStopTime = System.nanoTime();
        final long gatkSplitWordTotalTime_ns = gatkSplitWordStopTime - gatkSplitWordStartTime;
        System.out.println( " Done! (" + (gatkSplitWordTotalTime_ns / MS_TO_NS / S_TO_MS) + "s)" );

        //-------------------------------------------------
        // Now we do single char words:
        System.out.print("Testing GATK's Utils::split: on single characters..."); System.out.flush();
        final long gatkSplitSingleCharStringStartTime = System.nanoTime();

        for ( int i = 0; i < numIterations; ++i ) {
            for ( final String word : singleCharStringsToSplitOn ) {
                Utils.split( stringToSplit, word );
            }
        }
        final long gatkSplitSingleCharStringStopTime = System.nanoTime();
        final long gatkSplitSingleCharStringTotalTime_ns = gatkSplitSingleCharStringStopTime - gatkSplitSingleCharStringStartTime;
        System.out.println( " Done! (" + (gatkSplitSingleCharStringTotalTime_ns / MS_TO_NS / S_TO_MS) + "s)" );

        //------------------------------------------------------------------------------------------------------------------

        final long overallEndTime = System.nanoTime();
        final long overallElapsedTime_ns = overallEndTime - overallStartTime;
        final double overallElapsedTime_ms = overallElapsedTime_ns / MS_TO_NS;

        //------------------------------------------------------------------------------------------------------------------

        // Print results:
        printTimingTable(javaSplitWordTotalTime_ns, javaSplitSingleCharStringTotalTime_ns, htsjdkSplitSingleCharTotalTime_ns, gatkSplitWordTotalTime_ns, gatkSplitSingleCharStringTotalTime_ns, overallElapsedTime_ns, overallElapsedTime_ms);

        System.out.println("");
        System.out.println("====");
        System.out.println("");

        printTimingTableMarkdown(javaSplitWordTotalTime_ns, javaSplitSingleCharStringTotalTime_ns, htsjdkSplitSingleCharTotalTime_ns, gatkSplitWordTotalTime_ns, gatkSplitSingleCharStringTotalTime_ns );
    }

}
