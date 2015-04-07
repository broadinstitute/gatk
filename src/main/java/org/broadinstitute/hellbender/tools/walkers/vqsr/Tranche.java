/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.text.XReadLines;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 10, 2011
 */

public class Tranche {
    private static final int CURRENT_VERSION = 5;

    public double ts, minVQSLod, knownTiTv, novelTiTv;
    public int numKnown,numNovel;
    public String name;
    public VariantRecalibratorArgumentCollection.Mode model;

    int accessibleTruthSites = 0;
    int callsAtTruthSites = 0;

    public Tranche(double ts, double minVQSLod, int numKnown, double knownTiTv, int numNovel, double novelTiTv, int accessibleTruthSites, int callsAtTruthSites, VariantRecalibratorArgumentCollection.Mode model) {
        this(ts, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, model, "anonymous");
    }

    public Tranche(double ts, double minVQSLod, int numKnown, double knownTiTv, int numNovel, double novelTiTv, int accessibleTruthSites, int callsAtTruthSites, VariantRecalibratorArgumentCollection.Mode model, String name ) {
        this.ts = ts;
        this.minVQSLod = minVQSLod;
        this.novelTiTv = novelTiTv;
        this.numNovel = numNovel;
        this.knownTiTv = knownTiTv;
        this.numKnown = numKnown;
        this.model = model;
        this.name = name;

        this.accessibleTruthSites = accessibleTruthSites;
        this.callsAtTruthSites = callsAtTruthSites;

        if ( ts < 0.0 || ts > 100.0)
            throw new UserException("Target FDR is unreasonable " + ts);

        if ( numKnown < 0 || numNovel < 0)
            throw new ReviewedGATKException("Invalid tranche - no. variants is < 0 : known " + numKnown + " novel " + numNovel);

        if ( name == null )
            throw new ReviewedGATKException("BUG -- name cannot be null");
    }

    private double getTruthSensitivity() {
        return accessibleTruthSites > 0 ? callsAtTruthSites / (1.0*accessibleTruthSites) : 0.0;
    }

    public static class TrancheTruthSensitivityComparator implements Comparator<Tranche>, Serializable {
        @Override
        public int compare(final Tranche tranche1, final Tranche tranche2) {
            return Double.compare(tranche1.ts, tranche2.ts);
        }
    }

    @Override
    public String toString() {
        return String.format("Tranche ts=%.2f minVQSLod=%.4f known=(%d @ %.4f) novel=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
                ts, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, name);
    }

    /**
     * Returns an appropriately formatted string representing the raw tranches file on disk.
     *
     * @param tranches
     * @return
     */
    public static String tranchesString( final List<Tranche> tranches ) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        final PrintStream stream = new PrintStream(bytes);

        if( tranches.size() > 1 )
            Collections.sort(tranches, new TrancheTruthSensitivityComparator());

        stream.println("# Variant quality score tranches file");
        stream.println("# Version number " + CURRENT_VERSION);
        stream.println("targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity");

        Tranche prev = null;
        for ( Tranche t : tranches ) {
            stream.printf("%.2f,%d,%d,%.4f,%.4f,%.4f,VQSRTranche%s%.2fto%.2f,%s,%d,%d,%.4f%n",
                    t.ts, t.numKnown, t.numNovel, t.knownTiTv, t.novelTiTv, t.minVQSLod, t.model.toString(),
                    (prev == null ? 0.0 : prev.ts), t.ts, t.model.toString(), t.accessibleTruthSites, t.callsAtTruthSites, t.getTruthSensitivity());
            prev = t;
        }

        return bytes.toString();
    }

    private static double getDouble(Map<String,String> bindings, String key, boolean required) {
        if ( bindings.containsKey(key) ) {
            String val = bindings.get(key);
            return Double.valueOf(val);
        }
        else if ( required ) {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
        else
            return -1;
    }

    private static int getInteger(Map<String,String> bindings, String key, boolean required) {
        if ( bindings.containsKey(key) )
            return Integer.valueOf(bindings.get(key));
        else if ( required ) {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
        else
            return -1;
    }

    /**
     * Returns a list of tranches, sorted from most to least specific, read in from file f
     *
     * @param f
     * @return
     */
    public static List<Tranche> readTranches(File f) {
        String[] header = null;
        List<Tranche> tranches = new ArrayList<Tranche>();

        try {
            for( final String line : new XReadLines(f) ) {
                if ( line.startsWith("#") )
                    continue;

                final String[] vals = line.split(",");
                if( header == null ) {
                    header = vals;
                    if ( header.length == 5 || header.length == 8 || header.length == 10 )
                        // old style tranches file, throw an error
                        throw new UserException.MalformedFile(f, "Unfortunately your tranches file is from a previous version of this tool and cannot be used with the latest code.  Please rerun VariantRecalibrator");
                    if ( header.length != 11 )
                        throw new UserException.MalformedFile(f, "Expected 11 elements in header line " + line);
                } else {
                    if ( header.length != vals.length )
                        throw new UserException.MalformedFile(f, "Line had too few/many fields.  Header = " + header.length + " vals " + vals.length + ". The line was: " + line);

                    Map<String,String> bindings = new HashMap<String, String>();
                    for ( int i = 0; i < vals.length; i++ ) bindings.put(header[i], vals[i]);
                    tranches.add(new Tranche(getDouble(bindings,"targetTruthSensitivity", true),
                            getDouble(bindings,"minVQSLod", true),
                            getInteger(bindings,"numKnown", false),
                            getDouble(bindings,"knownTiTv", false),
                            getInteger(bindings,"numNovel", true),
                            getDouble(bindings,"novelTiTv", true),
                            getInteger(bindings,"accessibleTruthSites", false),
                            getInteger(bindings,"callsAtTruthSites", false),
                            VariantRecalibratorArgumentCollection.parseString(bindings.get("model")),
                            bindings.get("filterName")));
                }
            }

            Collections.sort(tranches, new TrancheTruthSensitivityComparator());
            return tranches;
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(f, e);
        }
    }
}
