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

package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Option;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.tools.recalibration.BaseRecalibration;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Optional;


/**
 * Tool to analyze and evaluate base recalibration tables.
 * <p/>
 * It generates plots to assess the quality of a recalibration run.
 *
 * <h3>Input</h3>
 *
 * The tool can take up to three different sets of recalibration tables.
 * The resulting plots will be overlaid on top of each other to make
 * comparisons easy.
 *
 * <table style="text-align: left">
 *     <thead>
 *       <tr><th>Set</th><th>Argument</th><th>Label</th><th>Color</th><th>Description</th></tr>
 *     </thead>
 *     <tbody>
 *       <tr><td>Original</td><td>-before</td><td>BEFORE</td><td style="color: #ff34b3">Maroon1</td>
 *         <td>First pass recalibration
 *             tables obtained from applying {@link BaseRecalibration}
 *             on the original alignment.</td></tr>
 *       <tr><td>Recalibrated</td><td>-after</td><td>AFTER</td><td style="color: #0000ff">Blue</td>
 *         <td>Second pass recalibration tables
 *             results from the application of {@link BaseRecalibration}
 *             on the alignment recalibrated using the first pass tables</td></tr>
 *       <tr><td>Input</td><td>-BQSR</td><td>BQSR</td><td style="color: #000000">Black</td>
 *           <td>Any recalibration table without a specific role</td></tr>
 *     </tbody>
 * </table>
 * <br/>
 *
 * You need to specify one set at least. Multiple sets need to have the same values for the following parameters:
 * <br/></br>
 * <i>covariate (order is not important), no_standard_covs, run_without_dbsnp, solid_recal_mode,
 * solid_nocall_strategy, mismatches_context_size, mismatches_default_quality, deletions_default_quality,
 * insertions_default_quality, maximum_cycle_value, low_quality_tail, default_platform, force_platform,
 * quantizing_levels</i> and <i>binary_tag_name</i>
 * <h3>Output</h3>
 *
 * Currently this tool generates two outputs:
 *
 * <dl>
 *   <dt style="font-weight: normal">-plots <i>my-report.pdf</i></dt>
 *   <dd>A pdf document that encloses plots to assess the quality of the recalibration.</dd>
 *   <dt style="font-weight: normal">-csv <i>my-report.csv</i></dt>
 *   <dd>A csv file that contains a table with all the data required to generate those plots.</dd>
 * </dl>
 *
 * You need to specify at least one of them.
 *
 * <h3>Other Arguments</h3>
 *
 * <h4>-ignoreLMT, --ignoreLastModificationTimes</h4>
 *
 * when set, no warning message will be displayed in the -before recalibration table file is older than the -after one.
 *
 * <h3>Examples</h3>
 *
 *
 * <h4>Plot a single recalibration table</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *      -T AnalyzeCovariates \
 *      -R myrefernce.fasta \
 *      -BQSR myrecal.table \
 *      -plots BQSR.pdf
 * </pre>
 *
 * <h4>Plot before (first pass) and after (second pass) recalibration table to compare them</h4>
 *
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *      -T AnalyzeCovariates \
 *      -R myrefernce.fasta \
 *      -before recal2.table \
 *      -after recal3.table \
 *      -plots recalQC.pdf
 * </pre>
 *
 * <h4>Plot up to three recalibration tables for comparison</h4>
 *
 * <pre>
 *
 * # You can ignore the before/after semantics completely if you like (if you do add -ignoreLMT
 * # to avoid a possible warning), but all tables should have been generated using the same parameters.
 *
 * java -jar GenomeAnalysisTK.jar \
 *      -T AnalyzeCovariates \
 *      -R myrefernce.fasta \
 *      -ignoreLMT \
 *      -BQSR recal1.table \   # you can discard any two
 *      -before recal2.table \
 *      -after recal3.table \
 *      -plots myrecals.pdf
 * </pre>
 *
 * <h4>Full BQSR quality assessment pipeline</h4>
 *
 * <pre>
 * # Generate the first pass recalibration table file.
 * java -jar GenomeAnalysisTK.jar \
 *      -T BaseRecalibrator \
 *      -R myreference.fasta \
 *      -I myinput.bam \
 *      -knownSites bundle/my-trusted-snps.vcf \ # optional but recommendable
 *      -knownSites bundle/my-trusted-indels.vcf \ # optional but recommendable
 *      ... other options
 *      -o firstpass.table
 *
 * # Generate the second pass recalibration table file.
 * java -jar GenomeAnalysisTK.jar \
 *      -T BaseRecalibrator \
 *      -BQSR firstpass.table \
 *      -R myreference.fasta \
 *      -I myinput.bam \
 *      -knownSites bundle/my-trusted-snps.vcf \
 *      -knownSites bundle/my-trusted-indels.vcf \
 *      ... other options \
 *      -o secondpass.table
 *
 * # Finally generate the plots and also keep a copy of the csv (optional).
 * java -jar GenomeAnalysisTK.jar \
 *      -T AnalyzeCovariates \
 *      -R myrefernce.fasta \
 *      -before firstpass.table \
 *      -after secondpass.table \
 *      -csv BQSR.csv \ # optional
 *      -plots BQSR.pdf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        usage = "TODO",
        usageShort = "Tool to analyze and evaluate base recalibration tables.",
        programGroup = QCProgramGroup.class
)
public final class AnalyzeCovariates extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(AnalyzeCovariates.class);

    // Constants on option short names that are used in some error/warning messages:

    static final String CSV_ARG_SHORT_NAME = "csv";
    static final String PDF_ARG_SHORT_NAME = "plots";
    static final String BEFORE_ARG_SHORT_NAME = "before";
    static final String AFTER_ARG_SHORT_NAME = "after";

    /**
     * File containing the recalibration tables from the first pass.
     */
    @Option(shortName=BEFORE_ARG_SHORT_NAME,fullName="beforeReportFile", doc = "file containing the BQSR first-pass report file",optional = true)
    protected File beforeFile = null;

    /**
     * File containing the recalibration tables from the second pass.
     */
    @Option(shortName=AFTER_ARG_SHORT_NAME, fullName="afterReportFile", doc = "file containing the BQSR second-pass report file",optional = true)
    protected File afterFile = null;

    /**
     * If true, it won't show a warning if the last-modification time of the before and after input files suggest that they have been reversed.
     */
    @Option(shortName="ignoreLMT", fullName="ignoreLastModificationTimes", doc= "do not emit warning messages related to suspicious last modification time order of inputs",optional = true)
    protected boolean ignoreLastModificationTime = false;

    /**
     * Output report file name.
     */
    @Option(shortName=PDF_ARG_SHORT_NAME, fullName="plotsReportFile" ,doc = "location of the output report",optional = true)
    protected File pdfFile = null;

    /**
     * Output csv file name.
     */
    @Option(shortName=CSV_ARG_SHORT_NAME,fullName="intermediateCsvFile" ,doc = "location of the csv intermediate file",optional = true)
    protected File csvFile = null;

    /**
     * Enables on-the-fly recalibrate of base qualities, intended primarily for use with BaseRecalibrator and PrintReads
     * (see Best Practices workflow documentation). The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Option(fullName="BQSR", shortName="BQSR", optional=true, doc="Input covariates table file for on-the-fly base quality score recalibration")
    public File BQSR_RECAL_FILE = null;

    /**
     * Checks inputs and argument values.
     * <p/>
     * Notice that this routine will not validate the content of files. It may have some minor side effects as
     * the output of warning messages back to the user.
     *
     * @throw IllegalStateException there is some required argument value that has not been loaded yet.
     * @throw UserException if there is some error caused by or under the end user's control.
     */
    private void checkArgumentsValues() {
        checkInputReportFile("BQSR",BQSR_RECAL_FILE);
        checkInputReportFile("before",beforeFile);
        checkInputReportFile("after",afterFile);
        if (BQSR_RECAL_FILE == null && beforeFile == null && afterFile == null) {
            throw new UserException("you must provide at least one recalibration report file "
                    + "(arguments -BQSR, -" + BEFORE_ARG_SHORT_NAME + " or -" + AFTER_ARG_SHORT_NAME);
        }

        checkOutputFile(PDF_ARG_SHORT_NAME,pdfFile);
        checkOutputFile(CSV_ARG_SHORT_NAME, csvFile);
        checkInputReportFileLMT(beforeFile,afterFile);
        checkOutputRequested();
    }

    /**
     * Checks whether the last-modification-time of the inputs is consistent with their relative roles.
     *
     * This routine does not thrown an exception but may output a warning message if inconsistencies are spotted.
     *
     * @param beforeFile the before report file.
     * @param afterFile  the after report file.
     */
    private void checkInputReportFileLMT(final File beforeFile, final File afterFile) {

        if (ignoreLastModificationTime  || beforeFile == null || afterFile == null) {
            return; // nothing to do here
        } else if (beforeFile.lastModified() > afterFile.lastModified()) {
            Utils.warnUser("Last modification timestamp for 'Before' and 'After'"
                    + "recalibration reports are in the wrong order. Perhaps, have they been swapped?");
        }
    }

    /**
     * Checks that at least one output was requested.
     *
     * @throw UserException if no output was requested.
     */
    private void checkOutputRequested() {
        if (pdfFile == null && csvFile == null) {
            throw new UserException("you need to request at least one output:"
                    + " the intermediate csv file (-" + CSV_ARG_SHORT_NAME + " FILE)"
                    + " or the final plot file (-" + PDF_ARG_SHORT_NAME + " FILE).");
        }
    }

    /**
     * Checks the value provided to input file arguments.
     *
     * @throw UserException if there is any problem cause by or under the end user's control
     *
     * @param name command line argument short name.
     * @param value the argument value.
     */
    private void checkInputReportFile(final String name,final File value) {
        if (value == null) {
            return;
        } else if (!value.exists()) {
            throw new UserException.BadArgumentValue(name, "input report '" +
                    value + "' does not exist or is unreachable");
        } else if (!value.isFile()) {
            throw new UserException.BadArgumentValue(name, "input report '" +
                    value + "' is not a regular file");
        } else if (!value.canRead()) {
            throw new UserException.BadArgumentValue(name, "input report '" +
                    value + "' cannot be read");
        }
    }

    /**
     * Checks the value provided for output arguments.
     *
     * @throw UserException if there is any problem cause by or under the end user's control
     *
     * @param name command line argument short name.
     * @param value the argument value.
     */
    private void checkOutputFile(final String name, final File value) {
        if (value == null) {
            return;
        }
        if (value.exists() && !value.isFile()) {
            throw new UserException.BadArgumentValue(name, "the output file location '"
                    + value + "' exists as not a file");
        }
        final File parent = value.getParentFile();
        if (parent == null) {
            return;
        }
        if (!parent.exists()) {
            throw new UserException.BadArgumentValue(name, "the output file parent directory '"
                    + parent + "' does not exists or is unreachable");
        } else if (!parent.isDirectory()) {
            throw new UserException.BadArgumentValue(name, "the output file parent directory '"
                    + parent + "' is not a directory");
        } else if (!parent.canWrite()) {
            throw new UserException.BadArgumentValue(name, "the output file parent directory '"
                    + parent + "' cannot be written");
        }

    }

    /**
     * Generates the plots using the external R script.
     *
     * <p/>
     * If <code>plotsFile</code> is <code>null</code>, it does not perform any plotting.
     *
     * @param csvFile the intermediary csv file.
     * @param plotsFile the output plot location.
     */
    private void generatePlots(final File csvFile, final Map<String,File> reportFiles, final File plotsFile) {

        if (plotsFile == null) {
            return;
        }
        logger.info("Generating plots file '" + plotsFile + "'");
        final File exampleReportFile = reportFiles.values().iterator().next();
        RecalUtils.generatePlots(csvFile,exampleReportFile,plotsFile);
    }

    @Override
    public Object doWork() {
        checkArgumentsValues();
        final Map<String, File> reportFiles = buildReportFileMap();
        final Map<String, RecalibrationReport> reports = buildReportMap(reportFiles);
        checkReportConsistency(reports);
        final File csvFile = resolveCsvFile();
        generateCsvFile(csvFile,reports);
        final File plotFile = resolvePlotFile();
        generatePlots(csvFile, reportFiles, plotFile);
        return Optional.empty();
    }

    /**
     * Returns the plot output file
     * @return might be <code>null</code> if the user has not indicated and output file.
     */
    private File resolvePlotFile() {
        return pdfFile;
    }

    /**
     * Generates the intermediary Csv file.
     *
     * @param csvFile where to write the file.
     * @param reports the reports to be included.
     */
    private void generateCsvFile(final File csvFile, final Map<String, RecalibrationReport> reports) {
        try {
            logger.info("Generating csv file '" + csvFile + "'");
            RecalUtils.generateCsv(csvFile, reports);
        } catch (FileNotFoundException e) {
            throw new UserException(
                    String.format("There is a problem creating the intermediary Csv file '%s': %s",
                            csvFile, e.getMessage()),e);
        }
    }

    /**
     * Checks whether multiple input recalibration report files argument values are consistent (equal).
     *
     * @param reports map with report to verify.
     *
     * @throw UserException if there is any inconsistency.
     */
    private void checkReportConsistency(final Map<String, RecalibrationReport> reports) {
        @SuppressWarnings({"unchecked", "rawtypes"})
        final Map.Entry<String,RecalibrationReport>[] reportEntries =
                reports.entrySet().toArray((Map.Entry<String,RecalibrationReport>[]) new Map.Entry[reports.size()]);

        final Map.Entry<String,RecalibrationReport> exampleEntry = reportEntries[0];

        for (int i = 1; i < reportEntries.length; i++) {
            final Map<String,? extends CharSequence> diffs = exampleEntry.getValue().getRAC().compareReportArguments(
                    reportEntries[i].getValue().getRAC(),exampleEntry.getKey(),reportEntries[i].getKey());
            if (diffs.size() != 0) {
                throw new UserException.IncompatibleRecalibrationTableParameters("There are differences in relevant arguments of"
                        + " two or more input recalibration reports. Please make sure"
                        + " they have been created using the same recalibration parameters."
                        + " " + Utils.join("// ", reportDifferencesStringArray(diffs)));
            }
        }
    }


    /**
     * Creates a map with all input recalibration files indexed by their "role".
     * <p/>
     * The key is the role and the value the corresponding report file.
     * <p/>
     * Roles: "Before" (recalibration), "After" (recalibration), "BQSR" (the tool standard argument recalibration file)
     *
     * @return never <code>null</code>
     */
    private Map<String, File> buildReportFileMap() {
        final Map<String,File> reports = new LinkedHashMap<>(3);
        if (BQSR_RECAL_FILE != null) {
            reports.put("BQSR",BQSR_RECAL_FILE);
        }
        if (beforeFile != null) {
            reports.put("Before",beforeFile);
        }
        if (afterFile != null) {
            reports.put("After",afterFile);
        }
        return reports;
    }

    /**
     * Transforms a recalibration file map into a report object map.
     *
     * @param reportFileMap the file map to transforms.
     * @return never <code>null</code>, a new map with the same size as
     *  <code>reportFileMap</code> and the same key set.
     */
    private Map<String, RecalibrationReport> buildReportMap(final Map<String, File> reportFileMap) {
        final Map<String,RecalibrationReport> reports = new LinkedHashMap<>(reportFileMap.size());
        for (final Map.Entry<String,File> e : reportFileMap.entrySet()) {
            reports.put(e.getKey(),new RecalibrationReport(e.getValue()));
        }
        return reports;
    }

    /**
     * Generates a flatter String array representation of recalibration argument differences.
     * @param diffs the differences to represent.
     *
     * @return never <code>null</code>, an array of the same length as the size of the input <code>diffs</code>.
     */
    private String[] reportDifferencesStringArray(final Map<String, ? extends CharSequence> diffs) {
        final String[] result = new String[diffs.size()];
        int i = 0;
        for (final Map.Entry<String, ? extends CharSequence> e : diffs.entrySet()) {
            result[i++] = capitalize(e.getKey()) + ": " + e.getValue();
        }
        return result;
    }

    /**
     * Returns the input string capitalizing the first letter.
     *
     * @param str the string to capitalize
     * @return never <code>null</code>.
     */
    private String capitalize(final String str) {
        if (str.isEmpty()) {
            return str;
        } else {
            return Character.toUpperCase(str.charAt(0)) + str.substring(1);
        }
    }

    /**
     * Returns the csv file to use.
     * <p/>
     * This is the the one specified by the user if any or a temporary file
     * that will be deleted as soon as the VM exists by default.
     *
     * @return never <code>null</code>.
     */
    private File resolveCsvFile() {
        if (csvFile != null)  {
            return csvFile;
        } else {
            try {
              final File result = File.createTempFile("AnalyzeCovariates", ".csv");
              result.deleteOnExit();
              return result;
            } catch (IOException e) {
                throw new UserException("Could not create temporary Csv file",e);
            }
        }
    }

}


