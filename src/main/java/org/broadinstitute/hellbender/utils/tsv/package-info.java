/**
 * Utility classes to read and write tab separated value (tsv) files.
 * <h3>File format description<a name="format"></a></h3>
 *
 * <p>
 * Tab separated values may contain any number of <i>comment lines</i> (started with {@value org.broadinstitute.hellbender.utils.tsv.TableUtils#COMMENT_PREFIX}),
 * a column name containing line (aka. the <i>header line</i>) and any number of <i>data lines</i> one per record.
 * </p>
 * <p>While comment lines can contain any sequence of characters, the header and data lines are divided in
 * columns using exactly one {@value org.broadinstitute.hellbender.utils.tsv.TableUtils#COLUMN_SEPARATOR_STRING} character.</p>
 * <p>Blank lines are treated as having a single column with the empty string as the only value (or column name)</p>
 * <p>
 * The header line is the first non-comment line, whereas any other non-comment line after that is
 * considered a data line. Comment lines can appear anywhere in the file and their
 * present is ignored by the reader ({@link org.broadinstitute.hellbender.utils.tsv.TableReader TableReader} implementations).
 * </p>
 * <p>
 * The header line values, the column names, must all be different (otherwise a formatting exception will be thrown), and
 * all data lines have to have as many values as there are columns in the header line.
 * </p>
 * <p>Values can be quoted using {@value org.broadinstitute.hellbender.utils.tsv.TableUtils#QUOTE_STRING}. This becomes necessary when the value contain
 * any special formatting characters like a new-line, the quote character itself, the column separator character or
 * the escape character {@value org.broadinstitute.hellbender.utils.tsv.TableUtils#ESCAPE_STRING}.</p>
 * <p>Within quotes, especial characters must be escaped using the {@value org.broadinstitute.hellbender.utils.tsv.TableUtils#ESCAPE_STRING}</p>
 * <p>Examples 1:</p>
 * <pre>
 *     # comment 1
 *     # comment 2
 *     CONTIG   START   END     NAME    SAMPLE1 SAMPLE2
 *     # comment 3
 *     chr1     123100  123134 tgt_0    100.0   102.0
 *     chr1     134012  134201 tgt_1    50      12
 *     # comment 4
 *     chr2     ...
 * </pre>
 * <h3>Reading tsv files</h3>
 * You will need to extend class
 * {@link org.broadinstitute.hellbender.utils.tsv.TableReader TableReader}, either using
 * a top- or inner class and overriding {@link org.broadinstitute.hellbender.utils.tsv.TableReader#createRecord(DataLine) createRecord}
 * method to map input data-lines, wrapped into a {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine}, to
 * your row element class of choice.
 * <p>
 * Example, a SimpleInterval reader from a tsv file with three columns, CONTIG, START and END:
 * </p>
 * <pre>
 *
 *     ...
 *
 *     public void doWork(final File inputFile) throws IOException {
 *
 *         final TableReader&lt;SimpleInterval&gt; reader = new TableReader&lt;SimpleInterval&gt(inputFile) {
 *
 *            // Optional (but recommended) check that the columns in the file are the ones expected:
 *            &#64;Override
 *            protected void processColumns(final TableColumns columns) {
 *                  if (!columns.containsExactly("CONTIG","START","END"))
 *                      throw formatException("Bad column names");
 *            }
 *
 *            &#64;Override
 *            protected TableCounts createRecord(final DataLine dataLine) {
 *                return new SimpleInterval(dataLine.get("CONTIG"),
 *                                       dataLine.getInt("START"),
 *                                       dataLine.getInt("END"));
 *            }
 *         };
 *
 *         for (final SimpleInterval interval : reader) {
 *             // whatever you wanna do per interval.
 *         }
 *         reader.close();
 *         ...
 *
 *     }
 * </pre>
 * <h3>Writing tsv files</h3>
 * You will need to extend class
 * {@link org.broadinstitute.hellbender.utils.tsv.TableWriter TableWriter}, either using
 * a top- or inner class and overriding {@link org.broadinstitute.hellbender.utils.tsv.TableWriter#composeLine composeLine}
 * method to map your record object type to a output line, represented by a {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine}.
 * <p>
 * Instances of {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine} can be obtained by calling {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine}
 * can be obtained by calling the writers protected parameter-less method {@link org.broadinstitute.hellbender.utils.tsv.TableWriter#composeLine composeLine}.
 * </p>
 * <p>
 * The column names are passed in order to the constructor.
 * </p>
 * <p>
 * Example:
 * </p>
 * <pre>
 *     public void doWork(final File outputFile) throws IOException {
 *
 *         final TableWriter&lt;SimpleInterval&gt; writer =
 *              new TableWriter&lt;SimpleInterval&gt(outputFile, new TableColumns("CONTIG","START","END")) {
 *            &#64;Override
 *            protected void composeLine(final SimpleInterval interval, final DataLine dataLine) {
 *                // we can use append with confidence because we know the column order.
 *                dataLine
 *                    .append(interval.getContig())
 *                    .append(interval.getStart(),interval.getEnd());
 *            }
 *         };
 *
 *         for (final SimpleInterval interval : intervalsToWrite) {
 *             writer.writeRecord(interval);
 *         }
 *         writer.close();
 *         ...
 *
 *     }
 * </pre>
 * <h3>Readers and Writers using function composition</h3>
 * {@link org.broadinstitute.hellbender.utils.tsv.TableUtils TableUtils} contains methods to create
 * readers and writers without the need to explicitly extending {@link org.broadinstitute.hellbender.utils.tsv.TableReader TableReader}
 * or {@link org.broadinstitute.hellbender.utils.tsv.TableWriter TableWriter} but by specifying their behaviour through
 * lambda functions.
 * <p>Example of a reader:
 * <pre>
 *     final TableReader&lt;SimpleInterval&gt; reader = TableUtils.reader(inputFile,
 *                (columns,formatExceptionFactory) -> {
 *                   // we check the columns is what we except them to be:
 *                   if (!columns.matchesExactly("CONTIG","START","END"))
 *                      throw formatExceptionFactory.apply("Bad header");
 *                   // we return the lambda to translate dataLines into intervals.
 *                   return (dataLine) -> new SimpleIntervals(dataLine.get(0),dataLine.getInt(1),dataLine.getInt(2));
 *                });
 * </pre>
 * <p>
 * The lambda that you need to indicate seems a bit complicate but is not so... basically it receives the
 * columns in the input and it must return another lambda that will translate data-lines into records considering
 * those columns.
 * </p>
 * <p>
 * Before doing that, it check whether the columns are the excepted ones and int the correct order (always recommended).
 * </p><p>
 *     The additional formatExceptionFactory parameter allows the reader implementation to correctly report formatting issues.
 * </p></p>
 * <p>Example of a writer:</p>
 * <pre>
 *     final TableWriter&lt;SimpleInterval&gt; reader = TableUtils.reader(outputFile,
 *                new TableColumnCollection("CONTIG","START","END"),
 *                (interval,dataLine) -> {
 *                  dataLine.append(interval.getContig()
 *                          .append(interval.getStart(),interval.getEnd());
 *                });
 * </pre>
 * <p>
 *     The case of the writer is far more simple as there is no need to report formatting errors as we are
 *     the ones producing the file.
 * </p>
 */
package org.broadinstitute.hellbender.utils.tsv;