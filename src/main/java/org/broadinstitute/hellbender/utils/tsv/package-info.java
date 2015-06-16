/**
 * Utility classes to read and write tab separated value (tsv) files.
 * <h3>File format description</h3>
 * <p>
 * Tab separated values may contain any number of <i>comment lines</i> (started with {@value org.broadinstitute.hellbender.utils.tsv.TableConstants#COMMENT_PREFIX}),
 * a column name containing line (aka. the <i>header line</i>) and any number of <i>data lines</i> one per record.
 * </p>
 * <p>While comment lines can contain any sequence of characters, the header and data lines are divided in
 * columns using exactly one {@value org.broadinstitute.hellbender.utils.tsv.TableConstants#COLUMN_SEPARATOR_STRING} character.</p>
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
 * <p>Values can be quoted using {@value org.broadinstitute.hellbender.utils.tsv.TableConstants#QUOTE_STRING}. This becomes necessary when the value contain
 * any special formatting characters like a new-line, the quote character itself, the column separator character or
 * the escape character {@value org.broadinstitute.hellbender.utils.tsv.TableConstants#ESCAPE_STRING}.</p>
 * <p>Within quotes, especial characters must be escaped using the {@value org.broadinstitute.hellbender.utils.tsv.TableConstants#ESCAPE_STRING}</p>
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
 * a top- or inner class and overriding {@link org.broadinstitute.hellbender.utils.tsv.TableReader#record(DataLine) record}
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
 *            protected void processTableColumns(final TableColumns columns) {
 *                  if (!columns.containsExactly("CONTIG","START","END"))
 *                      throw formatException("Bad column names");
 *            }
 *
 *            &#64;Override
 *            protected TableCounts record(final DataLine dataLine) {
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
 * a top- or inner class and overriding {@link org.broadinstitute.hellbender.utils.tsv.TableWriter#dataLine(java.lang.Object) dataLine}
 * method to map your record object type to a output line, represented by a {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine}.
 * <p>
 * Instances of {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine} can be obtained by calling {@link org.broadinstitute.hellbender.utils.tsv.DataLine DataLine}
 * can be obtained by calling the writers protected parameter-less method {@link org.broadinstitute.hellbender.utils.tsv.TableWriter#dataLine() dataLine()}.
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
 *            protected DataLine dataLine(final SimpleInterval interval) {
 *                // we can use append with confidence because we know the column order.
 *                return dataLine()
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
 */
package org.broadinstitute.hellbender.utils.tsv;