package main

import (
	"bufio"
	"context"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	"cloud.google.com/go/bigquery"
	"github.com/carbocation/genomisc"
	"google.golang.org/api/iterator"
)

const (
	// pheno file field columns
	eid = iota
)

type SamplePheno struct {
	Column
	SampleID string `db:"sample_id"`
	Value    string `db:"value"`
}

type Column struct {
	FieldID      string             `db:"FieldID"`
	Instance     string             `db:"instance"`
	ArrayIDX     string             `db:"array_idx"`
	CodingFileID bigquery.NullInt64 `db:"coding_file_id"`
}

type WrappedBigQuery struct {
	Context  context.Context
	Client   *bigquery.Client
	Project  string
	Database string
}

var (
	BufferSize = 4096 * 8
	STDOUT     = bufio.NewWriterSize(os.Stdout, BufferSize)
)

func main() {
	defer STDOUT.Flush()

	var (
		phenoPaths  flagSlice
		acknowledge bool
		BQ          = &WrappedBigQuery{}
	)
	flag.StringVar(&BQ.Project, "project", "broad-ml4cvd", "Name of the Google Cloud project that hosts your BigQuery database instance")
	flag.StringVar(&BQ.Database, "bigquery", "", "BigQuery phenotype database name")
	flag.Var(&phenoPaths, "pheno", "phenotype file for the UKBB. Pass this flag once per file if you want to process multiple files at once. Every FieldID that is seen in an earlier file will be ignored in later files.")
	flag.BoolVar(&acknowledge, "ack", false, "Acknowledge the limitations of the tool")
	flag.Parse()

	if phenoPaths.String() == "" || BQ.Database == "" {
		flag.PrintDefaults()
		os.Exit(1)
	}

	for _, phenoPath := range phenoPaths {
		if strings.HasSuffix(phenoPath, ".gz") {
			log.Printf("\n**\n**\nWARNING This tool does not currently operate on gzipped files. Based on your filename, this will likely crash. Please gunzip %s\n**\n**\n", phenoPath)
		}

		if _, err := os.Stat(phenoPath); os.IsNotExist(err) {
			log.Fatalf("Fatal error: %v does not exist\n", phenoPath)
		} else if err != nil {
			log.Fatalf("Fatal error: %v (possibly disk or permissions issues?): %v\n", phenoPath, err)
		}
	}

	if !acknowledge {
		fmt.Fprintln(os.Stderr, "!! -- !!")
		fmt.Fprintln(os.Stderr, "NOTE")
		fmt.Fprintln(os.Stderr, "!! -- !!")
		fmt.Fprintf(os.Stderr, `This tool checks the %s:%s.phenotype table to avoid duplicating data. 
However, it only looks at data that is already present in the BigQuery table. 
Any data that has been emitted to disk but not yet loaded into BigQuery cannot be seen by this tool.
If you pass multiple files to this tool, it will behave as expected.
If you run this tool multiple times (once per file) without first loading the output from the last round into BigQuery, you may get duplicated data.
!! -- !!
Please re-run this tool with the --ack flag to demonstrate that you understand this limitation.%s`, BQ.Project, BQ.Database, "\n")
		os.Exit(1)
	}

	log.Println("Processing", len(phenoPaths), "files:", phenoPaths)

	// Map out known FieldIDs so we don't add duplicate phenotype records
	// connect to BigQuery
	var err error
	BQ.Context = context.Background()
	BQ.Client, err = bigquery.NewClient(BQ.Context, BQ.Project)
	if err != nil {
		log.Fatalln("Connecting to BigQuery:", err)
	}
	defer BQ.Client.Close()
	knownFieldIDs, err := FindExistingFields(BQ)
	if err != nil {
		log.Fatalln(err)
	}
	log.Println("Found", len(knownFieldIDs), "FieldIDs already in the database, will ignore those FieldIDs in the new data, if present.")

	fmt.Printf("sample_id\tFieldID\tinstance\tarray_idx\tvalue\tcoding_file_id\n")

	// Process each file. KnownFieldIDs is modified within.
	for _, phenoPath := range phenoPaths {
		if err := ProcessOnePath(BQ, knownFieldIDs, phenoPath); err != nil {
			log.Fatalln(err)
		}
	}
}

func ProcessOnePath(BQ *WrappedBigQuery, knownFieldIDs map[string]struct{}, phenoPath string) error {
	headers := make(map[int]Column)

	phenoPath = genomisc.ExpandHome(phenoPath)
	log.Printf("Converting %s\n", phenoPath)

	// phenotype file

	f, err := os.Open(phenoPath)
	if err != nil {
		log.Fatalln(err)
	}
	defer f.Close()

	delim := genomisc.DetermineDelimiter(f)

	f.Seek(0, 0)

	// Buffered reader
	br := bufio.NewReaderSize(f, BufferSize)
	fileCSV := csv.NewReader(br)
	fileCSV.Comma = delim

	log.Printf("Determined phenotype delimiter to be \"%s\"\n", string(delim))

	// Map the headers
	headRow, err := fileCSV.Read()
	if err != nil {
		log.Fatalln("Header parsing error:", err)
	}
	if err := parseHeaders(BQ, headRow, headers); err != nil {
		log.Fatalln(err)
	}

	// Track which fields are new
	newFieldIDs := make(map[string]struct{})
	for col := range headRow {
		if col == 0 {
			continue
		}

		if _, exists := knownFieldIDs[headers[col].FieldID]; !exists {
			// This FieldID is new to us -- store it in the DB later
			newFieldIDs[headers[col].FieldID] = struct{}{}
		}
	}

	summarizeImport(knownFieldIDs, newFieldIDs)

	if len(newFieldIDs) < 1 {
		log.Println("No new fields were found in file", phenoPath, " -- skipping")
		return nil
	}

	var addedPhenoCount int
	sample := 1
	for ; ; sample++ {
		// Log to screen more frequently during profiling.
		if sample%1e10 == 0 {
			log.Println("Inserting sample", sample)
		}
		if sample%10000 == 0 {
			log.Println("Saw sample", sample)
		}

		row, err := fileCSV.Read()
		if err != nil && err == io.EOF {
			break
		} else if err != nil {
			log.Fatalln(err)
		}

		// Handle each sample

		sampleID := row[0]
		pheno := SamplePheno{}
		for col := range row {
			if col == 0 {
				// Skip the sample ID
				continue
			}

			pheno.Column = headers[col]
			pheno.SampleID = sampleID
			pheno.Value = row[col]

			if _, exists := knownFieldIDs[pheno.FieldID]; exists {
				// Skip any FieldIDs that were inserted on previous runs. Note
				// that this precludes updates on the database where the update
				// consists of adding new entries, such as a new array_idx to a
				// previously known FieldID.
				continue
			}

			if strings.TrimSpace(pheno.Value) == "" {
				// Skip blanks
				continue
			}

			fmt.Fprintf(STDOUT, "%s\t%s\t%s\t%s\t%s\t%v\n", pheno.SampleID, pheno.FieldID, pheno.Instance, pheno.ArrayIDX, pheno.Value, pheno.CodingFileID)

			addedPhenoCount++
		}
	}

	log.Println("Populated the table with", addedPhenoCount, "new records from", len(newFieldIDs), "previously unseen FieldIDs the phenotype file")

	// Make the tool aware of new fields so it won't import them from the next
	// file.
	for key := range newFieldIDs {
		knownFieldIDs[key] = struct{}{}
	}

	return nil
}

func summarizeImport(knownFieldIDs, newFieldIDs map[string]struct{}) {
	knownSlice := make([]int, 0, len(knownFieldIDs))
	for v := range knownFieldIDs {
		if _, exists := newFieldIDs[v]; !exists {
			// If we aren't trying to add this field, ignore
			continue
		}

		if intVal, err := strconv.Atoi(v); err == nil {
			knownSlice = append(knownSlice, intVal)
		}
	}
	sort.IntSlice(knownSlice).Sort()

	newSlice := make([]int, 0, len(newFieldIDs))
	for v := range newFieldIDs {
		if intVal, err := strconv.Atoi(v); err == nil {
			newSlice = append(newSlice, intVal)
		}
	}
	sort.IntSlice(newSlice).Sort()

	fmt.Fprintf(os.Stderr, "Total known fields: %v. %v fields in the file are duplicate and will *not* be added: %v\n", len(knownFieldIDs), len(knownSlice), splitToString(knownSlice, ","))
	fmt.Fprintf(os.Stderr, "%d fields are new and *will* be added: %v\n", len(newFieldIDs), splitToString(newSlice, ","))
}

// Yields a map that notes which columns are associated with which field
func fieldMap(headers map[int]Column) map[string][]int {
	// Map FieldID => ith, jth, kth, (etc) Column(s) in the CSV
	out := make(map[string][]int)

	for i, v := range headers {
		out[v.FieldID] = append(out[v.FieldID], i)
	}

	return out
}

type codingLookup struct {
	FieldID      int64              `bigquery:"FieldID"`
	CodingFileID bigquery.NullInt64 `bigquery:"coding_file_id"`
}

func parseHeaders(wbq *WrappedBigQuery, row []string, headers map[int]Column) error {
	for col, header := range row {
		dash := strings.Index(header, "-")
		dot := strings.Index(header, ".")

		// The very first field is just "eid"
		if dash == -1 || dot == -1 {
			if col > 0 {
				return fmt.Errorf("All columns after the 0th are expected to be of format XX-XX.XX, but column %d is %s", col, header)
			}
			headers[col] = Column{
				FieldID: header,
			}
			continue
		}

		// Other fields have all 3 parts
		headers[col] = Column{
			FieldID:  header[:dash],
			Instance: header[dash+1 : dot],
			ArrayIDX: header[dot+1:],
		}
	}

	// Now make it easy to figure out if we have a coding file to map to

	// Get the FieldID => coding_file_id map
	codinglookup := []codingLookup{}
	// if err := db.Select(&codinglookup, "SELECT FieldID, coding_file_id FROM dictionary WHERE coding_file_id IS NOT NULL"); err != nil {
	// 	return err
	// }
	query := wbq.Client.Query(fmt.Sprintf(`SELECT FieldID, coding_file_id
FROM %s.dictionary
WHERE coding_file_id IS NOT NULL`, wbq.Database))
	itr, err := query.Read(wbq.Context)
	if err != nil {
		return err
	}
	for {
		var values codingLookup

		err := itr.Next(&values)
		if err == iterator.Done {
			break
		}
		if err != nil {
			return err
		}
		codinglookup = append(codinglookup, values)
	}

	// Annotate each header column with the CodingFileID that applies to
	// the FieldID in that column.
	fm := fieldMap(headers)
	for _, lookup := range codinglookup {
		relevantHeaders := fm[strconv.FormatInt(lookup.FieldID, 10)]
		for _, v := range relevantHeaders {
			whichHeader := headers[v]
			whichHeader.CodingFileID = lookup.CodingFileID
			headers[v] = whichHeader
		}
	}

	return nil
}

func FindExistingFields(wbq *WrappedBigQuery) (map[string]struct{}, error) {
	// Map out known FieldIDs so we don't add duplicate phenotype records
	knownFieldIDs := make(map[string]struct{})

	query := wbq.Client.Query(fmt.Sprintf(`SELECT DISTINCT p.FieldID
	FROM %s.phenotype p
	WHERE p.FieldID IS NOT NULL
`, wbq.Database))
	itr, err := query.Read(wbq.Context)
	if err != nil && strings.Contains(err.Error(), "Error 404") {
		// Not an error; the table just doesn't exist yet
		return knownFieldIDs, nil
	} else if err != nil {
		return nil, err
	}
	for {
		var values struct {
			FieldID int64 `bigquery:"FieldID"`
		}
		err := itr.Next(&values)
		if err == iterator.Done {
			break
		}
		if err != nil {
			return nil, err
		}
		knownFieldIDs[strconv.FormatInt(values.FieldID, 10)] = struct{}{}
	}

	return knownFieldIDs, nil
}

// Via https://stackoverflow.com/a/42159097/199475
func splitToString(a []int, sep string) string {
	if len(a) == 0 {
		return ""
	}

	b := make([]string, len(a))
	for i, v := range a {
		b[i] = strconv.Itoa(v)
	}
	return strings.Join(b, sep)
}
