package main

import (
	"bufio"
	"context"
	"flag"
	"log"
	"os"
	"path"
	"strings"

	"cloud.google.com/go/bigquery"
)

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

var materializedDB string

func main() {
	defer STDOUT.Flush()

	var BQ = &WrappedBigQuery{
		Context: context.Background(),
	}
	var tabfile string
	var displayQuery bool
	var override bool
	var diseaseName string

	flag.StringVar(&BQ.Project, "project", "", "Google Cloud project you want to use for billing purposes only")
	flag.StringVar(&BQ.Database, "database", "", "BigQuery source database name (note: must be formatted as project.database, e.g., broad-ml4cvd.ukbb7089_201904)")
	flag.StringVar(&tabfile, "tabfile", "", "Tabfile-formatted phenotype definition")
	flag.StringVar(&materializedDB, "materialized", "broad-ml4cvd.ukbb7089_201904", "project.database storing materialized view tables")
	flag.BoolVar(&displayQuery, "display-query", false, "Display the constructed query and exit?")
	flag.BoolVar(&override, "override", false, "Force run, even if this tool thinks your tabfile is inadequate?")
	flag.StringVar(&diseaseName, "disease", "", "If not specified, the tabfile will be parsed and become the disease name.")
	flag.Parse()

	if BQ.Project == "" || BQ.Database == "" || tabfile == "" || materializedDB == "" {
		flag.PrintDefaults()
		os.Exit(1)
	}

	tabs, err := ParseTabFile(tabfile)
	if err != nil {
		log.Fatalln(err)
	}

	if diseaseName == "" {
		diseaseName = path.Base(tabfile)
		if parts := strings.Split(diseaseName, "."); len(parts) > 1 {
			diseaseName = strings.Join(parts[0:len(parts)-1], ".")
		}
	}

	log.Println("Processing disease", diseaseName)

	missingFields, err := tabs.CheckSensibility()
	if err != nil && !override {
		log.Println(err)
		log.Fatalf("%s: Add the missing fields to your tabfile, or re-run with the -override flag to process anyway.\n", diseaseName)
	} else if err != nil && override {
		log.Println(diseaseName, err)
		log.Printf("%s: Overriding error check for missing fields and continuing.\n", diseaseName)
	}

	BQ.Client, err = bigquery.NewClient(BQ.Context, BQ.Project)
	if err != nil {
		log.Fatalln("Connecting to BigQuery:", err)
	}
	defer BQ.Client.Close()

	query, err := BuildQuery(BQ, tabs, displayQuery)
	if err != nil {
		log.Fatalln(diseaseName, err)
	}

	if displayQuery {
		return
	}

	if err := ExecuteQuery(BQ, query, diseaseName, missingFields); err != nil {
		log.Fatalln(diseaseName, err)
	}
}
