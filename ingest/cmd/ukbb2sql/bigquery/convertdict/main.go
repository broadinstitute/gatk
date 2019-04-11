package main

import (
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	"strings"
)

const (
	ExpectedRows = 17
)

func main() {
	var (
		dictPath string
	)

	flag.StringVar(&dictPath, "dict", "https://biobank.ndph.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv", "URL to CSV file with the UKBB data dictionary")
	flag.Parse()

	if dictPath == "" {
		flag.PrintDefaults()
		log.Fatalln()
	}

	if err := ImportDictionary(dictPath); err != nil {
		log.Fatalln(err)
	}
}

func ImportDictionary(url string) error {
	log.Printf("Importing from %s\n", url)

	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	reader := csv.NewReader(resp.Body)
	reader.Comma = ','
	reader.LazyQuotes = true

	header := make([]string, 0)
	j := 0
	for ; ; j++ {
		row, err := reader.Read()
		if err != nil && err == io.EOF {
			resp.Body.Close()
			break
		} else if err != nil {
			buf := bytes.NewBuffer(nil)
			io.Copy(buf, resp.Body)
			if strings.Contains(buf.String(), "internal error") {
				log.Println("Dictionary File is not permitted to be downloaded from the UKBB")
				continue
			}
		}

		// Handle the header
		if j == 0 {
			log.Printf("Header (%d elements): %+v\n", len(row), row)
			header = append(header, row...)
			for k, v := range header {
				if v == "Coding" {
					header[k] = "coding_file_id"
					break
				}
			}

			if nCols := len(header); nCols != ExpectedRows {
				return fmt.Errorf("Expected a CSV with %d columns; got one with %d", ExpectedRows, nCols)
			}

			fmt.Println(strings.Join(header, "\t"))

			continue
		}

		// Handle the entries
		if len(row) == ExpectedRows {
			fmt.Println(strings.Join(row, "\t"))
		}
	}

	log.Println("Created dictionary file with", j, "entries")

	return nil
}
