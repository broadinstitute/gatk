package main

import (
	"context"
	"fmt"
	"log"
	"time"

	"cloud.google.com/go/bigquery"
)

func Censor(BQ *WrappedBigQuery, deathCensorDateString, phenoCensorDateString string) error {
	var err error
	var res map[int64]string
	var out = make(map[int64]CensorResult)

	// connect to BigQuery
	BQ.Context = context.Background()
	BQ.Client, err = bigquery.NewClient(BQ.Context, BQ.Project)
	if err != nil {
		log.Fatalln("Connecting to BigQuery:", err)
	}
	defer BQ.Client.Close()

	deathCensorDate, err := time.Parse("2006-01-02", deathCensorDateString)
	if err != nil {
		log.Fatalln(err)
	}

	phenoCensorDate, err := time.Parse("2006-01-02", phenoCensorDateString)
	if err != nil {
		log.Fatalln(err)
	}

	log.Printf("Censoring death at %v and missing phenotypes at %v\n", deathCensorDate, phenoCensorDate)

	// Find everyone and assign their enrollment date
	// Must come first.
	res, err = BigQuerySingleFieldFirst(BQ, 53)
	if err != nil {
		log.Fatalln(err)
	}
	now := time.Now()
	N := len(res)
	log.Println("Found", N, "enrollment results")
	if N == 0 {
		log.Fatalln("No enrollment dates found. Incidence cannot be computed without an enrollment date (FieldID 53).")
	}
	for k, v := range res {
		entry := out[k]
		entry.SampleID = k
		entry.computed = now
		entry.deathCensored = deathCensorDate
		entry.phenoCensored = phenoCensorDate
		entry.enrolled, err = time.Parse("2006-01-02", v)
		if err != nil {
			entry.Missing = append(entry.Missing, "enroll_date[53]")
			log.Println("Enrollment date parsing issue for", entry)
		}
		out[k] = entry
	}

	N, err = BQ.AddBirthYear(out)
	log.Println("Found", N, "birth year results")
	if N == 0 {
		log.Fatalln("No birth years found. Ages cannot be computed without at least a birth year (FieldID 34).")
	}
	if err != nil {
		log.Fatalln(err)
	}

	N, err = BQ.AddBirthMonth(out)
	log.Println("Found", N, "birth month results")
	if N == 0 {
		for k := range out {
			entry := out[k]
			entry.Missing = append(entry.Missing, "birth_month[52]")
			out[k] = entry
		}
		log.Println("Warning: 0 birth months found. Are you missing FieldID 52?")
	}
	if err != nil {
		log.Fatalln(err)
	}

	N, err = BQ.AddLostDate(out)
	log.Println("Found", N, "lost results")
	if N == 0 {
		for k := range out {
			entry := out[k]
			entry.Missing = append(entry.Missing, "lost_to_followup[191]")
			out[k] = entry
		}
		log.Println("Warning: 0 dates for loss to followup found. Are you missing FieldID 191?")
	}
	if err != nil {
		log.Fatalln(err)
	}

	N, err = BQ.AddDiedDate(out)
	log.Println("Found", N, "died results")
	if N == 0 {
		log.Println("Warning: 0 deaths found. Are you missing FieldID 40000?")
	}
	if err != nil {
		log.Fatalln(err)
	}

	// Print
	fmt.Printf("sample_id\tbirthdate\tenroll_date\tenroll_age\tdeath_date\tdeath_age\tdeath_censor_date\tdeath_censor_age\tphenotype_censor_date\tphenotype_censor_age\tlost_to_followup_date\tlost_to_followup_age\tcomputed_date\tmissing_fields\n")
	for _, v := range out {

		// Some samples have been removed
		if v.SampleID <= 0 {
			continue
		}

		fmt.Printf("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			v.SampleID,
			TimeToUKBDate(v.Born()),
			TimeToUKBDate(v.enrolled),
			TimesToFractionalYears(v.Born(), v.enrolled),
			TimeToUKBDate(v.died.Time),
			TimesToFractionalYears(v.Born(), v.died.Time),
			TimeToUKBDate(v.DeathCensored()),
			TimesToFractionalYears(v.Born(), v.DeathCensored()),
			TimeToUKBDate(v.PhenoCensored()),
			TimesToFractionalYears(v.Born(), v.PhenoCensored()),
			TimeToUKBDate(v.lost.Time),
			TimesToFractionalYears(v.Born(), v.lost.Time),
			TimeToUKBDate(v.computed),
			v.MissingToString(),
		)
	}

	return nil
}

func (BQ *WrappedBigQuery) AddBirthYear(out map[int64]CensorResult) (int, error) {
	// Birth year
	res, err := BigQuerySingleFieldFirst(BQ, 34)
	if err != nil {
		return 0, err
	}
	N := len(res)

	for k, v := range res {
		entry := out[k]
		entry.bornYear = v
		out[k] = entry
	}

	return N, nil
}

func (BQ *WrappedBigQuery) AddBirthMonth(out map[int64]CensorResult) (int, error) {
	// Birth month ==> we don't have it in every releas
	res, err := BigQuerySingleFieldFirst(BQ, 52)
	if err != nil {
		return 0, err
	}
	N := len(res)

	for k, v := range res {
		entry := out[k]
		entry.bornMonth = v
		out[k] = entry
	}

	return N, nil
}

func (BQ *WrappedBigQuery) AddLostDate(out map[int64]CensorResult) (int, error) {

	// Lost to follow-up
	res, err := BigQuerySingleFieldFirst(BQ, 191)
	if err != nil {
		return 0, err
	}
	N := len(res)

	for k, v := range res {
		entry := out[k]

		lostDate, err := time.Parse("2006-01-02", v)
		if err != nil {
			log.Println("Lost to followup date parsing issue for", entry)
		}
		if err := entry.lost.Scan(lostDate); err != nil {
			log.Println("Lost to followup date parsing issue for", entry)
		}

		out[k] = entry
	}

	return N, nil
}

func (BQ *WrappedBigQuery) AddDiedDate(out map[int64]CensorResult) (int, error) {

	// Died
	res, err := BigQuerySingleFieldFirst(BQ, 40000)
	if err != nil {
		return 0, err
	}
	N := len(res)

	for k, v := range res {
		entry := out[k]

		diedDate, err := time.Parse("2006-01-02", v)
		if err != nil {
			log.Println("Died date parsing issue for", entry)
		}
		if err := entry.died.Scan(diedDate); err != nil {
			log.Println("Died date parsing issue for", entry)
		}

		out[k] = entry
	}

	return N, nil
}
