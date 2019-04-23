package main

import (
	"fmt"
	"strings"
	"time"

	"gopkg.in/guregu/null.v3"
)

type CensorResult struct {
	SampleID int64

	// In the database, populate with a list of fields that we would like to
	// have (e.g., month of birth, lost to followup) but which were not present,
	// so we know if the table was constructed from incomplete data
	Missing []string

	// Guaranteed
	enrolled time.Time
	computed time.Time // Date this was computed

	// May be null appropriately
	died null.Time
	lost null.Time

	// Unsure
	phenoCensored time.Time
	deathCensored time.Time

	// convenience / not exported
	bornYear  string
	bornMonth string
}

func (s CensorResult) Born() time.Time {
	// If we know year + month, then neutral assumption is that birthday is on
	// the middle day of the month. If we just know year, then assumption is
	// being born midway through the year (July 2).
	month := s.bornMonth
	day := "15"

	if month == "" {
		month = "7"
		day = "02"
	}

	dt, err := time.Parse("2006-01-02", fmt.Sprintf("%04s-%02s-%02s", s.bornYear, month, day))
	if err != nil {
		return time.Time{}
	}

	return dt
}

func (s CensorResult) DiedString() string {
	if !s.died.Valid {
		return NullMarker
	}

	return TimeToUKBDate(s.died.Time)
}

func (s CensorResult) DeathCensored() time.Time {
	if s.died.Valid {
		return s.died.Time
	}

	if s.lost.Valid {
		return s.lost.Time
	}

	return s.deathCensored
}

func (s CensorResult) PhenoCensored() time.Time {
	if s.died.Valid {
		return s.died.Time
	}

	if s.lost.Valid {
		return s.lost.Time
	}

	return s.phenoCensored
}

func (s CensorResult) MissingToString() string {
	if res := strings.Join(s.Missing, "|"); len(res) > 0 {
		return res
	}

	return "NA"
}
