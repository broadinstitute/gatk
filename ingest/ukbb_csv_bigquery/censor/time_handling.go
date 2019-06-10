package main

import (
	"fmt"
	"time"
)

func TimeToUKBDate(t time.Time) string {
	if t.Equal(time.Time{}) {
		return NullMarker
	}

	return t.Format("2006-01-02")
}

func TimesToFractionalYears(earlier, later time.Time) string {
	if later.Before(earlier) {
		return NullMarker
	}
	y, m, d, h, min, sec := time_diff(earlier, later)

	return fmt.Sprintf("%.6f", float64(y)+float64(m)/12+float64(d)/(12*30)+float64(h)/(24*365)+float64(min)/(60*24*365)+float64(sec)/(60*60*24*365))
}

// Taken directly from https://stackoverflow.com/a/36531443/199475
func time_diff(a, b time.Time) (year, month, day, hour, min, sec int) {
	if a.Location() != b.Location() {
		b = b.In(a.Location())
	}
	if a.After(b) {
		a, b = b, a
	}
	y1, M1, d1 := a.Date()
	y2, M2, d2 := b.Date()

	h1, m1, s1 := a.Clock()
	h2, m2, s2 := b.Clock()

	year = int(y2 - y1)
	month = int(M2 - M1)
	day = int(d2 - d1)
	hour = int(h2 - h1)
	min = int(m2 - m1)
	sec = int(s2 - s1)

	// Normalize negative values
	if sec < 0 {
		sec += 60
		min--
	}
	if min < 0 {
		min += 60
		hour--
	}
	if hour < 0 {
		hour += 24
		day--
	}
	if day < 0 {
		// days in month:
		t := time.Date(y1, M1, 32, 0, 0, 0, 0, time.UTC)
		day += 32 - t.Day()
		month--
	}
	if month < 0 {
		month += 12
		year--
	}

	return
}
