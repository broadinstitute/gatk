package main

import (
	"math/rand"
	"time"
)

func init() {
	// Ensure different folder names on each run
	rand.Seed(time.Now().UTC().UnixNano())
}

// RandOrthoglyphs produces a string of length n randomly.
func RandOrthoglyphs(n int) string {
	var letters = []rune("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
	lenLetters := len(letters)
	b := make([]rune, n)
	for i := range b {
		b[i] = letters[rand.Intn(lenLetters)]
	}
	return string(b)
}
