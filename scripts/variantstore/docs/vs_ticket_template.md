# VS Board Ticket Template for LLMs

## Instructions for LLM

When creating a Jira ticket for the VS (Variant-Team) board, follow the conversational style used by the team. Tickets are written in plain prose, not with formal markdown headers.

---

## Format

**Summary:** `[Optional prefix] Concise action phrase`
- Use `[Spike]` prefix for investigations
- Use imperative voice for tasks: "Update X to do Y", "Add support for Z"
- For bugs, state the problem: "X does not respect Y"

**Description:**
```
One or two paragraphs explaining:
- What the problem or need is
- Context the reader needs (what currently happens, why it matters)
- Any approach discussed in refinement (if applicable)

Notes: (optional—use if there are technical details, links, contacts, or caveats)
• Note 1
• Note 2

AC: (or A/C:)
• Concrete, testable outcome 1
• Concrete, testable outcome 2
```

---

## Key Patterns

- **Prose first**: Start with 1–2 paragraphs of context. No headers like "## Summary" or "## Background".
- **Notes section**: Use a simple "Notes:" label followed by bullets when you need to capture links, contacts, or technical caveats.
- **AC format**: Use `AC:` or `A/C:` followed by bullet points. Not a markdown header.
- **Testable bullets**: AC items should be concrete and verifiable, not vague.
- **Spike timeboxing**: For investigations, note if it should be timeboxed (e.g., "Since this is a nice to have, timebox this investigation. If this is a large lift, we will want to discuss further with Lee.")
- **Links inline**: Include hyperlinks and mailto: addresses directly in the text where relevant.

---

## Do / Don't

| Do | Don't |
|---|---|
| Write in plain conversational prose | Use formal markdown headers (##) |
| Start with context paragraphs | Jump straight to bullet lists |
| Use `AC:` or `A/C:` for acceptance criteria | Use "Acceptance Criteria" as a header |
| Write AC as testable outcomes | Write vague AC ("should work well") |
| Include links, contacts, references in Notes | Omit context that helps future readers |
| For spikes, mention what artifact or decision is expected | Leave the output of a spike ambiguous |

---

## Example 1: Bug/Feature (from VS-1747)

**Summary:** GvsExtractCohortFromSampleNames does not respect the withdrawn field

**Description:**
```
If a user passes in a sample id that is in the BQ tables but marked as withdrawn, we will currently include the sample in the extracted vcf.

In refinement we discussed adding a flag to respectWithdrawn and set the default to the current behavior while we investigate if we want that behavior to continue. This will enable us to respect the withdrawn flag for subsets like lrWGS in the meantime, and easily update the default later.

AC:
• When we are respecting the flag, the workflow checks early if a sample is marked as withdrawn or is not in the dataset, and then throws a clear error promptly.
• Samples that are marked as withdrawn are not exported in vcfs produced by GvsExtractCohortFromSampleNames, even if explicitly listed in the sample names.
```

---

## Example 2: Technical Fix (from VS-1772)

**Summary:** Update prepare code to properly pad intervals

**Description:**
```
Extract inserts padding before (and possibly after?) intervals in order to catch overlapping ref blocks and indels. The prepare logic doesn't currently do so, resulting in inaccurate results. Update the prepare logic to add padding to intervals that is consistent with what we do in extract.

AC:
• The prepare logic that uses intervals correctly inserts padding around them in ways that are consistent with extract.
• Verify that the results of extracting after running prepare on a small region are identical to running prepare on the full data and only extracting that small interval.
```

---

## Example 3: Spike with Notes (from VS-1520)

**Summary:** [Spike] LOFtee annotation for VAT

**Description:**
```
A P2 for Foxtrot is to add loss of function annotations to the VAT. We should investigate if this is feasible as Lee warns it could be complicated.

Notes:
• LOFtee does not have a canonical docker image. This image may work, but it will likely need some upgrades/changes, since it is already about 3 years old.
• LOFtee is written by Konrad Karczewski of ATGU.
• gnomAD/Hail has cached all possible SNP LOFtee results. This can be used for faster annotation. From Konrad: "gnomAD has pre-computed LOFTEE (and VEP) for all variants in gnomAD. but they've also pre-computed all POSSIBLE SNVs (3 billion x 3 possible mutations at each base = 9 billion rows), which can be used as a lookup table. then all you need to VEP [and LOFtee] is the indels"
• If we can get the information into a custom table, we should be able to feed it into NIRVANA per Lee.

A/C:
• Since this is a nice to have, timebox this investigation. If this is a large lift, we will want to discuss further with Lee.
• Make a VAT for the 3k dataset and share with scientific owners for review.
```

---

## Example 4: Research/Scoping Task (from VS-1745)

**Summary:** Spike: Propose which annotations we will add to the VAT post v9

**Description:**
```
This spreadsheet has a list of annotations and priority. Add a column and cross check which are already in the VAT. Of the ones that are not already in the VAT, look at the remaining high and medium priority annotations and scope out the effort for each. From there, have a team discussion and we will make a proposal list for Lee on which annotations we will add between v9 and v10.

Possibly helpful: VAT design doc, Epic for prep for foxtrot epic.
```
