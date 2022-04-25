// Update GATK User Stats Spreadsheet
// Author: Jonn Smith
// Date: 2022 Jan 04
//
// Record the number of pulls for our GATK docker container(s).
// Designed to be run periodically.
//
// Inspired by: https://faun.pub/track-how-many-times-docker-image-was-pulled-55cafe67ac97 
//
// Script URL: https://script.google.com/home/projects/1UQnN7Y7Ci8bP2QSczLdNerLUKh6vAmh1gCc40xEqiIGiPngfCLYeImbq/edit

const OVERALL_STATS_SHEET_NAME = "Overall Stats";
const DOCKER_STATS_SHEET_NAME = "Docker Repo Pulls";
const GITHUB_STATS_SHEET_NAME = "Github Stats";
const REQUIRED_SHEETS = [OVERALL_STATS_SHEET_NAME, DOCKER_STATS_SHEET_NAME, GITHUB_STATS_SHEET_NAME];

function updateSpredsheet() {

  // Set up our variables for this instance of this script:

  // NOTE: THIS ORDER MATTERS!  APPEND!  DO NOT INSERT!
  const dockerHubImages = ['broadinstitute/gatk', 'broadinstitute/gatk3'];
  const gitHubUrl = "https://api.github.com/repos/broadinstitute/gatk/releases";

  // Auth token to pull down the info on the releases so we don't get rate-limited: 
  // This is OPTIONAL, but HIGHLY RECOMMENDED!
  const GITHUB_AUTH_TOKEN = "";

  // ============================================

  // Make sure the spreadsheet is ready:
  setupSpreadsheetIfNecessary(dockerHubImages);

  // ============================================

  // Update the docker image downloads:
  var dockerImageDownloadCountList = recordDockerImagePullCount(dockerHubImages);

  // Update our github release downloads:
  var githubDownloads = recordGitHubReleaseDownloadCount(gitHubUrl, GITHUB_AUTH_TOKEN);

  // ============================================

  // Update our summary sheet:
  var spreadsheet = SpreadsheetApp.getActive().getSheetByName(OVERALL_STATS_SHEET_NAME);
  
  var downloadCounts = [];
  downloadCounts = downloadCounts.concat(dockerImageDownloadCountList);
  downloadCounts.push(githubDownloads);
  var totalDownloads = 0;
  downloadCounts.forEach(function (count) {
    totalDownloads += count;
  });

  // This is the last row with data in it:
  var row = [new Date()];
  row = row.concat(downloadCounts);
  row.push(totalDownloads);

  spreadsheet.appendRow(row);

  // Now make all the columns auto-fit so we can view them correctly:
  REQUIRED_SHEETS.forEach(function (sheetName) {
    const sheet = SpreadsheetApp.getActiveSpreadsheet().getSheetByName(sheetName);
    sheet.autoResizeColumns(1, sheet.getLastColumn());
  })

  // Finally set the active sheet to the overall stats page:
  SpreadsheetApp.getActiveSpreadsheet().setActiveSheet(
    SpreadsheetApp.getActiveSpreadsheet().getSheetByName(OVERALL_STATS_SHEET_NAME)
  );
}

function setupSpreadsheetIfNecessary(dockerHubImages) {
  // Set up the sheets / tabs, headers, etc. required by this script.

  const spreadsheet = SpreadsheetApp.getActiveSpreadsheet();

  // Check to see if we have to set up our sheets:
  var mustInitializeSheets = false;
  for (let i = 0; i < REQUIRED_SHEETS.length; i++) {
    // Create the sheet itself:
    if (spreadsheet.getSheetByName(REQUIRED_SHEETS[i]) == null) { 
      mustInitializeSheets = true;
      break;
    }
  }

  if (!mustInitializeSheets) {
    console.log("Sheet already set up for data.");
    return;
  }
  else {
    console.log("Must initialize sheets!");
  }

  // Make sure all our required sheets exist:
  REQUIRED_SHEETS.forEach( function (sheetName) {
    console.log("Adding sheet: " + sheetName);
    // Create the sheet itself:
    if (spreadsheet.getSheetByName(sheetName) == null) { 
      const newSheet = spreadsheet.insertSheet(sheetName);

      // Set up the top left cell as our date:
      const dateRange = newSheet.getRange("A1");
      dateRange.setValue("Date / Time");
      dateRange.setFontWeight('bold');
    }
  })

  // Set up the headers in the overall stats sheet next:
  const overallStatsSheet = spreadsheet.getSheetByName(OVERALL_STATS_SHEET_NAME);

  // Start with docker image headers:
  var headerNum = 0;
  dockerHubImages.forEach(function (dockerImage) {
    console.log("Setting header for docker image: " + dockerImage);
    var headerCell = overallStatsSheet.getRange(1, 2 + headerNum, 1, 1);

    headerCell.setValue(dockerImage + " Docker Pulls");
    headerCell.setFontWeight('bold');

    headerNum += 1;
  })

  // Now the github header:
  console.log("Setting header for GitHub releases");
  var headerCell = overallStatsSheet.getRange(1, 2 + headerNum, 1, 1);
  headerCell.setValue("Github Repo Pulls");
  headerCell.setFontWeight('bold');
  headerNum += 1;

  // Now our total downloads:
  console.log("Setting header for Total Downloads");
  headerCell = overallStatsSheet.getRange(1, 2 + headerNum, 1, 1);
  headerCell.setValue("Total Overall Downloads");
  headerCell.setFontWeight('bold');
  headerNum += 1;

  // Clean up the default sheet if it's still there:
  const defaultSheet = spreadsheet.getSheetByName("Sheet1");
  if (defaultSheet != null) {
    console.log("Removing default sheet.");
    spreadsheet.deleteSheet(defaultSheet);
  }
}

function cleanGithubUrl(url) {
  // Convert a regular github url to an API url so that
  // we can receive JSON results.

  const urlRegex = /http.*?github.com\/(.*)/;

  // Make sure the interpreter knows URL should be a string:
  url = String(url);

  // Set up a default return value:
  var newUrl = url;

  if (url.startsWith("https://github.com") || url.startsWith("http://github.com") || url.startsWith("https://www.github.com") || url.startsWith("http://www.github.com")) {
    const match = url.match(urlRegex);
    newUrl = "https://api.github.com/repos/" + match[1]
  }

  if (!newUrl.endsWith("releases") && !newUrl.endsWith("releases/")) {
    if (newUrl.endsWith("/")) {
      newUrl = newUrl + "releases"
    }
    else{
      newUrl = newUrl + "/releases"
    }
  }

  return newUrl;
}

function recordGitHubReleaseDownloadCount(githubUrl, authToken) {

  // Header for authorization to pull down the info on the releases so we don't get rate-limited: 
  // This is OPTIONAL, but HIGHLY RECOMMENDED!
  var REQUEST_HEADERS = {
    "headers" : {
      "Authorization" : "token " + authToken
    }
  };

  var NUM_FIELDS_PER_RELEASE = 2;
  var COUNT_START_COLUMN = 2;

  var numResultsPerPage = 100;
  var githubUrl = cleanGithubUrl(githubUrl);

  var row = [new Date()];
  
  // Get the spreadsheet with our docker info in it:
  var spreadsheet = SpreadsheetApp.getActive().getSheetByName(GITHUB_STATS_SHEET_NAME);

  // This is the last row with data in it:
  var lastRow = spreadsheet.getLastRow();

  // Get info on all our releases.
  // We'll need to paginate here, so let's do something very simple 
  // when we have to get "more" rows:
  var releaseInfo = [];
  var pageNum = 1;
  while (true) {
    var response = "";
    if (authToken.length > 0) {
      response = UrlFetchApp.fetch(githubUrl + "?per_page=" + numResultsPerPage + "&page=" + pageNum, REQUEST_HEADERS); 
    }
    else {
      response = UrlFetchApp.fetch(githubUrl + "?per_page=" + numResultsPerPage + "&page=" + pageNum); 
    }
    var githubJsonResponse = JSON.parse(response.getContentText());

    console.log("Releases found on page %d: %d", pageNum, githubJsonResponse.length);

    // If we've reached the last page, we won't have anything here:
    if (githubJsonResponse.length == 0) {
      break;
    }

    releaseInfo = releaseInfo.concat(githubJsonResponse);
    pageNum++;
  }
  console.log("Total number of releases found: %d", releaseInfo.length);

  // Sort our releases so we can put them into the spreadsheet 
  // and we don't have to worry about new releases later:
  releaseInfo.sort((a,b) => {
    if (a['published_at'] < b['published_at']) {
      return -1;
    }
    else if (a['published_at'] > b['published_at']) {
      return 1;
    }
    else {
      return 0;
    }
  });

  // Get the headers of our spreadsheet so we can check to see if our releases are all in there.
  // Otherwise we may have to add in a new header for a new release.

  var lastCol = spreadsheet.getLastColumn();
  var releaseVersions = [];
  if ((lastCol != 0) && (lastCol != 1)) { 
    releaseVersions = spreadsheet.getRange(1,2,1, lastCol-1).getValues()[0];
  }

  // Populate our new row of data:
  var releaseIndex = 0;
  var totalDownloads = 0;
  releaseInfo.forEach( function (release) {
    var tagName = release['tag_name'];
    var downloadCount = 0;
    if (release['assets'].length != 0) {
      var downloadCount = release['assets'][0]['download_count'];
    }

    // Add our new download count to our row:
    row.push(downloadCount);

    // Now check if we have this entry already in our headers.
    // We do this with the assumption that releases are never 
    // removed from the repo, so we can use simple counts to track
    // when new releases are added:
    var oldDownloadCount = 0;

    // Need to divide the release versions by 2 
    // because they comprise all the headers (including the new downloads)
    if ((releaseIndex + 1) > (releaseVersions.length/2)) {

      // New release!
      console.log("New release detected: %s", tagName)

      // We need to add this new release to our headers!
      console.log("Setting header cell: %d, %d) to track %s total downloads", 1, COUNT_START_COLUMN + (releaseIndex * NUM_FIELDS_PER_RELEASE), tagName)
      var newReleaseDownloadCountHeaderCell = spreadsheet.getRange(1, COUNT_START_COLUMN + (releaseIndex * NUM_FIELDS_PER_RELEASE), 1, 1);
      newReleaseDownloadCountHeaderCell.setValue(tagName);
      newReleaseDownloadCountHeaderCell.setFontWeight('bold');

      console.log("Setting header cell: %d, %d) to track %s new downloads", 1, COUNT_START_COLUMN + (releaseIndex * NUM_FIELDS_PER_RELEASE) + 1, tagName)
      var newReleaseDownloadDiffCountHeaderCell = spreadsheet.getRange(1, COUNT_START_COLUMN + (releaseIndex * NUM_FIELDS_PER_RELEASE) + 1, 1, 1);
      newReleaseDownloadDiffCountHeaderCell.setValue(tagName + " New Downloads");
      newReleaseDownloadDiffCountHeaderCell.setFontWeight('bold');

      // no old data so use today's downloads
      row.push(downloadCount);
    }
    else {
      // We've seen this release before.  Let's get some diff numbers:
      oldDownloadCount = spreadsheet.getRange(lastRow, COUNT_START_COLUMN + (releaseIndex * NUM_FIELDS_PER_RELEASE)).getValue();
      row.push(downloadCount - oldDownloadCount);
    }

    console.log("Release: %s: %d (delta=%d)",tagName, downloadCount, downloadCount - oldDownloadCount);
    releaseIndex++;

    totalDownloads += downloadCount;
  });

  // Update the row for our docker image counts:
  spreadsheet.appendRow(row);

  return totalDownloads;
}

// ========================================================================

function recordDockerImagePullCount(images) {

  var NUM_FIELDS_PER_IMAGE = 2;
  var COUNT_START_COLUMN = 2;

  var row = [new Date()];
  
  // Get the spreadsheet with our docker info in it:
  var spreadsheet = SpreadsheetApp.getActive().getSheetByName(DOCKER_STATS_SHEET_NAME);

  // This is the last row with data in it:
  var lastRow = spreadsheet.getLastRow();

  // Get the headers of our docker spreadsheet so we can update as necessary:
  var lastCol = spreadsheet.getLastColumn();
  var imageVersions = [];
  if (lastCol != 1) { 
    imageVersions = spreadsheet.getRange(1,2,1, lastCol-1).getValues()[0];
  }

  var totalPullCountList = [];
  var imageNum = 0;
  images.forEach(function (image) {

    // Get the new count:
    var pull_count = get_image_pull_count(image);
    row.push(pull_count);

    var new_pulls = 0;

    // Check if we need to add a new header:
    if ((imageNum + 1) > (imageVersions.length/2)) {
      // New Docker Image to Track!
      console.log("New docker image detected: %s", image)

      // We need to add this new release to our headers!
      console.log("Setting header cell: %d, %d) to track %s total pulls", 1, COUNT_START_COLUMN + (imageNum * NUM_FIELDS_PER_IMAGE), image)
      var newReleaseDownloadCountHeaderCell = spreadsheet.getRange(1, COUNT_START_COLUMN + (imageNum * NUM_FIELDS_PER_IMAGE), 1, 1);
      newReleaseDownloadCountHeaderCell.setValue("Pulls (" + image + ")");
      newReleaseDownloadCountHeaderCell.setFontWeight('bold');

      console.log("Setting header cell: %d, %d) to track %s new pulls", 1, COUNT_START_COLUMN + (imageNum * NUM_FIELDS_PER_IMAGE) + 1, image)
      var newReleaseDownloadDiffCountHeaderCell = spreadsheet.getRange(1, COUNT_START_COLUMN + (imageNum * NUM_FIELDS_PER_IMAGE) + 1, 1, 1);
      newReleaseDownloadDiffCountHeaderCell.setValue("New Pulls (" + image + ")");
      newReleaseDownloadDiffCountHeaderCell.setFontWeight('bold');
    }
    else {
      // Get the old count as well:
      var oldPullCount = spreadsheet.getRange(lastRow, COUNT_START_COLUMN + (imageNum * NUM_FIELDS_PER_IMAGE)).getValue();
      new_pulls = pull_count - oldPullCount;
    }

    // Add the count diff to our row:
    row.push(new_pulls);
    imageNum++;

    totalPullCountList.push(pull_count);
  });
  
  // Update the row for our docker image counts:
  spreadsheet.appendRow(row);

  return totalPullCountList;
}

function get_image_pull_count(image) {
  var response = UrlFetchApp.fetch("https://hub.docker.com/v2/repositories/" + image);
  var imageStats = JSON.parse(response.getContentText());
  return imageStats['pull_count'];
}

