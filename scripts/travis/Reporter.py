# Usage: python Reporter.py <log url>
#
# This script creates a comment on github using the API token provided through GITHUB_API_TOKEN
# It first checks if a comment already exists which was made with the given identity
# and which mentions the same travis build ID as the current running one in the first line of the comment.
# If no comment is found it creates a new one with a header and the log url.
# If a matching comment is found this appends it's log information to the existing comment.
#
# There exist race conditions when running this in parallel which can result in shards being missed if this is called
# multiple times simultaneously.

import os
import sys

from github import Github

key = os.environ["GITHUB_API_TOKEN"]
test_type = os.getenv("TEST_TYPE", "")

g = Github(key)

github_login = g.get_user().login

pull_number = os.environ["TRAVIS_PULL_REQUEST"]
repo_name = os.environ["TRAVIS_REPO_SLUG"]
job_number = os.environ["TRAVIS_JOB_NUMBER"]
job_id = os.environ["TRAVIS_JOB_ID"]
build_number = os.environ["TRAVIS_BUILD_NUMBER"]
build_id = os.environ["TRAVIS_BUILD_ID"]
jdk_version = os.getenv("TRAVIS_JDK_VERSION","")
travis_page_url = os.getenv("TRAVIS_BUILD_WEB_URL")
job_page_url = os.getenv("TRAVIS_JOB_WEB_URL")

repo = g.get_repo(repo_name)
pull = repo.get_pull(int(pull_number))
comments = pull.get_issue_comments()

log_url = sys.argv[1]
message = "| %s | %s | [%s](%s) | [logs](%s) |" % (test_type, jdk_version, job_number, job_page_url, log_url)


def matches_build_id(comment, build_id):
    body = comment.body
    lines = body.splitlines()
    return build_id in lines[0]


def is_my_comment(comment):
    return comment.user.login == github_login


def update(comment):
    comment.edit(comment.body + "\n" + message)


def new_comment(pull):
    pull.create_issue_comment("Travis reported job failures from build [%s](%s)" % (build_number, travis_page_url)
                              + "\nFailures in the following jobs:"
                              + "\n"
                              + "\n | Test Type | JDK | Job ID | Logs |"
                              + "\n | --------- |---- | ------ | ---- |"
                              + "\n" + message)


# Update the comment or create a new one
for comment in comments:
    if is_my_comment(comment) and matches_build_id(comment, build_id):
        print("found matching comment")
        update(comment)
        break
else:
    print("creating new comment")
    new_comment(pull)
