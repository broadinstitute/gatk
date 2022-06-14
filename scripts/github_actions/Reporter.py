# Usage: python Reporter.py <log url>
#
# This script creates a comment on github using the API token provided through GITHUB_API_TOKEN
# It first checks if a comment already exists which was made with the given identity
# and which mentions the same github actions build ID as the current running one in the first line of the comment.
# If no comment is found it creates a new one with a header and the log url.
# If a matching comment is found this appends it's log information to the existing comment.
#
# There exist race conditions when running this in parallel which can result in shards being missed if this is called
# multiple times simultaneously.

import os
import sys

from github import Github

key = os.environ["GITHUB_TOKEN"]
test_type = os.getenv("TEST_TYPE", "")

g = Github(key)

github_login = g.get_user().login

pull_number = os.environ["GITHUB_EVENT_NUMBER"]
repo_name = os.environ["GITHUB_REPOSITORY"]
job_number = os.environ["JOB_MATRIX_ID"]
job_id = os.environ["GITHUB_RUN_ID"]
build_id = os.environ["GITHUB_RUN_ID"]
jdk_version = os.getenv("JDK_VERSION","")
job_page_url = os.getenv("GITHUB_LOGS_URL")
workflow_url = os.getenv("ACTIONS_JOB_WEB_URL")

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
    pull.create_issue_comment("Github actions tests reported job failures from actions build [%s](%s)" % (build_id, workflow_url)
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
