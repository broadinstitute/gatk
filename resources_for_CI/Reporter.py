import os
import sys

from github import Github

key = os.environ["GITHUB_TOKEN"]
g = Github(key)

ME = g.get_user().login

pull_number = os.environ["TRAVIS_PULL_REQUEST"]
repo_name = os.environ["TRAVIS_REPO_SLUG"]
job_number = os.environ["TRAVIS_JOB_NUMBER"]
job_id = os.environ["TRAVIS_JOB_ID"]
build_number = os.environ["TRAVIS_BUILD_NUMBER"]
build_id = os.environ["TRAVIS_BUILD_ID"]
jdk_version = os.getenv("TRAVIS_JDK_VERSION","")
test_type = os.getenv("TEST_TYPE", "")

repo = g.get_repo(repo_name)
pull = repo.get_pull(int(pull_number))
comments = pull.get_issue_comments()


def matches_build_number(comment, build_number):
    body = comment.body
    lines = body.splitlines()
    return build_number in lines[0]


def is_my_comment(comment):
    return comment.user.login == ME


travis_page_url = "https://travis-ci.com/" + repo_name + "/builds/" + build_id
job_page_url = "https://travis-ci.com/" + repo_name + "/jobs/" + job_id

log_url = sys.argv[1]
message = "| %s | %s | [%s](%s) | [logs](%s) |" % (test_type, jdk_version, job_number, job_page_url, log_url)


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
    if is_my_comment(comment) and matches_build_number(comment, build_number):
        print("found matching comment")
        update(comment)
        break
else:
    print("creating new comment")
    new_comment(pull)
