#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
   ## NOTE ##
   This module is one of several for which the Master copy is in the
   UM repository. When making changes, please ensure the changes are made in
   the UM repository or they will be lost during the release process when the UM
   copy is copied over.

   Script to process the results of a suite and write a summary to file. The
   summary is in Trac wiki mark-up. Any projects that do not have a local
   mirror repository are assumed not to be used at that site and are
   excluded from the report.

   Owner: Scientific Software Development and Deployment team
          (formerly : UM System Development Team)
   Cylc Suite Syntax: shutdown handler = "suite_report.py"
   Command Line syntax:
       suite_report.py -S <suite_dir> [-v] [-q] [-N] [-L <log_dir>]
"""

import glob
import os
import re
import sqlite3
import sys
import traceback
import time
import subprocess
from optparse import OptionParser, OptionGroup
from collections import defaultdict
from fcm_bdiff import get_branch_diff_filenames


CYLC_SUITE_ENV_FILE = "cylc-suite-env"
PROCESSED_SUITE_RC = "suite.rc.processed"
ROSE_SUITE_RUN_CONF = "rose-suite-run.conf"
SUITE_DB_FILENAME = "cylc-suite.db"
TRAC_LOG_FILE = "trac.log"
DEFAULT_VERBOSITY = 3

PINK_FAIL_TEXT = "'''[[span(style=color: #FF00FF, pink failure )]]'''"
DESIRED_ORDER = [PINK_FAIL_TEXT, "failed", "succeeded"]

FCM = {
    "meto": "fcm",
    "ecmwf": "fcm",
    "nci": "fcm",
    "bom": "fcm",
    "uoe": "fcm",
    "niwa": "fcm",
    "kma": "fcm",
    "vm": "fcm",
    "jasmin": "fcm",
    "cehwl1": "fcm",
    "mss": "fcm",
    "ncas": "fcm",
    "psc": "fcm",
    "Unknown": "true",
}
RESOURCE_MONITORING_JOBS = {
    "meto": [
        "atmos-xc40_cce_um_fast_omp-seukv-4x9-noios-2t",
    ],
    "ecmwf": [],
    "nci": [],
    "bom": [],
    "uoe": [],
    "niwa": [],
    "kma": [],
    "vm": [],
    "jasmin": [],
    "cehwl1": [],
    "mss": [],
    "ncas": [],
    "psc": [],
    "Unknown": [],
}
CYLC_REVIEW_URL = {
    "meto": "http://fcm1/cylc-review",
    "ecmwf": "Unavailable",
    "nci": "http://accessdev.nci.org.au/cylc-review",
    "bom": "http://scs-watchdog-dev/rose-bush",
    "uoe": "Unavailable",
    "niwa": "http://w-rose01.maui.niwa.co.nz/cylc-review",
    "kma": "Unavailable",
    "vm": "http://localhost/cylc-review",
    "jasmin": "Unavailable",
    "cehwl1": "Unavailable",
    "mss": "Unavailable",
    "ncas": "http://puma.nerc.ac.uk/cylc-review",
    "psc": "Unavailable",
    "Unknown": "Unavailable",
}
HIGHLIGHT_ROSE_ANA_FAILS = [
    "_vs_",
    "lrun_crun_atmos",
    "proc",
    "atmos_omp",
    "atmos_nruncrun",
    "atmos_thread",
]
COMMON_GROUPS = {
    "meto": [
        "all",
        "nightly",
        "developer",
        "xc40",
        "ex1a",
        "spice",
        "xc40_nightly",
        "ex1a_nightly",
        "spice_nightly",
        "xc40_developer",
        "ex1a_developer",
        "spice_developer",
        "ukca",
        "recon",
        "jules",
        "xc40_ukca",
        "ex1a_ukca",
        "spice_ukca",
        "xc40_jules",
        "ex1a_jules",
        "spice_jules",
    ],
    "ecmwf": [],
    "nci": [],
    "bom": [],
    "uoe": [],
    "niwa": [],
    "kma": [],
    "vm": [],
    "jasmin": [],
    "cehwl1": [],
    "mss": [],
    "ncas": [],
    "psc": [],
    "Unknown": [],
}


def _read_file(filename):
    """Takes filename (str)
    Return contents of a file, as list of strings."""
    if os.path.exists(filename):
        with open(filename, "r") as filehandle:
            lines = filehandle.readlines()
    else:
        print('[ERROR] Unable to find file :\n    "{0:s}"'.format(filename))
        raise IOError(
            '_read_file got invalid filename : "{0:s}"'.format(filename)
        )
    return lines


def _write_file(filename, lines, newline=False):
    """Takes filemname and list of strings and opt newline boolean.
    Writes array to file of given name.
    The optional newline argument adds a newline at the end of each
    element of the list.
    Returns None"""
    retn = "\n" if newline else ""
    with open(filename, "w") as filehandle:
        for line in lines:
            filehandle.write("{0:s}{1:s}".format(line, retn))


def _run_command(command, ignore_fail=False):
    """Takes command and command line options as a list.
    Runs the command with subprocess.Popen.
    Returns the exit code, standard out and standard error as list.
    """
    pobj = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    pobj.wait()
    retcode, stdout, stderr = (
        pobj.returncode,
        pobj.stdout.read().decode("utf-8"),
        pobj.stderr.read().decode("utf-8"),
    )
    if retcode != 0 and not ignore_fail:
        print("[ERROR] running {0:s}".format(command))
        print("[INFO] RC: {0:}".format(retcode))
        print("[INFO] Stdout: {0:s}".format(stdout))
        print("[INFO] Stderr: {0:s}".format(stderr))
        raise IOError("run_command")
    # Reformat stdout into a list
    stdout = "".join(stdout)
    stdout = stdout.split("\n")

    return retcode, stdout, stderr


def _remove_quotes(string):
    """Takes, modifies and returns string.
    Removes all quotes from the string.
    None input results in None output"""
    if string is not None:
        string = re.sub(r'"', r"", string)
        string = re.sub(r"'", r"", string)
    return string


def _dict_merge(main_dict, addon_dict, force=False):
    """Merge addon dictionary into main dictionary.
    Takes main_dict, addon_dict and optional bool 'force'
    Returns new merged dictionary.
    Optional argument force=True allows forced overwrite of existing
    value with None from the addon dictionary. Otherwise original
    value is preserved when value in addon dict is None.
    This preserving behaviour differentiates it from main.update(addon)"""
    merged_dict = main_dict.copy()
    for key, value in addon_dict.items():
        if isinstance(value, dict):
            if key not in merged_dict:
                merged_dict[key] = {}
            merged_dict[key] = _dict_merge(merged_dict[key], value)
        else:
            # Switch to Force main to take whatever addon has
            # No matching key in main - take whatever addon has including None
            # Or
            # Override main with contents of addon
            if force or key not in merged_dict or value is not None:
                merged_dict[key] = value
    return merged_dict


def _select_preferred(option_list):
    """Takes a list of strings, returns the fist one that is not None.
    If the strings are report text in preffered order it essentially
    ensures you get the preffered option from a list of choices."""
    pref_opt = None
    for choice in option_list:
        if choice is not None:
            pref_opt = choice
            break
    return pref_opt


def _escape_svn(url):
    """Takes and returns url as string.
    Escape 'svn:' urls as Trac tries to convert them to links."""
    if not re.search(r"!svn://", url):  # Make sure it's not already escaped.
        url = re.sub(r"svn://", r"!svn://", url)
    return url


def _get_current_head_revision(mirror_url, fcm_exec):
    """Given a mirror repository (local) url, uses fcm branch-info to
    retrieve and append the head revision number.
    Requires url and fcm exec path (strings)
    Returns revision number as string"""
    revision = ""
    _, stdout, _ = _run_command([fcm_exec, "branch-info", mirror_url])
    find_last_changed_rev = re.compile(r"Last Changed Rev:\s*(\d+)")
    for line in stdout:
        result = find_last_changed_rev.search(line)
        if result:
            revision = str(result.group(1))
            break
    return revision


def _url_to_trac_link(url):
    """Takes a URL as string, edits text to resemble a Trac link for code
    on the SRS.
    Returns Trac link form of URL or None if 'svn' was absent from the url.
    """
    if re.search(r"/svn/", url):
        link_2_url = re.sub(r"svn", r"trac", url)
        elements = link_2_url.split("/")
        elements.insert(elements.index("trac") + 2, "browser")
        link_2_url = "/".join(elements)
        link_2_url = re.sub(r"@", r"?rev=", link_2_url)
    else:
        link_2_url = None
    return link_2_url


def _parse_string(
    varname,
    lines,
    remove_quotes=True,
    split_on_comma=False,
    default_unknown=False,
):
    """Given a variable name in the rose-suite-run.conf file, return its
    value."""
    find_var = re.compile(r"{0}\s*=\s*(.*)".format(varname))
    if split_on_comma:
        value = [None]
    elif default_unknown:
        value = "Unknown"
    else:
        value = None
    for line in lines:
        result = find_var.search(line)
        if result:
            value = result.group(1)
            if remove_quotes:
                value = _remove_quotes(value.rstrip())
            if split_on_comma:
                # Remove brackets and split on comma
                value = re.sub(r"\[", r"", value)
                value = re.sub(r"\]", r"", value)
                value = value.split(",")
    return value


class SuiteReport(object):
    """Object to hold data and methods required to produce a suite report
    from a rose-stem suite output."""

    def __init__(
        self,
        suite_path,
        log_path=None,
        verbosity=DEFAULT_VERBOSITY,
        sort_by_name=False,
    ):
        """Requires a path to the suite output directory.
        Takes optional arguments for log_path (output dir), and
        verbosity which dictates what tasks are omitted from the log.
        also the boolean sort_by_name to force sorting by task name over
        status when generating the task table in the report.
        """
        suite_path = re.sub(r"/$", "", suite_path)  # remove trailing '/'
        self.suite_path = suite_path
        self.log_path = log_path
        self.sort_by_name = sort_by_name
        self.verbosity = verbosity
        self.creation_time = time.strftime("%Y/%m/%d %X")
        self.uncommitted_changes = 0
        self.suitename = os.path.basename(self.suite_path)

        self.site = "Unknown"
        self.rose_orig_host = None
        self.suite_owner = None
        self.groups = []
        self.job_sources = {}
        self.primary_project = ""
        self.projects = {}
        self.status_counts = defaultdict(int)
        self.status_counts["failed"] = 0

        self.ascertain_suite_owner()
        self.parse_rose_suite_run()
        self.initialise_projects()
        self.parse_suite_rc_processed()
        projects = self.check_versions_files()
        self.job_sources = _dict_merge(self.job_sources, projects)

        # Work out which project this suite is run as - if UM in project sources
        # implies it is the UM. Otherwise it will be whichever is present of
        # Jules and UKCA. Assumes we only run suite_report for UM, Jules, UKCA
        # Logic will require extending once running as part of lfric_apps
        if "UM" in self.job_sources.keys():
            self.primary_project = "UM"
        else:
            if "JULES" in self.job_sources.keys():
                self.primary_project = "JULES"
            elif "UKCA" in self.job_sources.keys():
                self.primary_project = "UKCA"

        self.groups = [_remove_quotes(group) for group in self.groups]

        fcm_exec = FCM[self.site]
        invalid = []
        for project in self.job_sources:
            proj_dict = self.job_sources[project]
            if "repo loc" in proj_dict:
                proj_dict["repo loc"] = self.convert_to_srs(
                    proj_dict["repo loc"], self.projects
                )
            else:
                proj_dict["repo loc"] = self.convert_to_srs(
                    proj_dict["tested source"], self.projects
                )

            proj_dict["repo mirror"] = self.convert_to_mirror(
                proj_dict["repo loc"], self.projects
            )
            # If the mirror doesn't exist, move on to the next project.
            if not self.check_repository(fcm_exec, proj_dict["repo mirror"]):
                invalid.append(project)
                continue

            proj_dict["parent mirror"] = self.set_parent(
                fcm_exec, proj_dict["repo mirror"]
            )

            proj_dict["parent loc"] = self.convert_to_srs(
                proj_dict["parent mirror"], self.projects
            )
            # Check "repo loc" and "parent loc" have revisions,
            # and if not, try to get a head of 'trunk' one for them.
            for location in ("repo", "parent"):
                url = proj_dict[location + " loc"]
                mirror_url = proj_dict[location + " mirror"]
                if url is None or mirror_url is None:
                    continue
                if ":" in url and "@" not in url:
                    revision = _get_current_head_revision(mirror_url, fcm_exec)
                    proj_dict[location + " loc"] = url + "@" + revision
                    proj_dict[location + " mirror"] = (
                        mirror_url + "@" + revision
                    )
            proj_dict["repo link"] = self.generate_link(proj_dict["repo loc"])
            proj_dict["parent link"] = self.generate_link(
                proj_dict["parent loc"]
            )
            # If those attempts to generate links didn't work, try the hope
            # and guess approach.
            if proj_dict["repo link"] is None:
                proj_dict["repo link"] = self.link_from_loc_layout(
                    proj_dict["repo link"], proj_dict["repo mirror"], fcm_exec
                )
            if proj_dict["parent link"] is None:
                proj_dict["parent link"] = self.link_from_loc_layout(
                    proj_dict["parent loc"],
                    proj_dict["parent mirror"],
                    fcm_exec,
                )
            # Final attempt to ensure the links have revision numbers and not
            # keywords which aren't evaluated in the browser.
            if proj_dict["repo link"] is not None and re.search(
                r"rev=[a-zA-z]", proj_dict["repo link"]
            ):
                revision = self.revision_from_loc_layout(
                    proj_dict["repo mirror"], fcm_exec
                )
                proj_dict["repo link"] = re.sub(
                    r"rev=[a-zA-z0-9.]+",
                    "rev=" + revision,
                    proj_dict["repo link"],
                )
            proj_dict["human repo loc"] = self.convert_to_keyword(
                proj_dict["repo loc"], self.projects
            )
            proj_dict["human parent"] = self.convert_to_keyword(
                proj_dict["parent loc"], self.projects
            )
            proj_dict["ticket no"] = self.ascertain_ticket_number(
                proj_dict["repo mirror"], fcm_exec
            )
            proj_dict["bdiff_files"] = self.get_altered_files_list(
                proj_dict["repo mirror"]
            )

        # Check to see if ALL the groups being run fall into the "common groups"
        # category. This is used to control automatic hiding of successful tasks
        # later.
        self.only_common_groups = all(
            [group.strip() in COMMON_GROUPS[self.site] for group in self.groups]
        )

        # Finally, remove any projects which were deemed invalid.
        for project in invalid:
            del self.job_sources[project]

    def debug_print_obj(self):
        """Debug print method.
        Prints everything in the SuiteReport object."""
        print("-" * 80 + "\nSet up SuiteReport object\n" + "-" * 80 + "\n\n")
        for key, value in self.__dict__.items():
            if key == "projects":
                print('{0:s} contains "{1:d}" entries.'.format(key, len(value)))
            elif key == "sort_by_name":
                if value:
                    print('{0:s} is :"True"'.format(key))
                else:
                    print('{0:s} is :"False"'.format(key))
            elif key == "only_common_groups":
                if value:
                    print('{0:s} is :"True"'.format(key))
                else:
                    print('{0:s} is :"False"'.format(key))
            elif key == "verbosity":
                text = "Verbosity level is set to : "
                if value >= 4:
                    print(
                        text
                        + "Hide Housekeeping, Gatekeeping and Successful tasks"
                    )
                elif value >= 3:
                    print(
                        text
                        + "Hide Housekeeping, Gatekeeping and if all groups run"
                        'were "common" groups also hide Successful tasks'
                    )
                elif value >= 2:
                    print(text + "Hide Housekeeping and Gatekeeping tasks")
                elif value >= 1:
                    print(text + "Hide Housekeeping tasks")
                else:
                    print(text + "Forcibly Print Everything.")
            elif key == "job_sources":
                self.print_job_sources(value)
            else:
                print('{0:s} is :"{1:}"'.format(key, value))
        print(
            "\n" + "-" * 80 + "\nEnd of SuiteReport object\n" + "-" * 80 + "\n"
        )

    @staticmethod
    def print_job_sources(job_srcs_dict):
        """Debug print method.
        Prints everything in projects dictionary."""
        for key, value in job_srcs_dict.items():
            print("    {0:s} :".format(key))
            for sub_key, sub_value in value.items():
                if isinstance(sub_value, bool):
                    if sub_value:
                        print('        {0:s} is :"True"'.format(sub_key))
                    else:
                        print('        {0:s} is :"False"'.format(sub_key))
                else:
                    print('        {0:s} is :"{1:}"'.format(sub_key, sub_value))

    def parse_suite_rc_processed(self):
        """Parse the suite.rc.processed file.
        Extract all projects present that begin with a "SOURCE_".
        Allow SOURCE_<project> to override any SOURCE_<project>_<extension>
        entries. Creating a dictionary of format {<project> : <URL>,...}
        Also Extract the host machine rose was launched on.
        Takes full path for suite dir.
        Sets class variables"""

        suite_dir = self.suite_path

        rose_orig_host = "Unknown rose_orig_host"
        srp_file = os.path.join(suite_dir, PROCESSED_SUITE_RC)
        find_orig_host = re.compile(r"ROSE_ORIG_HOST\s*=\s*(.*)")
        # in pattern below, need to include "_REV" after the project name and
        # before the " *=" and then exclude lines with "_REV" later as
        # otherwise the search will identify PROJ_REV as a unique project
        # name. The other option would be to have an alternate 3rd group match
        # of "_.*?" but that would exclude any future project names that might
        # have an underscore in them.
        find_sources = re.compile(
            r"\s*(?:HOST_SOURCE_|SOURCE_)(.*?)(|_BASE|_MIRROR|_REV)\s*=\s*(.*)"
        )
        sources = {}
        multiple_branches = {}
        for line in _read_file(srp_file):
            # check for ROSE_ORIG_HOST
            result = find_orig_host.search(line)
            if result:
                rose_orig_host = result.group(1).rstrip()
            # check for SOURCE_.*
            result = find_sources.match(line)
            # Discard the ones which were SOURCE_PROJ_REV
            if result and result.group(2) != "_REV":
                # Allow SOURCE_PROJ to override any existing entries
                # Otherwise only add new entries
                if result.group(1) not in sources or result.group(2) == "":
                    sources[result.group(1)] = {}
                    if " " in result.group(3):
                        multiple_branches[(result.group(1))] = result.group(3)
                        sources[result.group(1)][
                            "tested source"
                        ] = result.group(3).split()[0]
                    else:
                        sources[result.group(1)][
                            "tested source"
                        ] = result.group(3)

        self.rose_orig_host = rose_orig_host
        self.job_sources = sources
        self.multi_branches = multiple_branches
        return

    def ascertain_suite_owner(self):
        """Find suite owner.
        To retrieve suite owner :
           Take default from ENV['USER'] or 'Unknown Suite Owner' if none.
           Use ENV['CYLC_SUITE_OWNER'] if present.
        Sets Class variables"""

        # Get a default
        suite_owner = os.environ.get("USER", "Unknown Suite Owner")
        suite_owner = os.environ.get("CYLC_SUITE_OWNER", suite_owner)
        self.suite_owner = suite_owner
        return

    def parse_rose_suite_run(self):
        """Parse rose-suite-run.conf file.
        Takes full path for suite dir.
        Sets class variables"""

        suite_dir = self.suite_path
        rsr_file = os.path.join(suite_dir, "log", ROSE_SUITE_RUN_CONF)
        lines = _read_file(rsr_file)
        self.site = _parse_string("SITE", lines, default_unknown=True)
        self.groups = _parse_string(
            "RUN_NAMES", lines, split_on_comma=True, remove_quotes=False
        )
        self.fcm = _parse_string("FCM_VERSION", lines)
        self.cylc = _parse_string("CYLC_VERSION", lines)
        self.rose = _parse_string("ROSE_VERSION", lines)

        # This test is a little problematic when running this script on a JULES
        # rose-stem suite as JULES has no 'need' of the two compare variables
        # and to prevent the warning their absence would produce from occuring
        # unnecessarily in JULES they have been added to rose-suite.conf for now
        compare_output = _parse_string("COMPARE_OUTPUT", lines)
        compare_wallclock = _parse_string("COMPARE_WALLCLOCK", lines)
        self.required_comparisons = (
            compare_output == "true" and compare_wallclock == "true"
        )

        self.trustzone = None
        if "TRUSTZONE" in os.environ:
            self.trustzone = os.environ["TRUSTZONE"]

        self.host_xcs = False
        if self.site == "meto":
            for line in lines:
                if "HOST_XC40='xcsr'" in line:
                    self.host_xcs = True
        return

    def initialise_projects(self):
        """Uses fcm kp to initialise a directory containing project keywords
        linked to SVN URLS. Format {<project> : <URL>,...}
        Takes full path for suite dir.
        Sets class variable"""

        fcm_exec = FCM[self.site]

        projects = {}
        _, stdout, _ = _run_command([fcm_exec, "kp"])
        find_primary_loc = re.compile(r"location{primary}")
        find_projects = re.compile(r"\[(.*)\]\s*=\s*(.*)")
        find_x_keyword = re.compile(r".x$")
        find_xm_keyword = re.compile(r".xm$")
        find_srs_url = re.compile(r"https://code.metoffice")
        find_mirror_url = re.compile(r"svn:|https://")
        for line in stdout:
            if not find_primary_loc.search(line):
                continue
            result = find_projects.search(line)
            if result:
                project = result.group(1)
                url = result.group(2)
                # Check for keywords conforming to the meto prescribed pattern
                # of ending in '.x' for the external repo and '.xm' for the
                # local mirror.
                if (
                    find_x_keyword.search(project) and find_srs_url.match(url)
                ) or (
                    find_xm_keyword.search(project)
                    and find_mirror_url.match(url)
                ):
                    projects[project] = url
        self.projects = projects
        return

    def check_versions_files(self):
        """Locate the log/*.version files.
        Call parse_versions_fileto parse the contents of each file.
        Recover which projects are being augmented by branch or WC
        Takes full path for suite dir.
        Returns dictionary of project dictionares and number of projects
        with uncommitted changes."""
        projects = {}
        self.uncommitted_changes = 0
        find_proj_name = re.compile(r"/(\w+)-\d+.version")

        version_files = glob.glob(
            "{0:s}/*.version".format(os.path.join(self.suite_path, "log"))
        )

        for vfile in version_files:
            if "rose-suite-run.version" in vfile:
                continue
            result = find_proj_name.search(vfile)
            if result:
                project = result.group(1).upper()
                projects[project] = {}
                (url, revision, wc_changes) = self.parse_versions_file(vfile)
                projects[project]["last changed rev"] = revision
                projects[project]["working copy changes"] = wc_changes
                projects[project]["version file"] = os.path.basename(vfile)
                if wc_changes:
                    self.uncommitted_changes += 1
                if url is not None:
                    if revision is not None:
                        ending = "@" + revision
                    else:
                        ending = ""
                    projects[project]["repo loc"] = url + ending
        return projects

    @staticmethod
    def parse_versions_file(vfile):
        """Parse a versions file to extract the url and revision for
        the branches behind any working copies, plus any uncommitted
        changes.
        Takes full path to a .version file.
        Returns url and revision as strings plus wc changes as boolean."""
        url = None
        revision = None
        working_copy_changes = False
        find_svn_status = re.compile(r"SVN STATUS", re.IGNORECASE)
        find_url = re.compile(r"URL:\s*")
        find_last_changed_rev = re.compile(r"Last Changed Rev:\s*")
        for line in _read_file(vfile):
            if find_svn_status.search(line):
                working_copy_changes = True
            if find_url.match(line):
                url = find_url.sub(r"", line).rstrip()
            if find_last_changed_rev.match(line):
                revision = find_last_changed_rev.sub(r"", line).rstrip()
        return url, revision, working_copy_changes

    @staticmethod
    def set_parent(fcm_exec, mirror_url):
        """For given URL, on the internal mirror repository, use
        'fcm branch-info' to try and ascertain the branch parent, if any.
        Takes fcm_exec path and mirror_url as strings.
        Returns parent URL or None"""

        parent = None
        stdout = ""
        command = [fcm_exec, "branch-info", mirror_url]
        _, stdout, _ = _run_command(command, ignore_fail=True)
        find_branch_parent = re.compile(r"Branch Parent:\s*(.*)")
        for line in stdout:
            result = find_branch_parent.search(line)
            if result:
                parent = result.group(1).rstrip()
        return parent

    @staticmethod
    def check_repository(fcm_exec, url):
        """Checks whether a given repository is accessible or not.
        Takes fcm_exec path and a url (either SRS or mirror) as strings.
        Returns True if the repository exists, False otherwise."""
        retcode = 0
        command = [fcm_exec, "info", url]
        retcode, _, _ = _run_command(command, ignore_fail=True)
        if retcode == 0:
            return True
        return False

    @staticmethod
    def export_file(repo_url, fname, outname="~/temp.txt"):
        """
        Runs an fcm export on a file and saves it as outname
        Attempts to check it out 5 times to account for any network glitches.
        Returns None if all attempts fail, otherwise the user expanded path to
        the file
        Inputs: repo_url, eg. fcm:um.xm_tr
                fname: the path of the file in the repo
                outname: the path to the output file. Default ~/temp.txt
        """

        fname = fname.lstrip("/")
        outname = os.path.expanduser(outname)

        # Try 5 times, if all fail then use working copy version
        for _ in range(5):
            try:
                subproc = "fcm export -q {}/{} {} --force".format(
                    repo_url, fname, outname
                )
                subprocess.check_output(subproc, shell=True)
                return outname
            except subprocess.CalledProcessError as error:
                print(error)
        else:
            return None

    @staticmethod
    def clean_tempfile(fname="~/temp.txt"):
        """
        Clean up a temp file exported
        """

        try:
            os.remove(os.path.expanduser(fname))
        except EnvironmentError:
            pass

    def generate_owner_dictionary(self, mode):
        """
        Function that parses an owners file to create a dictionary of owners,
        mapping a configuration/section to an owner

        Input:
            mode - either config or code depending on dictionary created
        """

        if mode == "config":
            fname = "ConfigOwners.txt"
            sep = "Configuration"
        elif mode == "code":
            fname = "CodeOwners.txt"
            sep = "Owner"
        else:
            return None

        # Export the Owners file from the HOT
        exported_file = "~/tmp_owners.txt"
        file_path = self.export_file("fcm:um.xm_tr", fname, exported_file)
        if file_path is None:
            # Couldn't check out file - use working copy Owners file instead
            wc_path = get_working_copy_path(
                self.job_sources["UM"]["tested source"]
            )
            if not wc_path:
                wc_path = ""
            file_path = os.path.join(wc_path, fname)
            print("Using the checked out version of Owners file")

        # Read through file and generate dictionary
        try:
            with open(file_path, "r") as inp_file:
                owners_dict = {}
                inside_listing = False
                for line in inp_file:
                    if "{{{" in line:
                        inside_listing = True
                        continue
                    if "}}}" in line:
                        inside_listing = False
                        continue
                    if inside_listing:
                        if line != "\n" and sep not in line:
                            dummy_list = line.split()
                            section = dummy_list[0].strip()
                            owners = dummy_list[1].strip()
                            if "umsysteam" in owners:
                                owners = "!umsysteam@metoffice.gov.uk"
                            try:
                                others = dummy_list[2].replace("\n", "")
                                if others == "--":
                                    others = ""
                            except IndexError:
                                others = ""
                            owners_dict.update(
                                {section.lower(): [owners, others]}
                            )
        except EnvironmentError:
            print("Can't find working copy for Owners File")
            return None

        # Clean up the checked out copy of the owners file
        self.clean_tempfile(exported_file)

        return owners_dict

    def create_approval_table(self, needed_approvals, mode):
        """
        Function to write out the trac.log table for config and CO approvals
        Input: needed_approvals - dictionary with keys as owners and values,
                                  a list of configs or code sections
               mode - either "config" or "code" depending on which type of table
                      is being created
        """

        table = ["'''Required " + mode.capitalize() + " Owner Approvals'''"]

        if mode == "config":
            table += [" || '''Owner''' || '''Approval''' || '''Configs''' || "]
        else:
            table += [
                " || '''Owner (Deputy)''' || '''Approval''' || "
                + "'''Code Section''' || "
            ]

        if needed_approvals is None:
            table += [
                " |||||| No UM "
                + mode.capitalize()
                + " Owner Approvals Required || "
            ]
        else:
            for owner in needed_approvals.keys():
                row = " || " + owner + " || Pending || "
                count = 0
                for val in needed_approvals[owner]:
                    if count == 3:
                        row += "[[br]]"
                        count = 0
                    row += "{{{" + val + "}}} "
                    count += 1
                row += " || "
                table.append(row)

        table.append("")

        return table

    def get_config_owners(self, failed_configs, config_owners):
        """
         Function that reads through a list of failed rose-ana jobs and records
         owners for each job thathas failed.

        Input:
                 failed_configs - list of failed rose-ana tasks
                 config_owners - dictionary mapping config owners to configs
        """

        # Dictionary to store needed approvals
        needed_approvals = defaultdict(set)

        for job in failed_configs:
            job = job.lower()

            # Split the job name into it's various 'parts'
            parts = job.split("-")

            # Get the Config Name
            # try-except as mule rose-ana names follow different format
            try:
                config = parts[2]
            except IndexError:
                if "mule" in job:
                    config = "mule"
                else:
                    config = ""

            # Get the config owner + others to notify
            try:
                owners = config_owners[config]
            except KeyError:
                owners = ["Unknown", ""]
            owner = owners[0]
            notify = owners[1]

            # If others to notify, append names to config name
            if notify != "":
                config += "(" + notify + ")"

            # Record Owner and Config
            needed_approvals[owner].add(config)

        return needed_approvals

    def required_config_approvals(self, failed_configs):
        """
        Calls functions to create a table of required config approvals
        """

        config_owners = self.generate_owner_dictionary("config")
        if config_owners is None:
            return None

        config_approvals = self.get_config_owners(failed_configs, config_owners)

        if len(config_approvals.keys()) == 0:
            config_approvals = None

        approval_table = self.create_approval_table(config_approvals, "config")

        return approval_table

    def get_code_owners(self, code_owners):
        """
        Function to get required code owner approvals based on fcm_bdiff
        - code_owners - dict returning code owners for a given code section
        """

        # Get list of altered files and exit if no files changed
        # 'UM' used here and just below as this function is currently only valid
        # for the UM. Hopefully lfric_apps will be able to use similar in the
        # future - at this point we can change 'UM' to self.primary_project
        bdiff_files = self.job_sources["UM"]["bdiff_files"]
        if len(bdiff_files) == 0:
            return None

        # Get the mirror repo and remove the @REVISION part
        repo_loc = self.job_sources["UM"]["repo mirror"].split("@")[0]

        # Dictionary to store needed approvals
        needed_approvals = defaultdict(set)

        # Get Owners for each file changed
        for fle in bdiff_files:
            if ".." in fle:
                # This is to fix an invalid path returned by fcm_bdiff in the
                # case a branch has been reversed off trunk (see comments in
                # get_branch_diff_filenames() for detail)
                # The file path (fle) is split by the first example of the
                # branch_name and then the file path as we expect is the last
                # value of that list. We then remove any trailing '/'.
                branch_name = repo_loc.split("/")[-1]
                fle = fle.split(branch_name, 1)[1].strip("/")
            fpath = fle
            fle = fle.lower()

            # Manually sort directories with known sections
            if "configowners.txt" in fle or "codeowners.txt" in fle:
                continue
            if fle.startswith("admin"):
                needed_approvals["!umsysteam@metoffice.gov.uk"].add("admin")
                continue
            if fle.startswith("bin"):
                needed_approvals["!umsysteam@metoffice.gov.uk"].add("bin")
                continue
            if fle.startswith("fcm-make"):
                section = "fcm-make_um"
            elif fle.startswith("fab"):
                section = "fab"
            elif fle.startswith("rose-stem"):
                if "umdp3_check" in fle:
                    section = "umdp3_checker"
                elif "run_cppcheck" in fle:
                    section = "run_cppcheck"
                elif "rose-stem/bin" in fle:
                    section = "rose_bin"
                else:
                    section = "rose_stem"
            elif fle.startswith("rose-meta"):
                if "versions.py" in fle:
                    section = "upgrade_macros"
                elif "rose-meta.conf" in fle:
                    section = "rose-meta.conf"
                else:
                    section = "stash"
            else:
                # Find area of files in other directories
                file_path = self.export_file("fcm:um.xm_tr", fpath)
                if file_path is None:
                    file_path = ""

                try:
                    with open(file_path, "r") as inp_file:
                        for line in inp_file:
                            if "file belongs in" in line:
                                section = line.strip("\n")
                                break
                        else:
                            section = ""
                except EnvironmentError:
                    section = ""

                # Clean up checked out file
                self.clean_tempfile()

                # Get code area name out
                section = re.sub(r"/\*", "", section)
                section = re.sub(r"\*/", "", section)
                try:
                    section = section.split(":")[1].strip().lower()
                except IndexError:
                    section = ""

            # Compare area name to code owners list
            try:
                owners = code_owners[section]
                owner, deputy = owners
                if len(deputy) > 0:
                    owner += " (" + deputy + ")"
                needed_approvals[owner].add(section)
            except KeyError:
                needed_approvals[
                    "Unknown - ensure section " "is in CodeOwners.txt"
                ].add(section)

        return needed_approvals

    def required_co_approvals(self):
        """
        Calls functions to create table of code owner approvals
        """

        code_owners = self.generate_owner_dictionary("code")
        if code_owners is None:
            return None

        code_approvals = self.get_code_owners(code_owners)
        if code_approvals is False:
            return None

        approval_table = self.create_approval_table(code_approvals, "code")

        return approval_table

    @staticmethod
    def parse_lfric_extract_list(fpath="~/temp.txt"):
        """
        Read through the lfric_extract list and get a list of files and dirs.
        Return a dictionary with keys 'files' and 'dirs'
        """

        files = []
        dirs = []
        in_include_section = False

        with open(os.path.expanduser(fpath)) as input_file:
            for line in input_file:
                if in_include_section:
                    item = line.rstrip("\\").strip()
                    if "." in item.split("/")[-1]:
                        files.append(item)
                    else:
                        dirs.append(item)
                    if not line.endswith("\\"):
                        in_include_section = False
                if "extract.path-incl" in line:
                    in_include_section = True

        return {"files": files, "dirs": dirs}

    def get_lfric_interactions(self, extract_list):
        """
        Function to count the number of project sources with modified files
        extracted by lfric. Takes in a dict with keys files and dirs and values
        are those extracted by lfric
        """

        num_interactions = 0

        for project in self.job_sources.keys():
            for mod_file in self.job_sources[project]["bdiff_files"]:
                if "trunk" in mod_file:
                    mod_file = mod_file.split("trunk/")[-1]
                # Check modified file isn't in extracted files list
                matching_item = False
                if mod_file in extract_list["files"]:
                    matching_item = True
                # Loop over directories extracted
                # Check that the directory doesn't contain the modified file
                for drc in extract_list["dirs"]:
                    if drc in mod_file:
                        matching_item = True
                        break
                if matching_item:
                    num_interactions += 1
                    break

        return num_interactions

    @staticmethod
    def write_lfric_testing_message(num_interactions):
        """
        Based on no. projects with lfric interaction write a message stating
        lfric testing requirements
        """

        message = []

        if num_interactions > 0:
            if num_interactions > 1:
                message += ["There were {} projects ".format(num_interactions)]
            else:
                message += ["There was 1 project "]
            message += [
                "with LFRic interaction.[[br]]LFRic testing is "
                + "'''required''' before this ticket is submitted for review."
            ]
        else:
            message += [
                "No files shared with LFRic have been modified.[[br]]LFRic "
                + "testing is not required for this ticket."
            ]

        message.append("")
        return message

    def check_lfric_extract_list(self):
        """
        Determine whether any files modified in source branches are extracted by
        lfric.
        Return a trac formatted string stating whether LFRic testing is required
        """

        return_message = ["'''LFRic Testing Requirements'''"]
        return_message.append("")

        # Export the extract list from the lfric trunk
        exported_extract_file = "~/tmp_extract.txt"
        extract_list_path = self.export_file(
            "fcm:lfric.xm_tr",
            "lfric_atm/fcm-make/extract.cfg",
            exported_extract_file,
        )

        try:
            extract_list_dict = self.parse_lfric_extract_list(extract_list_path)
        except EnvironmentError or TypeError or AttributeError:
            # Potential error here changed type between python2 and 3
            extract_list_path = None

        # If the path returned is None, the extract list failed, most likely as
        # the user doesn't have lfric access. In this case return a warning.
        if extract_list_path is None:
            return_message += [
                "Unable to export the lfric extract_list. "
                + "LFRic testing may be required."
            ]
            return return_message

        # Clean up the checked out copy of the extract list file
        self.clean_tempfile(exported_extract_file)

        num_interactions = self.get_lfric_interactions(extract_list_dict)

        return_message += self.write_lfric_testing_message(num_interactions)

        return return_message

    def generate_task_table(
        self,
        data,
        common_groups=False,
        verbosity=DEFAULT_VERBOSITY,
        sort_by_name=False,
        status_counts=defaultdict(int),
    ):
        """Returns a trac-formatted table of the tasks run in this suite.
        Tasks are provided in a dictionary of format {"task" : "status",...}
        verbosity (int) sets verbosity level. In practice, as the number
        increases, verbosity decreases, and the range of tasks hidden from
        the report increases.
        sort_by_name (bool) sorts by task name when true, otherwise sorting
        is done by status.
        Potentially 2 summary tables are also produced. One details the
        number of tasks found with each status type. The 2nd is only present
        if not empty and indicates how many tasks of the relevant types have
        been hidden"""

        def key_by_name_or_status(task_item):
            """A key generating function for use by sorted.
            task_item is a tuple of (name, status).
            If sorting by name, return a tuple of (name, status),
            otherwise return (status, name) for use as the sorting key"""
            if sort_by_name:
                return task_item
            return (task_item[1], task_item[0])

        def forced_status_sort(item_tuple):
            """A key generating function for use by sorted.
            item_tuple is a tuple of the ("key", "value") pair from a
            dictionary.
            If the 'key' is in the DESIRED_ORDER list then return the key's
            index from that list, otherwise return the key.
            This forces the items in the DESIRED order list to be listed
            first, in the order they appear in the list, followed by all the
            other keys.
            caveat, it relies on numbers preceeding alphabetic characters
            and all the status keys starting with a letter."""
            key = item_tuple[0]
            if key in DESIRED_ORDER:
                return str(DESIRED_ORDER.index(key))
            return key

        hidden_counts = defaultdict(int)
        lines = [" || '''Task''' || '''State''' || "]
        find_housekeep = re.compile(r"housekeep")
        find_monitor = re.compile(r"monitor")
        find_gatekeeper = re.compile(r"gatekeeper")
        failed_configs = []
        for task, state in sorted(
            list(data.items()), key=key_by_name_or_status
        ):
            # Count the number of times task have any given status.
            status_counts[state] += 1
            if (verbosity >= 1) and find_housekeep.match(task):
                hidden_counts["''Housekeeping''"] += 1
                continue
            if (verbosity >= 2) and find_gatekeeper.match(task):
                hidden_counts["''Gatekeeping''"] += 1
                continue
            if find_monitor.match(task):
                hidden_counts["''Monitoring''"] += 1
                continue
            highlight_start = "'''"
            highlight_end = "'''"
            if "succeeded" in state:
                # Omit printing the task if verbosity level is set to omit all
                # successful tasks, or omit them only when all groups run are
                # "common groups"
                highlight_start = ""
                highlight_end = ""
                if (verbosity >= 4) or (verbosity >= 3 and common_groups):
                    hidden_counts["'''Succeeded'''"] += 1
                    continue
            elif "rose_ana" in task and "failed" in state:
                # Check if task requires extra care
                for extra_care_string in HIGHLIGHT_ROSE_ANA_FAILS:
                    if extra_care_string in task:
                        highlight_start = (
                            "'''[[span(style=color: #FF00FF, *****"
                        )
                        highlight_end = "***** )]]'''"
                        status_counts[PINK_FAIL_TEXT] += 1
                        status_counts[state] -= 1
                        break
                else:
                    # Record this as a failed config
                    failed_configs.append(task)

            lines.append(
                " || {0:s} || {2:s}{1:s}{3:s} || ".format(
                    task, state, highlight_start, highlight_end
                )
            )
        if len(lines) == 1:
            lines.append(
                " |||| This table is deliberately empty as all tasks "
                "are hidden || "
            )

        status_summary = ["'''Suite Output'''"]
        status_summary += self.gen_resources_table()
        status_summary += [""]
        status_summary.append(" |||| '''All Tasks''' || ")
        status_summary.append(" || '''Status''' || '''No. of Tasks''' || ")
        for status, count in sorted(
            status_counts.items(), key=forced_status_sort
        ):
            status_summary.append(
                " || {0:s} || {1:6d} || ".format(status, count)
            )
        status_summary.append("")
        if len(hidden_counts) > 0:
            status_summary.append(" |||| '''Hidden Tasks''' || ")
            status_summary.append(
                " || '''Type''' || '''No. of Tasks Hidden''' || "
            )
            for task_type, count in hidden_counts.items():
                status_summary.append(
                    " || {0:s} || {1:6d} || ".format(task_type, count)
                )
            status_summary.append("")

        # Check whether lfric shared files have been touched
        # Not needed if lfric the suite source
        if "LFRIC" not in self.primary_project:
            lfric_testing_message = self.check_lfric_extract_list()

        # Generate table for required config and code owners
        # Only run if UM is the primary project
        return_list = []
        if self.primary_project == "UM":
            co_approval_table = self.required_co_approvals()
            if co_approval_table:
                return_list += co_approval_table
            config_approval_table = self.required_config_approvals(
                failed_configs
            )
            if config_approval_table:
                return_list += config_approval_table
        return lfric_testing_message + return_list + status_summary + lines

    @staticmethod
    def convert_to_mirror(url, projects_dict):
        """Take a URL as a string, and a dictionary of {project : url, ...}
        If url is a shared repository URL in the projects dictionary convert
        to an internal mirror URL if also available.
        Otherwise return the original URL
        Assumes mirror loc of proj with url svn:something/somewhere/project
        is given as svn:something/somewhere/projectm
        """
        if url is None:
            return None
        mirror_url = url
        for proj, proj_url in projects_dict.items():
            # checking given url against urls in the projects_dict
            if proj_url in url:
                new_proj = proj + "m"
                if new_proj in projects_dict:
                    old_proj_url = proj_url
                    new_proj_url = projects_dict[new_proj]
                    mirror_url = re.sub(
                        old_proj_url, new_proj_url, url, count=1
                    )
                    break
            # checking given url against keywords in the projects_dict
            elif proj in url:
                new_proj = proj + "m"
                if new_proj in projects_dict:
                    mirror_url = re.sub(proj, new_proj, url, count=1)
                    break
        return mirror_url

    @staticmethod
    def convert_to_srs(url, projects_dict):
        """Take a URL as a string, and a dictionary of {project : url, ...}
        If url is a mirror repository URL in the projects dictionary convert
        to an SRS URL if also availble.
        Otherwise return the original URL
        """
        if url is None:
            return None
        srs_url = url
        for proj, proj_url in projects_dict.items():
            # Only check for keywords which correspond to mirror or SRS format
            if re.search(r".x(|m)$", proj):
                if re.search(proj_url, url):
                    # checking given url against urls in the projects_dict
                    shared_project = re.sub(r"m$", r"", proj)
                    if shared_project in projects_dict:
                        mirror_url = proj_url
                        shared_url = projects_dict[shared_project]
                        srs_url = re.sub(mirror_url, shared_url, url, count=1)
                        break
                elif re.search("fcm:" + proj + r"[^m]", url):
                    # Looking for an fcm: shorthand notation based on keyword.
                    shared_project = re.sub(r"m$", r"", proj)
                    if shared_project in projects_dict:
                        # if the fcm keyword ends in '_tr' it's on the trunk
                        if re.match(r"fcm:" + proj + r"_tr", url):
                            srs_url = re.sub(
                                r"fcm:" + proj + r"_tr",
                                projects_dict[shared_project] + r"/trunk",
                                url,
                                count=1,
                            )
                        # if the fcm keyword ends in '_br' it's from branches
                        elif re.match(r"fcm:" + proj + r"_br", url):
                            srs_url = re.sub(
                                r"fcm:" + proj + r"_br",
                                projects_dict[shared_project] + r"/branches",
                                url,
                                count=1,
                            )
                        # maintain keyword style, but convert to srs.
                        else:
                            srs_url = re.sub(proj, shared_project, url, count=1)
                        break
        return srs_url

    @staticmethod
    def convert_to_keyword(url, projects):
        """Takes url and project dictionary.
        Convert a the URL to a keyword based version, if a keyword exists
        in the project dictionary provied.
        Returns None if no keyword is defined.
        """
        if url is None:
            return None
        keyword_url = None
        for proj, proj_url in projects.items():
            if proj_url in url:
                old_proj_url = projects[proj]
                new_proj_url = "fcm:{0:s}".format(proj)
                keyword_url = re.sub(old_proj_url, new_proj_url, url, count=1)
                keyword_url = re.sub(r"/trunk", "_tr", keyword_url, count=1)
                keyword_url = re.sub(r"/branches", "_br", keyword_url, count=1)
                break
            if "fcm:" + proj in url:
                keyword_url = url
                break
        return keyword_url

    def generate_link(self, url):
        """Given a URL, see if it can be made into a shared repository link.
        Returns a link as a str or None"""
        link = None
        if url is not None:
            # Look for a matching part of the URL in the list of projects
            for _, svn in self.projects.items():
                if re.search(svn, url):
                    link = _url_to_trac_link(url)
                    break
        return link

    @staticmethod
    def link_from_loc_layout(url, mirror_url, fcm_exec):
        """Attempt to generate a link to a url using a bunch of assumptions
        we know mostly hold due to working practices at the Met Office.
        Takes url, mirror url and fcm exec path as strings.
        Returns Link as string or None"""
        link = None
        if url is None or mirror_url is None or re.search(r"^file:/", url):
            return None
        _, stdout, _ = _run_command([fcm_exec, "loc-layout", mirror_url])
        path = None
        root = None
        lproject = None
        revision = None
        find_path = re.compile(r"^path:\s*")
        find_root = re.compile(r"^root:\s*")
        find_project = re.compile(r"^project:\s*")
        find_peg_rev = re.compile(r"^peg_rev:\s*")
        for line in stdout:
            if find_path.match(line):
                path = find_path.sub(r"", line)
                continue
            if find_root.match(line):
                root = find_root.sub(r"", line)
                continue
            if find_project.match(line):
                lproject = find_project.sub(r"", line)
                continue
            if find_peg_rev.match(line):
                revision = find_peg_rev.sub(r"", line)
        if root is not None and lproject is not None and path is not None:
            # Convert to a trac url.
            if re.search(r"/svn/", url):
                url = re.sub(r"svn", r"trac", url)
                elements = url.split("/")
                elements.insert(elements.index("trac") + 2, "browser")
                url = "/".join(elements)
                if revision is not None:
                    link = url + "?rev={0:s}".format(revision)
                else:
                    link = url
        return link

    @staticmethod
    def revision_from_loc_layout(mirror_url, fcm_exec):
        """Attempt to recover a revision number using a url to the mirror
        repository. Also used to translate vn4.3 into 1234"""
        if mirror_url is None:
            return None
        _, stdout, _ = _run_command([fcm_exec, "loc-layout", mirror_url])
        revision = None
        find_peg_rev = re.compile(r"^peg_rev:\s*")
        for line in stdout:
            if find_peg_rev.match(line):
                revision = find_peg_rev.sub(r"", line)
                break
        return revision

    def generate_project_table(self):
        """Returns a trac-formatted table containing the project source
        trees used in this suite.
        Method of SuiteReport object.
        Returns list of trac formatted table rows."""

        lines = [
            " || '''Project''' || '''Tested Source Tree''' || "
            + "'''Repository Location''' || '''Branch Parent''' || "
            + "'''Ticket number''' || '''Uncommitted Changes''' ||"
        ]

        def gen_table_element(text_list, link, bold=False):
            """Takes list of items (strings or Nones) in preference order.
            Calls _select_preferred to get the first non None entry in list.
            Optional Bool "bold" turns on Trac formatting of bold text.
            Formats text as a Trac link if link is not None.
            Returns 'highlighted' text/link and "||" to form a single
            element of a trac table."""
            text = _select_preferred(text_list)
            highlight = "'''" if bold else ""
            if text is not None and link is not None:
                element = " {2:s}[{0:s} {1:s}]{2:s} || ".format(
                    link, text, highlight
                )
            elif text is not None:
                element = " {1:s}{0:s}{1:s} || ".format(
                    _escape_svn(text), highlight
                )
            else:
                element = " || "
            return element

        for project in sorted(self.job_sources):
            proj_dict = self.job_sources[project]

            line = " || {0:s} || ".format(project)
            line += gen_table_element([proj_dict["tested source"]], None)
            line += gen_table_element(
                [proj_dict["human repo loc"], proj_dict["repo loc"]],
                proj_dict["repo link"],
            )
            line += gen_table_element(
                [proj_dict["human parent"], proj_dict["parent loc"]],
                proj_dict["parent link"],
            )
            if "ticket no" in proj_dict:
                if proj_dict["ticket no"] is not None:
                    project_ticket_link = [
                        "{0:s}:{1:s}".format(project, proj_dict["ticket no"])
                    ]
                else:
                    project_ticket_link = [None]
            else:
                project_ticket_link = [None]
            line += gen_table_element(project_ticket_link, None)
            wc_link = None
            wc_text = None
            if "working copy changes" in proj_dict:
                if proj_dict["working copy changes"]:
                    wc_text = "YES"
                    wc_link = r"{0:s}/{1:s}/{2:s}/{3:s}?path=log/{4:s}".format(
                        CYLC_REVIEW_URL[self.site],
                        "view",
                        self.suite_owner,
                        self.suitename,
                        proj_dict["version file"],
                    )
            line += gen_table_element([wc_text], wc_link, bold=True)
            lines.append(line)
        return lines

    def gen_resources_table(self):
        """Loops over the RESOURCE_MONITORING_JOBS and returns a
        trac-formatted table of resources used by those jobs in this suite."""
        lines = ["", "  No resource monitoring jobs run"]
        found_nothing = True
        for job in RESOURCE_MONITORING_JOBS[self.site]:
            filename = os.path.join(
                self.suite_path, "log", "job", "1", job, "NN", "job.out"
            )
            if os.path.isfile(filename):
                wallclock, memory = self.get_wallclock_and_memory(filename)
                if wallclock and memory:
                    if found_nothing:
                        lines = [
                            "",
                            " |||| '''Resource Monitoring Task''' || ",
                            " || '''Task''' || '''Wallclock''' || "
                            + "'''Total Memory''' ||",
                        ]
                        found_nothing = False
                    lines.append(
                        " || {0:s} || {1:} || {2:} || ".format(
                            job, wallclock, memory
                        )
                    )
        lines.append("")
        return lines

    @staticmethod
    def get_wallclock_and_memory(filename):
        """Given an output filename read and parse for the wallclock
        and memory."""
        wallclock = "Unavailable"
        memory = "Unavailable"
        find_wallclock = re.compile(
            r"PE\s*0\s*Elapsed Wallclock Time:\s*(\d+(\.\d+|))"
        )
        find_total_mem = re.compile(r"Total Mem\s*(\d+)")
        find_um_atmos_exe = re.compile(r"um-atmos.exe")
        check_for_percentage = re.compile("[0-9]+[%]")
        find_mem_n_units = re.compile(
            r"(?P<num>[0-9]*\.[0-9]*)(?P<unit>[A-Za-z])"
        )
        try:
            for line in _read_file(filename):
                result = find_wallclock.search(line)
                if result:
                    wallclock = int(round(float(result.group(1))))
                result = find_total_mem.search(line)
                if result:
                    memory = int(result.group(1))
                if find_um_atmos_exe.match(line):
                    split_line = line.split()
                    if check_for_percentage.search(split_line[6]):
                        mem = find_mem_n_units.search(split_line[5])
                        memory = float(mem.group("num"))
                        if mem.group("unit") == "G":
                            memory *= 1000000
                        elif mem.group("unit") == "M":
                            memory *= 1000
                    else:
                        memory = int(line.split()[6])
        except Exception as err:
            wallclock = "Failure processing EOJ"
            memory = "Failure processing EOJ"
            stacktr = traceback.format_exc()
            print(
                "[ERROR] Processing wallclock and memory use :\n{0:s}".format(
                    stacktr
                )
            )
            print("Error type : {0:s}".format(type(err)))
            print(err)
        return wallclock, memory

    @staticmethod
    def query_database(suite_db_file):
        """Query the database and return a dictionary of states."""
        database = sqlite3.connect(suite_db_file)
        cursor = database.cursor()
        cursor.execute("select name, status from task_states;")
        data = {}
        for row in cursor:
            data[row[0]] = row[1]
        database.close()
        return data

    @staticmethod
    def ascertain_ticket_number(mirror_url, fcm_exec):
        """Try and work out the ticket number from the Trac log.
        Takes URL on local (mirror) repository and fcm_exec path.
        Uses 'fcm log'
        Relies on commit line starting with '#[0-9]+' - meto working
        practices for commit says "start with ticket number"
        Returns ticket number as string or None."""
        ticket_number = None
        if re.search("/trunk[/@$]", mirror_url) or re.search(
            r"[fs][cv][mn]:\w+(.xm|.x|)_tr[/@$]", mirror_url
        ):
            return ticket_number
        _, stdout, _ = _run_command([fcm_exec, "log", "-l", "1", mirror_url])
        for line in stdout:
            result = re.search(r"^\s*(#\d+)", line)
            if result:
                ticket_number = result.group(1)
        return ticket_number

    @staticmethod
    def get_altered_files_list(mirror_loc):
        """
        Use the get_branch_diff_filenames function from fcm_bdiff to get a list
        of files edited on a branch. Remove any entry that doesn't contain a
        file extension as these are likely directories.
        """

        # Get a list of altered files from the fcm mirror url
        bdiff_files = get_branch_diff_filenames(mirror_loc, path_override="")

        # If '.' is in the files list remove it
        try:
            bdiff_files.remove(".")
        except ValueError:
            pass

        # Remove any item that is not a file - decide based on whether the item
        # has a file extension
        for item in bdiff_files:
            # The last part of the file path
            file_name = item.split("/")[-1]
            # If the final part doesn't have a file extension remove this item
            # as it is a directory
            if "." not in file_name:
                bdiff_files.remove(item)

        return bdiff_files

    @staticmethod
    def generate_groups(grouplist):
        """Convert the list of groups run into a trac-formatted string."""
        output = ""
        for group in grouplist[:-1]:
            output += "{0:s} [[br]]".format(_remove_quotes(group))
        output += "{0:s}".format(_remove_quotes(grouplist[-1]))
        return output

    def print_report(self):
        """'Prints a Trac formatted report of the suite_report object"""
        try:
            trac_log = []
            ticket_nos = ""

            # Check to see if any of the soucres have associated tickets and
            # put links to them in the header if so.
            for project, url_dict in self.job_sources.items():
                try:
                    if url_dict["ticket no"] is not None:
                        ticket_nos += "{0:s}:{1:s} ".format(
                            project, url_dict["ticket no"]
                        )
                except KeyError:
                    pass
            trac_log.append('{{{#!div style="background : #FFFFBF"')
            if ticket_nos != "":
                trac_log.append(
                    " = Ticket {0:s} ".format(ticket_nos)
                    + "Testing Results - rose-stem output = "
                )
            else:
                trac_log.append(" = Testing Results - rose-stem output = ")
            trac_log.append("")

            trac_log.append(
                " || Suite Name: || {0:s} || ".format(self.suitename)
            )

            trac_log.append(
                " || Suite Owner: || {0:s} || ".format(self.suite_owner)
            )

            if self.trustzone:
                trac_log.append(
                    " || Trustzone: || {0:s} || ".format(self.trustzone)
                )

            if self.fcm:
                trac_log.append(
                    " || FCM version: || {0:s} || ".format(self.fcm)
                )

            if self.rose:
                trac_log.append(
                    " || Rose version: || {0:s} || ".format(self.rose)
                )

            if self.cylc:
                trac_log.append(
                    " || Cylc version: || {0:s} || ".format(self.cylc)
                )

            trac_log.append(
                " || Report Generated: || {0:s} || ".format(self.creation_time)
            )

            trac_log.append(
                " || Cylc-Review: || {0:s}/{1:s}/{2:s}/{3:s} || ".format(
                    CYLC_REVIEW_URL[self.site],
                    "taskjobs",
                    self.suite_owner,
                    self.suitename,
                )
            )

            trac_log.append(" || Site: || {0:s} || ".format(self.site))
            trac_log.append(
                " || Groups Run: || {0:s} || ".format(
                    self.generate_groups(self.groups)
                )
            )
            if self.rose_orig_host is not None:
                trac_log.append(
                    " || ''ROSE_ORIG_HOST:'' || {0:s} || ".format(
                        self.rose_orig_host
                    )
                )
            if self.host_xcs:
                trac_log.append(" || HOST_XCS || True || ")
            trac_log.append("")

            if self.uncommitted_changes:
                trac_log.append("")
                trac_log.append("-----")
                word = "changes" if self.uncommitted_changes > 1 else "change"
                trac_log.append(" = WARNING !!! = ")
                trac_log.append(
                    "This rose-stem suite included "
                    + "{0:d} uncommitted".format(self.uncommitted_changes)
                    + " project {0:s} and is therefore ".format(word)
                    + "'''not valid''' for review"
                )
                trac_log.append("-----")
                trac_log.append("")

            if not self.required_comparisons:
                trac_log.append("")
                trac_log.append("-----")
                trac_log.append(" = WARNING !!! = ")
                trac_log.append(
                    "This rose-stem suite did not run the "
                    + "required comparisons (COMPARE_OUTPUT "
                    + "and/or COMPARE_WALLCLOCK are not true) "
                    + "and is therefore '''not valid''' for "
                    + "review"
                )
                trac_log.append("-----")
                trac_log.append("")

            if len(self.multi_branches.keys()) > 0:
                trac_log.append("")
                trac_log.append("-----")
                trac_log.append(" = WARNING !!! = ")
                trac_log.append(
                    "This rose-stem suite included multiple "
                    + "branches in {0:d} projects:".format(
                        len(self.multi_branches.keys())
                    )
                )
                trac_log.append("")
                for project, branch_names in self.multi_branches.items():
                    trac_log.append("'''{0}'''".format(project))
                    for branch_name in "".join(branch_names).split():
                        trac_log.append(" * {0}".format(branch_name))
                trac_log.append("")
                trac_log.append("-----")
                trac_log.append("")

            trac_log.extend(self.generate_project_table())
            trac_log.append("")

            data = self.query_database(
                os.path.join(self.suite_path, SUITE_DB_FILENAME)
            )
            trac_log.extend(
                self.generate_task_table(
                    data,
                    common_groups=self.only_common_groups,
                    verbosity=self.verbosity,
                    sort_by_name=self.sort_by_name,
                    status_counts=self.status_counts,
                )
            )
            trac_log.append("")
            trac_log.append("}}}")
            trac_log.append("")

        except Exception as err:
            print(err)
            try:
                suite_dir = self.suite_path
            except:
                suite_dir = "--cylc_suite_dir--"
            trac_log.extend(
                [
                    "There has been an exception in "
                    + "SuiteReport.print_report()",
                    "See output for more information",
                    "rose-stem suite output will be in the files :\n",
                    "~/cylc-run/{0:s}/log/suite/log".format(suite_dir),
                ]
            )
        finally:
            # Pick up user specified log path if available,
            # otherwise default to cyclc suite dir.
            trac_log_path = "/No/Path/Provided"
            if self.log_path:
                trac_log_path = os.path.join(self.log_path, TRAC_LOG_FILE)
            else:
                trac_log_path = os.path.join(self.suite_path, TRAC_LOG_FILE)

            # Attempt to provide user with some output,
            # even in event of serious exceptions
            try:
                _write_file(trac_log_path, trac_log, newline=True)
            except IOError:
                print(
                    "[ERROR] Writing to {0:s} file : {1:s}".format(
                        TRAC_LOG_FILE, trac_log_path
                    )
                )
                print(
                    "{0:s} to this point ".format(TRAC_LOG_FILE)
                    + "would have read as follows :\n"
                )
                print("----- Start of {0:s}.log -----".format(TRAC_LOG_FILE))
                for line in trac_log:
                    print(line)
                print(
                    "\n----- End of {0:s}.log -----\n\n".format(TRAC_LOG_FILE)
                )
                raise


# ==============================================================================
#    End of   "class.SuiteReport()"
# ==============================================================================


def get_working_copy_path(path):
    """
    Function that tries to find a working copy given a path to that copy.
    Input:
        - path: Path to a working copy. The path provided in this script
                contains the hostname in the format <hostname>:<path>.
                Python seems unable to parse this hence the code below.
    """
    if not os.path.exists(path):
        try:
            path = path.split(":")[1]
            if not os.path.exists(path):
                return None
        except IndexError:
            return None
    return path


def parse_options():
    """Use OptionParser to parse options from command line.
    Would prefer ArgParse but currently unavailable in default site
    Python installation at meto."""
    expected_args = 0
    parser = OptionParser()
    parser.add_option(
        "--increase-verbosity",
        "-v",
        dest="increase_verbosity",
        action="count",
        default=0,
        help="Increases Verbosity level. Default = "
        "{0:2d}".format(DEFAULT_VERBOSITY),
    )
    parser.add_option(
        "--decrease-verbosity",
        "-q",
        dest="decrease_verbosity",
        action="count",
        default=0,
        help="Decreases Verbosity level.                         "
        "  Levels are as follows :                          "
        "0: print everything.                               "
        "1: print except Housekeeping tasks.                "
        "2: print except above and Gatekeeping tasks.       "
        "3: print except above and Successful tasks if only "
        '.  running "common" groups.    [Default]         '
        "4: print except above and all Successful tasks     "
        ".  for all groups/tasks.                           ",
    )
    parser.add_option(
        "--name-sort",
        "-N",
        dest="sort_by_name",
        default=False,
        action="store_true",
        help="Sort task table by task names",
    )
    parser.add_option(
        "--suite-path",
        "-S",
        dest="suite_path",
        action="store",
        help="Path to suite",
    )
    parser.add_option(
        "--log_path",
        "-L",
        dest="log_path",
        action="store",
        default=None,
        help="Output dir for {0:s}".format(TRAC_LOG_FILE),
    )
    # Deprecated options - the options below have been deprecated and should no
    # longer be used or have any effect.
    # "sort by status" is now the default action - "-N" sort by name re-instates
    # the original default sorting mechanism
    # "-H" to "include" the housekeeping tasks has been removed and is now
    # handled by the verbosity flags "-v" and "-q". As is including/excluding
    # the gatekeeper tasks which is new with the version that deprecated these
    # options.
    group = OptionGroup(
        parser,
        "Depricated Options",
        "The following options have been deprecated and no "
        "longer have any effect.  "
        "Please stop using them as they will cause a fatal "
        "unknown option error in the near future",
    )
    group.add_option(
        "--include-housekeeping",
        "-H",
        dest="include_housekeep",
        action="store_true",
        default=False,
        help="Include housekeeping tasks",
    )
    group.add_option(
        "--status-sort",
        "-s",
        dest="sort_by_status",
        default=False,
        action="store_true",
        help="Sort by task status is now default and this option"
        + " is deprecated",
    )
    parser.add_option_group(group)
    # -- End of deprecated options -- block to be removed around Dec 2023

    (opts, args) = parser.parse_args()

    # The calculation below seems counter-intuative. Lower 'verbosity' score
    # actually means a higher verbosity. So SUBTRACT the count for
    # increase_verbosity and ADD the count for decrease_verbosity.
    opts.verbosity = (
        DEFAULT_VERBOSITY - opts.increase_verbosity + opts.decrease_verbosity
    )

    # Find the suite database file
    if opts.suite_path is None:
        if "CYLC_SUITE_RUN_DIR" in os.environ:
            opts.suite_path = os.environ["CYLC_SUITE_RUN_DIR"]
            # cylc provides 3 args to the shutdown handler - so allow for them.
            expected_args = 3
        else:
            sys.exit("Path to suite not provided.")

    if len(args) != expected_args:
        parser.print_help()
        message = "got {0:d} extra arguments, expected {1:d}\n".format(
            len(args), expected_args
        )
        message += "Extra arguments are :\n{0:}".format(args)
        print(message)
        sys.exit(message)

    if opts.sort_by_status or opts.include_housekeep:
        message = []
        message.append("")
        message.append(
            "       ---------------------------------------------"
            "----------------"
        )
        message.append(
            """                            ###  WARNING !!!  ###
                       suite_report.py was called using a deprecated option.
                       Please remove this option from wherever the call was made
                       before this option becomes invalid and this warning
                       becomes a fatal error"""
        )
        message.append(
            "       ---------------------------------------------"
            "----------------"
        )
        message.append("")
        if opts.sort_by_status:
            message.append(
                """
         -s, --status-sort   Sort by task status is now default
                             Use -N to sort by task name (previous default)
            """
            )
        if opts.include_housekeep:
            message.append(
                """
         -H, --include-housekeeping         Include housekeeping tasks
            This is now handled by setting the verbosity level using -[qv]
            """
            )
        message.append(
            """
         Please use "suite_report.py -h" to see full help on options.
        """
        )
        message.append("")
        for line in message:
            print(line)

    return opts


def main():
    """Main program.
    Sets up a SuiteReport object and calls it's print_report method."""
    opts = parse_options()

    suite_report_obj = SuiteReport(
        suite_path=opts.suite_path,
        log_path=opts.log_path,
        verbosity=(opts.verbosity),
        sort_by_name=opts.sort_by_name,
    )

    suite_report_obj.print_report()


if __name__ == "__main__":
    main()
