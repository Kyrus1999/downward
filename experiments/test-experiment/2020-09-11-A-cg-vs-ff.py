#! /usr/bin/env python

import os

import project

import downward.reports as rep


REPO = project.get_repo_base()
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
if project.REMOTE:
    SUITE = project.SUITE_SATISFICING
    SUITE = ["caldera-sat18-adl",
             "caldera-split-sat18-adl",
             "citycar-sat14-adl",
             "organic-synthesis-split-sat18-strips",
             "parcprinter-sat11-strips",]

    ENV = project.BaselSlurmEnvironment(email="cyrill.imahorn@stud.unibas.ch")
else:
    SUITE = ["depot:p01.pddl", "grid:prob01.pddl", "gripper:prob01.pddl",]
    ENV = project.LocalEnvironment(processes=6)

CONFIGS = [
    #("add", ["--search", f"eager_greedy([add()])"]),
    ("ff", ["--search", f"eager_greedy([ff()])"]),
]

BUILD_OPTIONS = ["--debug"]
DRIVER_OPTIONS = ["--overall-time-limit", "30m", "--debug"]
REVS = [
    ("main", "base"),
    #("version1", "version1"),
    ("version1-simplify", "version1-simplify"),
    #("patrick3", "iter"),
    #("iter_test", "iter_ap"),
    #("iter_test_2", "iter_ap_reduced")
]
ATTRIBUTES = [
    "error",
    "run_dir",
    "search_start_time",
    "search_start_memory",
    "total_time",
    "search_time",
    "initial_h_value",
    "h_values",
    "coverage",
    "expansions",
    "memory",
    project.EVALUATIONS_PER_TIME,
]

exp = project.CommonExperiment(environment=ENV)
for config_nick, config in CONFIGS:
    for rev, rev_nick in REVS:
        algo_name = f"{rev_nick}:{config_nick}" if rev_nick else config_nick
        exp.add_algorithm(
            algo_name,
            REPO,
            rev,
            config,
            build_options=BUILD_OPTIONS,
            driver_options=DRIVER_OPTIONS,
        )
        algo = exp._algorithms[algo_name]
        if not project.REMOTE:
            algo.driver_options.remove("--validate")
            index = algo.driver_options.index("--overall-memory-limit")
            del algo.driver_options[index + 1]
            del algo.driver_options[index]
exp.add_suite(BENCHMARKS_DIR, SUITE)

project.add_absolute_report(
    exp, attributes=ATTRIBUTES, filter=[project.add_evaluations_per_time]
)

# attributes = ["expansions"]
# pairs = [
#     ("20.06:01-cg", "20.06:02-ff"),
# ]
# suffix = "-rel" if project.RELATIVE else ""
# for algo1, algo2 in pairs:
#     for attr in attributes:
#         exp.add_report(
#             project.ScatterPlotReport(
#                 relative=project.RELATIVE,
#                 get_category=None if project.TEX else lambda run1, run2: run1["domain"],
#                 attributes=[attr],
#                 filter_algorithm=[algo1, algo2],
#                 filter=[project.add_evaluations_per_time],
#                 format="tex" if project.TEX else "png",
#             ),
#             name=f"{exp.name}-{algo1}-vs-{algo2}-{attr}{suffix}",
#         )

class TranslatorDiffReport(rep.PlanningReport):
    def get_cell(self, run):
        return ";".join(run.get(attr) for attr in self.attributes)

    def get_text(self):
        lines = []
        for runs in self.problem_runs.values():
            #lhashes = [r.get("translator_output_sas_hash") for r in runs]
            #hashes = set(lhashes)
            h_inits = [r.get("initial_h_value") for r in runs]
            found = False
            reason = ""
            for i in range(0, len(h_inits), 1):
                if h_inits[i] != h_inits[0]:
                    found = True
                    reason += str(i) + ", "

            #for i in range(1, len(h_inits), 2):
            #    if h_inits[i] != h_inits[1]:
            #        found = True
            #        reason += str(i) + ", "


            # if None in hashes:
            #     reason = f"{len([h for h in lhashes if h is None])} failed + "
            # if len(hashes) > 1:
            #     reason += f"{len([h for h in lhashes if h is not None])} differ"
            # if len(reason):
            #     lines.append(reason + ";" + ";".join([self.get_cell(r) for r in runs]))
            if found:
                lines.append(reason + ";" + ";".join([self.get_cell(r) for r in runs]))
        return "\n".join(lines)



exp.add_report(TranslatorDiffReport(attributes=["domain", "problem", "algorithm", "run_dir"]), outfile="different_output_sas.csv")
exp.run_steps()
