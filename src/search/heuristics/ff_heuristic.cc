#include "ff_heuristic.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>

using namespace std;

namespace ff_heuristic {
// construction and destruction
FFHeuristic::FFHeuristic(const Options &opts)
    : AdditiveHeuristic(opts),
      relaxed_plan(task_proxy.get_operators().size(), false) {
    if (log.is_at_least_normal()) {
        log << "Initializing FF heuristic..." << endl;
    }
}

void FFHeuristic::mark_preferred_operators_and_relaxed_plan(
    const State &state, relaxation_heuristic::PropositionNode* goal) {
    if (!goal->marked) { // Only consider each subgoal once.
        goal->marked = true;
        OperatorNode* op_id = goal->reached_by;
        if (op_id != nullptr) { // We have not yet chained back to a start node.
            bool is_preferred = true;
            for (auto *precond : op_id->preconditions) {
                mark_preferred_operators_and_relaxed_plan(
                    state, precond);
                if (precond->reached_by != nullptr) {
                    is_preferred = false;
                }
            }
            int operator_no = op_id->operator_no;
            if (operator_no != -1) {
                // This is not an axiom.
                relaxed_plan[operator_no] = true;
                if (is_preferred) {
                    OperatorProxy op = task_proxy.get_operators()[operator_no];
                    assert(task_properties::is_applicable(op, state));
                    set_preferred(op);
                }
            }
        }
    }
}

int FFHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    int h_add = compute_add_and_ff(state);
    if (h_add == DEAD_END)
        return h_add;

    // Collecting the relaxed plan also sets the preferred operators.
    for (PropID goal_id : goal_propositions)
        mark_preferred_operators_and_relaxed_plan(state, get_proposition(goal_id));

    int h_ff = 0;
    for (size_t op_no = 0; op_no < relaxed_plan.size(); ++op_no) {
        if (relaxed_plan[op_no]) {
            relaxed_plan[op_no] = false; // Clean up for next computation.
            h_ff += task_proxy.get_operators()[op_no].get_cost();
        }
    }
    cout << "h_add " << h_add << endl;
    cout << "h_ff  " << h_ff  << endl;
    utils::exit_with(utils::ExitCode::SEARCH_UNSUPPORTED);
    return h_ff;
}


static shared_ptr<Heuristic> _parse(OptionParser &parser) {
    parser.document_synopsis("FF heuristic", "");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "supported");
    parser.document_language_support(
        "axioms",
        "supported (in the sense that the planner won't complain -- "
        "handling of axioms might be very stupid "
        "and even render the heuristic unsafe)");
    parser.document_property("admissible", "no");
    parser.document_property("consistent", "no");
    parser.document_property("safe", "yes for tasks without axioms");
    parser.document_property("preferred operators", "yes");

    Heuristic::add_options_to_parser(parser);
    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<FFHeuristic>(opts);
}

static Plugin<Evaluator> _plugin("ff", _parse);
}
