#include "additive_heuristic.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>
#include <vector>

using namespace std;

namespace additive_heuristic {
const int AdditiveHeuristic::MAX_COST_VALUE;

// construction and destruction
AdditiveHeuristic::AdditiveHeuristic(const Options &opts)
    : RelaxationHeuristic(opts),
      did_write_overflow_warning(false) {
    if (log.is_at_least_normal()) {
        log << "Initializing additive heuristic..." << endl;
    }
}

void AdditiveHeuristic::write_overflow_warning() {
    if (!did_write_overflow_warning) {
        // TODO: Should have a planner-wide warning mechanism to handle
        // things like this.
        if (log.is_warning()) {
            log << "WARNING: overflow on h^add! Costs clamped to "
                << MAX_COST_VALUE << endl;
        }
        cerr << "WARNING: overflow on h^add! Costs clamped to "
             << MAX_COST_VALUE << endl;
        did_write_overflow_warning = true;
    }
}

// heuristic computation
void AdditiveHeuristic::setup_exploration_queue(const State &state) {
    queue.clear();

    for (auto prop : propositions) {
        prop->cost = -1;
        prop->marked = false;
    }

    // Deal with operators and axioms without preconditions_props.
    for (auto &op : operators) {
        op->unsatisfied_preconditions = op->num_preconditions;
        op->cost = 0;
        if (op->unsatisfied_preconditions == 0) {
            update_precondition(queue, op);
        }
    }

    for (FactProxy fact : state) {
        Proposition* prop = propositions[get_prop_id(fact)];
        prop->cost = 0;
        prop->reached_by = nullptr;
        queue.push(0, prop);
    }
}

void AdditiveHeuristic::relaxed_exploration() {
    int unsolved_goals = goal_propositions.size();
    while (!queue.empty()) {
        pair<int, Proposition*> top_pair = queue.pop();
        int distance = top_pair.first;
        Proposition *prop = top_pair.second;
        int prop_cost = prop->cost;
        assert(prop_cost >= 0);
        assert(prop_cost <= distance);
        if (prop_cost < distance)
            continue;
        if (prop->is_goal && --unsolved_goals == 0)
            return;
        update_precondition(queue, prop);
    }
}

void AdditiveHeuristic::mark_preferred_operators(
    const State &state, Proposition* goal) {
    if (!goal->marked) { // Only consider each subgoal once.
        goal->marked = true;
        Operator* op_id = goal->reached_by;
        if (op_id != nullptr) { // We have not yet chained back to a start node.
            bool is_preferred = true;
            for (auto *precond : preconds_pool.get_slice(op_id->precondition_index, op_id->precondition_size)) {
                mark_preferred_operators(state, precond);
                if (precond->reached_by != nullptr) {
                    is_preferred = false;
                }
            }
            int operator_no = op_id->operator_no;
            if (is_preferred && operator_no != -1) {
                // This is not an axiom.
                OperatorProxy op = task_proxy.get_operators()[operator_no];
                assert(task_properties::is_applicable(op, state));
                set_preferred(op);
            }
        }
    }
}

int AdditiveHeuristic::compute_add_and_ff(const State &state) {
    setup_exploration_queue(state);
    relaxed_exploration();

    int total_cost = 0;
    for (PropID goal_id : goal_propositions) {
        const Proposition *goal = propositions[goal_id];
        int goal_cost = goal->cost;
        if (goal_cost == -1)
            return DEAD_END;
// TODO:        increase_cost(total_cost, goal_cost);
    total_cost += goal_cost;
    }
    return total_cost;
}

int AdditiveHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    int h = compute_add_and_ff(state);
    if (h != DEAD_END) {
        for (PropID goal_id : goal_propositions)
            mark_preferred_operators(state, get_proposition(goal_id));
    }
    return h;
}

void AdditiveHeuristic::compute_heuristic_for_cegar(const State &state) {
    compute_heuristic(state);
}

static shared_ptr<Heuristic> _parse(OptionParser &parser) {
    parser.document_synopsis("Additive heuristic", "");
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
        return make_shared<AdditiveHeuristic>(opts);
}

static Plugin<Evaluator> _plugin("add", _parse);
}
