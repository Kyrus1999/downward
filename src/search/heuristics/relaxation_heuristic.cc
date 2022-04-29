#include "relaxation_heuristic.h"

#include "../task_utils/task_properties.h"
#include "../utils/collections.h"
#include "../utils/logging.h"
#include "../utils/timer.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <unordered_map>
#include <vector>

using namespace std;
namespace relaxation_heuristic {

GraphNode::GraphNode() : cost(-1) {}
GraphNode::GraphNode(std::vector<GraphNode*> &&precondition_of)
    : cost(-1),
      precondition_of(move(precondition_of)) {}

PropositionNode::PropositionNode(PropID prop_id)
    : GraphNode(),
      prop_id(prop_id),
      reached_by(NO_OP),
      is_goal(false),
      marked(false),
      num_precondition_occurrences(-1) {
}

void PropositionNode::update_precondition(PropQueue &queue, GraphNode *predecessor) {
    auto *op_predecessor = static_cast<OperatorNode *>(predecessor);
    assert(prop_id != relaxation_heuristic::NO_PROP);
    assert(op_predecessor->cost >= 0);
    if (cost == -1 || cost > op_predecessor->cost) {
        cost = op_predecessor->cost;
        reached_by = op_predecessor->operator_no;
        queue.push(cost, this);
    }
    assert(cost != -1 && cost <= op_predecessor->cost);
}

void PropositionNode::update_precondition(PropQueue &queue) {
    for (auto node : precondition_of) {
        node->update_precondition(queue, this);
    }
}


OperatorNode::OperatorNode(int base_cost, int num_preconditions, int operator_no)
    : GraphNode(),
      base_cost(base_cost),
      num_preconditions(num_preconditions),
      operator_no(operator_no) {
}


void OperatorNode::update_precondition(PropQueue &queue, GraphNode *predecessor) {
    this->unsatisfied_preconditions--;
    assert(this->unsatisfied_preconditions >= 0);
    this->cost += predecessor->cost; // TODO: check for overflow
    if (this->unsatisfied_preconditions <= 0) {
        for (auto node : this->precondition_of) {
                node->update_precondition(queue, this);
            }
        }
    }


// construction and destruction
RelaxationHeuristic::RelaxationHeuristic(const options::Options &opts)
    : Heuristic(opts) {
    // Build propositions.
    int num_propositions = task_properties::get_num_facts(task_proxy);
//    propositions.reserve(num_propositions);
    for (PropID prop_id = 0; prop_id < num_propositions; ++prop_id) {
        PropositionNode* prop_node = new PropositionNode(prop_id);
        propositions.push_back(prop_node);
    }

    // Build proposition offsets.
    VariablesProxy variables = task_proxy.get_variables();
    proposition_offsets.reserve(variables.size());
    PropID offset = 0;
    for (VariableProxy var : variables) {
        proposition_offsets.push_back(offset);
        offset += var.get_domain_size();
    }
    assert(offset == static_cast<int>(propositions.size()));

    // Build goal propositions.
    GoalsProxy goals = task_proxy.get_goals();
    goal_propositions.reserve(goals.size());
    for (FactProxy goal : goals) {
        PropID prop_id = get_prop_id(goal);
        propositions[prop_id]->is_goal = true;
        goal_propositions.push_back(prop_id);
    }

    // Build unary operators for operators and axioms.
    //unary_operators.reserve(task_properties::get_num_total_effects(task_proxy));

    //unary_operators.reserve(task_properties::get_num_operators(task_proxy) + get_num_cond_effects());
    for (OperatorProxy op : task_proxy.get_operators())
        build_unary_operators(op);
    for (OperatorProxy axiom : task_proxy.get_axioms())
        build_unary_operators(axiom);

//    for (auto opn: operator_nodes) {
//        assert(opn->cost == -1);
//    }
//    for (PropositionNode* pop: propositions) {
//        for (auto node: pop->precondition_of){
//            //cout << node->cost << endl;
//            assert (node->cost == -1); // Why is this? g is an operator node, therefor should be 0?
//        }
//    }
    // Simplify unary operators.
    utils::Timer simplify_timer;
    simplify();
    if (log.is_at_least_normal()) {
        log << "time to simplify: " << simplify_timer << endl;
    }

    // Not required anymore
//    // Cross-reference unary operators.
//    vector<vector<OpID>> precondition_of_vectors(propositions.size());
//
//    int num_unary_ops = operator_nodes.size();
//    for (OpID op_id = 0; op_id < num_unary_ops; ++op_id) {
//        for (PropID precond : get_preconditions(op_id))
//            precondition_of_vectors[precond].push_back(op_id);
//    }
//
//    for (PropID prop_id = 0; prop_id < num_propositions; ++prop_id) {
//        const auto &precondition_of_vec = precondition_of_vectors[prop_id];
//        propositions[prop_id].precondition_of =
//            precondition_of_pool.append(precondition_of_vec);
//        propositions[prop_id].num_precondition_occurences = precondition_of_vec.size();
//    }
}

bool RelaxationHeuristic::dead_ends_are_reliable() const {
    return !task_properties::has_axioms(task_proxy);
}

int RelaxationHeuristic::get_num_cond_effects() {
    int counter = 0;
    for (OperatorProxy op : task_proxy.get_operators() ){
        for (EffectProxy effect: op.get_effects()) {
            EffectConditionsProxy eff_conds = effect.get_conditions();

            if (eff_conds.empty()) counter++;
        }
    }
    return counter;
}

//PropositionNode &RelaxationHeuristic::get_prop_node(int var, int value) const {
//    return &propositions[proposition_offsets[var] + value];
//}

//PropositionNode* RelaxationHeuristic::get_prop_node(const FactProxy &fact) const {
//    return get_prop_node(fact.get_variable().get_id(), fact.get_value());
//}

PropID RelaxationHeuristic::get_prop_id(int var, int value) const {
    return proposition_offsets[var] + value;
}

PropID RelaxationHeuristic::get_prop_id(const FactProxy &fact) const {
    return get_prop_id(fact.get_variable().get_id(), fact.get_value());
}
//
//const Proposition *RelaxationHeuristic::get_proposition(
//    int var, int value) const {
//    return &propositions[get_prop_id(var, value)];
//}
//
//Proposition *RelaxationHeuristic::get_proposition(int var, int value) {
//    return &propositions[get_prop_id(var, value)];
//}
//
//Proposition *RelaxationHeuristic::get_proposition(const FactProxy &fact) {
//    return get_proposition(fact.get_variable().get_id(), fact.get_value());
//}

void RelaxationHeuristic::build_unary_operators(const OperatorProxy &op) {
    int op_no = op.is_axiom() ? -1 : op.get_id();
    PreconditionsProxy preconditions = op.get_preconditions();
    vector<PropositionNode*> precondition_props;
    precondition_props.reserve(preconditions.size());
    for (FactProxy precondition : preconditions) {
        precondition_props.push_back(propositions[get_prop_id(precondition)]);
    }

    utils::sort_unique(precondition_props); // Why?

//    array_pool::ArrayPoolIndex precond_index =
//            preconditions_pool.append(preconditions_copy);
//    operator_nodes.emplace_back(preconditions_copy.size(), precond_index);
    OperatorNode* operator_node = new OperatorNode(op.get_cost(),
                                                   precondition_props.size(),
                                                   op_no);
    operator_nodes.push_back(operator_node);
    for (auto precondition: precondition_props) {
        precondition->precondition_of.push_back(operator_node);
        operator_node->preconditions.push_back(precondition);
    }
//    operator_node.precondition_of.reserve(op.get_effects().size());
    for (EffectProxy effect : op.get_effects()) {
        PropositionNode* effect_prop = propositions[get_prop_id(effect.get_fact())];
        EffectConditionsProxy eff_conds = effect.get_conditions();

        if (eff_conds.empty()) {
            operator_node->precondition_of.push_back(effect_prop);
        } else {
//            array_pool::ArrayPoolIndex effect_preconds_index =
//                    preconditions_pool.append(effect_preconds);
            vector<PropositionNode*> effect_conditions;
            effect_conditions.reserve(eff_conds.size());
            for (FactProxy eff_cond : eff_conds) {
                effect_conditions.push_back(propositions[get_prop_id(eff_cond)]);
            }
            //utils::sort_unique(effect_conditions);
            OperatorNode* conditional_effect = new OperatorNode(0,
                                                                effect_conditions.size() + 1,
                                                                op_no);
            operator_nodes.push_back(conditional_effect);
            for (auto precondition: effect_conditions) {
                precondition->precondition_of.push_back(conditional_effect);
                conditional_effect->preconditions.push_back(precondition);
            }
            for (auto precondition: precondition_props) {
                conditional_effect->preconditions.push_back(precondition);
            }
            utils::sort_unique(conditional_effect->preconditions);
            operator_node->precondition_of.push_back(conditional_effect);
            conditional_effect->precondition_of.push_back(effect_prop);
        }

    }
}


int RelaxationHeuristic::get_proposition_cost(int var, int value) const {
    PropID prop_id = get_prop_id(var, value);
    const PropositionNode *node = propositions[prop_id];
    return node->cost;
}

void RelaxationHeuristic::simplify() {
    /*
      Remove dominated unary operators, including duplicates.

      Unary operators with more than MAX_PRECONDITIONS_TO_TEST
      preconditions_props are (mostly; see code comments below for details)
      ignored because we cannot handle them efficiently. This is
      obviously an inelegant solution.

      Apart from this restriction, operator o1 dominates operator o2 if:
      1. eff(o1) = eff(o2), and
      2. pre(o1) is a (not necessarily strict) subset of pre(o2), and
      3. cost(o1) <= cost(o2), and either
      4a. At least one of 2. and 3. is strict, or
      4b. id(o1) < id(o2).
      (Here, "id" is the position in the unary_operators vector.)

      This defines a strict partial order.
    */
#ifndef NDEBUG
    int num_ops = operator_nodes.size();
    for (OpID op_id = 0; op_id < num_ops; ++op_id)
        assert(utils::is_sorted_unique(get_preconditions(op_id)));
#endif

    const int MAX_PRECONDITIONS_TO_TEST = 5;

    if (log.is_at_least_normal()) {
        log << "Simplifying " << operator_nodes.size() << " unary operators..." << flush;
    }

    /*
      First, we create a map that maps the preconditions_props and effect
      ("key") of each operator to its cost and index ("value").
      If multiple operators have the same key, the one with lowest
      cost wins. If this still results in a tie, the one with lowest
      index wins. These rules can be tested with a lexicographical
      comparison of the value.

      Note that for operators sharing the same preconditions_props and
      effect, our dominance relationship above is actually a strict
      *total* order (order by cost, then by id).

      For each key present in the data, the map stores the dominating
      element in this total order.
    */
    using Key = pair<vector<PropositionNode*>, vector<GraphNode*>>;
    using Value = pair<int, OpID>;
    using Map = utils::HashMap<Key, Value>;
    Map unary_operator_index;
    unary_operator_index.reserve(operator_nodes.size());
    for (size_t op_no = 0; op_no < operator_nodes.size(); ++op_no) {
        const OperatorNode *op = operator_nodes[op_no];
        /*
          Note: we consider operators with more than
          MAX_PRECONDITIONS_TO_TEST preconditions_props here because we can
          still filter out "exact matches" for these, i.e., the first
          test in `is_dominated`.
        */
        Key key(get_preconditions(op_no), op->precondition_of);
        Value value(op->base_cost, op_no);
        auto inserted = unary_operator_index.insert(
            make_pair(move(key), value));
        if (!inserted.second) {
            // We already had an element with this key; check its cost.
            Map::iterator iter = inserted.first;
            Value old_value = iter->second;
            if (value < old_value)
                iter->second = value;
        }
    }
    /*
      `dominating_key` is conceptually a local variable of `is_dominated`.
      We declare it outside to reduce vector allocation overhead.
    */
    Key dominating_key;
    /*
      is_dominated: test if a given operator is dominated by an
      operator in the map.
    */
    auto is_dominated = [&](const OperatorNode &op) {
            /*
              Check all possible subsets X of pre(op) to see if there is a
              dominating operator with preconditions_props X represented in the
              map.
            */

            OpID op_id = get_op_id(op);
            int cost = op.base_cost;

            const vector<PropositionNode*> precondition = get_preconditions(op_id);

            /*
              We handle the case X = pre(op) specially for efficiency and
              to ensure that an operator is not considered to be dominated
              by itself.

              From the discussion above that operators with the same
              precondition and effect are actually totally ordered, it is
              enough to test here whether looking up the key of op in the
              map results in an entry including op itself.
            */
            if (unary_operator_index[make_pair(precondition, op.precondition_of)].second != op_id)
                return true;

            /*
              We now handle all cases where X is a strict subset of pre(op).
              Our map lookup ensures conditions 1. and 2., and because X is
              a strict subset, we also have 4a (which means we don't need 4b).
              So it only remains to check 3 for all hits.
            */
            if (op.num_preconditions > MAX_PRECONDITIONS_TO_TEST) {
                /*
                  The runtime of the following code grows exponentially
                  with the number of preconditions_props.
                */
                return false;
            }

            vector<PropositionNode*> &dominating_precondition = dominating_key.first;
            dominating_key.second = op.precondition_of;

            // We subtract "- 1" to generate all *strict* subsets of precondition.
            int powerset_size = (1 << precondition.size()) - 1;
            for (int mask = 0; mask < powerset_size; ++mask) {
                dominating_precondition.clear();
                for (size_t i = 0; i < precondition.size(); ++i)
                    if (mask & (1 << i))
                        dominating_precondition.push_back(precondition[i]);
                Map::iterator found = unary_operator_index.find(dominating_key);
                if (found != unary_operator_index.end()) {
                    Value dominator_value = found->second;
                    int dominator_cost = dominator_value.first;
                    if (dominator_cost <= cost)
                        return true;
                }
            }
            return false;
        };

    int removed = 0;
    for (OperatorNode *np : operator_nodes) {
        if (np->operator_no == -1) { //axioms can not be optimized
            continue;
        }
        if (is_dominated(*np)) {
            for (PropositionNode* prop : np->preconditions) {
                prop->precondition_of.erase(std::remove(prop->precondition_of.begin(), prop->precondition_of.end(), np), prop->precondition_of.end());
            }
            removed++;
        }
    }
    cout << "Removed: " << removed << endl;
//    unary_operators.erase(
//            remove_if(
//                    unary_operators.begin(),
//                    unary_operators.end(),
//                    is_dominated),
//            unary_operators.end());

    if (log.is_at_least_normal()) {
        log << " done! [" << operator_nodes.size() << " unary operators]" << endl;
    }
}

RelaxationHeuristic::~RelaxationHeuristic() {
    for (auto ptr: operator_nodes) {
        delete ptr;
    }
    for (auto ptr: propositions) {
        delete ptr;
    }
}
}
