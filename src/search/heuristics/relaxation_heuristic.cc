#include "relaxation_heuristic.h"

#include "../task_utils/task_properties.h"

#include <algorithm>
#include <cassert>
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
      reached_by(nullptr),
      num_precondition_occurrences(-1),
      is_goal(false),
      marked(false) {
}

void PropositionNode::update_precondition(PropQueue &queue, GraphNode *predecessor) {
    auto *op_predecessor = static_cast<OperatorNode *>(predecessor);
    assert(prop_id != relaxation_heuristic::NO_PROP);
    assert(op_predecessor->cost >= 0);
    int newcost = op_predecessor->cost + op_predecessor->base_cost;
    if (cost == -1 || cost > newcost ) {
        cost = newcost;
        reached_by = op_predecessor;
        queue.push(cost, this);
    }
    assert(cost != -1 && cost <= newcost);
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
    for (OperatorProxy op : task_proxy.get_operators())
        build_unary_operators(op);
    for (OperatorProxy axiom : task_proxy.get_axioms())
        build_unary_operators(axiom);
}

RelaxationHeuristic::~RelaxationHeuristic() {
    for (auto ptr: operator_nodes) {
        delete ptr;
    }
    for (auto ptr: propositions) {
        delete ptr;
    }
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


PropID RelaxationHeuristic::get_prop_id(int var, int value) const {
    return proposition_offsets[var] + value;
}

PropID RelaxationHeuristic::get_prop_id(const FactProxy &fact) const {
    return get_prop_id(fact.get_variable().get_id(), fact.get_value());
}

void RelaxationHeuristic::build_unary_operators(const OperatorProxy &op) {
    int op_no = op.is_axiom() ? -1 : op.get_id();
    PreconditionsProxy preconditions = op.get_preconditions();
    vector<PropositionNode*> precondition_props;
    precondition_props.reserve(preconditions.size());
    for (FactProxy precondition : preconditions) {
        precondition_props.push_back(propositions[get_prop_id(precondition)]);
    }

    OperatorNode* operator_node = new OperatorNode(op.get_cost(),
                                                   precondition_props.size(),
                                                   op_no);
    operator_nodes.push_back(operator_node);
    for (auto precondition: precondition_props) {
        precondition->precondition_of.push_back(operator_node);
        operator_node->preconditions.push_back(precondition);
    }
    for (EffectProxy effect : op.get_effects()) {
        PropositionNode* effect_prop = propositions[get_prop_id(effect.get_fact())];
        EffectConditionsProxy eff_conds = effect.get_conditions();

        if (eff_conds.empty()) {
            operator_node->precondition_of.push_back(effect_prop);
        } else {
            vector<PropositionNode*> effect_conditions;
            effect_conditions.reserve(eff_conds.size());
            for (FactProxy eff_cond : eff_conds) {
                effect_conditions.push_back(propositions[get_prop_id(eff_cond)]);
            }
            OperatorNode* conditional_effect = new OperatorNode(op.get_cost(),
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

}
