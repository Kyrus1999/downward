#include "relaxation_heuristic.h"

#include "../task_utils/task_properties.h"

#include <algorithm>
#include <cassert>
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
      reached_by(nullptr),
      num_precondition_occurrences(-1),
      is_goal(false),
      marked(false) {
}

void PropositionNode::update_precondition(PropQueue &queue, GraphNode *predecessor) {
    auto *op_predecessor = static_cast<OperatorNode *>(predecessor);
    assert(op_predecessor->cost >= 0);
    int newcost = op_predecessor->cost + op_predecessor->base_cost;
    string s1 = to_string(op_predecessor->operator_no);
    string s2 = to_string(this->prop_id);
    printf("%5s: %7s; %5s: %7s; %5s: %7s\n", "OpID", s1.c_str() , "PropID", s2.c_str(), "Cost", to_string(newcost).c_str() );
    if (cost == -1 || cost > newcost) {
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
      operator_no(operator_no),
      parent_node(nullptr) {
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
    for (OperatorProxy op : task_proxy.get_operators())
        build_unary_operators(op);
    for (OperatorProxy axiom : task_proxy.get_axioms())
        build_unary_operators(axiom);

    // Simplify unary operators.
    utils::Timer simplify_timer;
    simplify();
    if (log.is_at_least_normal()) {
        log << "time to simplify: " << simplify_timer << endl;
    }
}

bool RelaxationHeuristic::dead_ends_are_reliable() const {
    return !task_properties::has_axioms(task_proxy);
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

    //utils::sort_unique(precondition_props);
    sort_vector_by_propid(precondition_props);

    OperatorNode* operator_node = new OperatorNode(op.get_cost(),
                                                   precondition_props.size(),
                                                   op_no);
    operator_nodes.push_back(operator_node);
    for (auto precondition: precondition_props) {
        precondition->precondition_of.push_back(operator_node);
        operator_node->preconditions.push_back(precondition);
    }
    vector<OperatorNode*> temp;
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
            conditional_effect->parent_node = operator_node;
            for (auto precondition: effect_conditions) {
                precondition->precondition_of.push_back(conditional_effect);
                conditional_effect->preconditions.push_back(precondition);
            }
            for (auto precondition: precondition_props) {
                conditional_effect->preconditions.push_back(precondition);
            }
            //utils::sort_unique(conditional_effect->preconditions);
            sort_vector_by_propid(conditional_effect->preconditions);
            //operator_node->precondition_of.push_back(conditional_effect);
            temp.push_back(conditional_effect);
            conditional_effect->precondition_of.push_back(effect_prop);
        }

    }
    for (OperatorNode* t : temp) {
        operator_node->precondition_of.push_back(t);
    }
}

void RelaxationHeuristic::sort_vector_by_propid(std::vector<PropositionNode*> &vector) {
    std::vector<PropID> id_vector;
    id_vector.reserve(vector.size());
    for (PropositionNode *p : vector) {
        id_vector.push_back(p->prop_id);
    }
    std::vector<PropID> sorted_id_vector(id_vector);
    utils::sort_unique(sorted_id_vector);
    std::vector<PropositionNode*> temp;
    for (PropID prop_id : sorted_id_vector) {
        auto iter = std::find(id_vector.begin(), id_vector.end(), prop_id);
        assert(iter != id_vector.end());
        temp.push_back(vector[iter - id_vector.begin()]);
    }
    vector.swap(temp);
}

bool RelaxationHeuristic::is_sorted_by_propid(std::vector<PropositionNode *> &vector) {
    std::vector<PropID> id_vector;
    id_vector.reserve(vector.size());
    for (PropositionNode *p : vector) {
        id_vector.push_back(p->prop_id);
    }
    return utils::is_sorted_unique(id_vector);
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
    for (OperatorNode *op: operator_nodes)
        assert(is_sorted_by_propid(op->preconditions));
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
    using Key = pair<vector<PropositionNode *>, GraphNode *>;
    using Value = pair<int, GraphNode *>;
    using Map = utils::HashMap<Key, Value>;
    Map unary_operator_index;
    unary_operator_index.reserve(operator_nodes.size());
    for (auto *op: operator_nodes) {
        /*
          Note: we consider operators with more than
          MAX_PRECONDITIONS_TO_TEST preconditions_props here because we can
          still filter out "exact matches" for these, i.e., the first
          test in `is_dominated`.
        */
        for (GraphNode *effect: op->precondition_of) {
            Key key(op->preconditions, effect);
            Value value(op->base_cost, op);
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
    }
    /*
      `dominating_key` is conceptually a local variable of `is_dominated`.
      We declare it outside to reduce vector allocation overhead.
    */
    Key dominating_key;


    int counter_deleted_effects = 0;
    int counter_deleted_nodes = 0;
    /*
      is_dominated: test if a given operator is dominated by an
      operator in the map.
    */
    for (int i = operator_nodes.size()-1; i >= 0; i--) {
        OperatorNode *op = operator_nodes.at(i);
        /*
          Check all possible subsets X of pre(op) to see if there is a
          dominating operator with preconditions_props X represented in the
          map.
        */

        int cost = op->base_cost;
        //const vector<PropositionNode *> precondition = op->preconditions;

        /*
          We handle the case X = pre(op) specially for efficiency and
          to ensure that an operator is not considered to be dominated
          by itself.

          From the discussion above that operators with the same
          precondition and effect are actually totally ordered, it is
          enough to test here whether looking up the key of op in the
          map results in an entry including op itself.
        */
        for (GraphNode *effect : op->precondition_of) {
            if (op != unary_operator_index[make_pair(op->preconditions, effect)].second) {
                for (unsigned long index = 0; index < op->precondition_of.size(); index++) {
                    if (op->precondition_of[index] == effect) {
                        op->precondition_of.erase(op->precondition_of.begin() + index);
                        break;
                    }
                }
                //std::remove(op->precondition_of.begin(), op->precondition_of.end(),
                //            effect);

                counter_deleted_effects++;
            }
        }
        /*
          We now handle all cases where X is a strict subset of pre(op).
          Our map lookup ensures conditions 1. and 2., and because X is
          a strict subset, we also have 4a (which means we don't need 4b).
          So it only remains to check 3 for all hits.
        */
        if (op->preconditions.size() > MAX_PRECONDITIONS_TO_TEST && !op->precondition_of.empty()) {
            /*
              The runtime of the following code grows exponentially
              with the number of preconditions_props.
            */
            continue;
        }
        vector<PropositionNode *> &dominating_precondition = dominating_key.first;
        for (GraphNode *effect : op->precondition_of) {
            dominating_key.second = effect;

            // We subtract "- 1" to generate all *strict* subsets of precondition.
            int powerset_size = (1 << op->preconditions.size()) - 1;
            for (int mask = 0; mask < powerset_size; ++mask) {
                dominating_precondition.clear();
                for (size_t i = 0; i < op->preconditions.size(); ++i)
                    if (mask & (1 << i))
                        dominating_precondition.push_back(op->preconditions[i]);
                Map::iterator found = unary_operator_index.find(dominating_key);
                if (found != unary_operator_index.end()) {
                    Value dominator_value = found->second;
                    int dominator_cost = dominator_value.first;
                    if (dominator_cost <= cost) {
                        //why is this not executed? (specially the remove part)
                        #ifndef NDEBUG
                            unsigned int size_before = op->precondition_of.size() ;
                        #endif
                        for (unsigned long index = 0; index < op->precondition_of.size(); index++) {
                            if (op->precondition_of[index] == effect) {
                                op->precondition_of.erase(op->precondition_of.begin() + index);
                                break;
                            }
                        }
                        //std::remove(op->precondition_of.begin(), op->precondition_of.end(), effect);
                        counter_deleted_effects++;
                        assert(size_before != op->precondition_of.size() );
                        break;
                    }
                }
            }
        }
        // a bit inefficient, since not all precondition propositions for a conditional node contains a link to the conditional node, but to its parent node.
        //could be made more efficient, by making the preconditions vector only containing the direct preconditions and call the preconditions via a function which also includes indirect ones.
        if (op->precondition_of.empty()) {
            if (op->parent_node) {
                for (unsigned long index = 0; index < op->parent_node->precondition_of.size(); index++) {
                    if (op->parent_node->precondition_of[index] == op) {
                        op->parent_node->precondition_of.erase(op->parent_node->precondition_of.begin() + index);
                        break;
                    }
                    //std::remove(op->parent_node->precondition_of.begin(),
                    //            op->parent_node->precondition_of.end(), static_cast<GraphNode*>(op));
                }
            }
            for (PropositionNode *precond : op->preconditions) {
                for (unsigned long index = 0; index < precond->precondition_of.size(); index++) {
                    if (precond->precondition_of[index] == op) {
                        precond->precondition_of.erase(precond->precondition_of.begin() + index);
                        break;
                    }
                }

                //std::remove(precond->precondition_of.begin(),
                //            precond->precondition_of.end(), static_cast<GraphNode*>(op));
            }
            for (unsigned long index = 0; index < operator_nodes.size(); index++) {
                if (operator_nodes[index] == op) {
                    operator_nodes.erase(operator_nodes.begin() + index);
                    break;
                }
            }
            //std::remove(operator_nodes.begin(), operator_nodes.end(), op);
            delete op;
            counter_deleted_nodes++;
        }
    }


    if (log.is_at_least_normal()) {
        log << " done! [" << operator_nodes.size() << " unary operators, Removed Nodes: " << counter_deleted_nodes << ", Removed Effects: " << counter_deleted_effects
            << "]" << endl;
    }

#ifndef NDEBUG
    //Check that there is no operator with no effect
    for (OperatorNode *op: operator_nodes)
        assert(!op->precondition_of.empty());
    //Check that no deleted node is still in reference
    std::unordered_set<GraphNode*> pointer_set;
    for (OperatorNode *op: operator_nodes)
        pointer_set.insert(op);
    for (PropositionNode *prop : propositions)
        for (GraphNode *child : prop->precondition_of)
            assert(pointer_set.find(child) != pointer_set.end());

    for (PropositionNode *prop : propositions)
        pointer_set.insert(prop);
    for (OperatorNode *op : operator_nodes)
        for (GraphNode *child : op->precondition_of)
            assert(pointer_set.find(child) != pointer_set.end());
#endif
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
