#ifndef HEURISTICS_RELAXATION_HEURISTIC_H
#define HEURISTICS_RELAXATION_HEURISTIC_H

#include "array_pool.h"

#include "../heuristic.h"

#include "../algorithms/priority_queues.h"
#include "../utils/collections.h"

#include <cassert>
#include <vector>

class FactProxy;
class OperatorProxy;

namespace relaxation_heuristic {

using PropID = int;
using OpID = int;

const OpID NO_OP = -1;
const PropID NO_PROP= -1;

class RelaxationHeuristic;
struct PropositionNode;
struct NodePool;

const int OPERATOR = 0;
const int PROPOSITION = 1;

struct SaveLocation {
    explicit SaveLocation(int vector, int index);
    const int vector;
    const int index;
};

using PropQueue = priority_queues::AdaptiveQueue<PropositionNode*>;

struct GraphNode {
    int cost; // Used for h^max cost or h^add cost;
    std::vector<SaveLocation> precondition_of;
    NodePool *np;
    explicit GraphNode();
    explicit GraphNode(NodePool* np);
    explicit GraphNode(std::vector<SaveLocation> &&precondition_of, NodePool* np);
    virtual ~GraphNode() = default;
    virtual void update_precondition(PropQueue &queue, GraphNode *predecessor)=0;
    virtual std::string myname() {return "GraphNode";}
};


struct OperatorNode : public GraphNode {
    const int base_cost;
    const int num_preconditions;
    int operator_no; // -1 for axioms; index into the task's operators otherwise
    int unsatisfied_preconditions;
// preconditions_props(preconditions),
//    PropID effect;

    explicit OperatorNode(int base_cost, int num_preconditions, int operator_no, NodePool *np);
    virtual ~OperatorNode() = default;
    //TODO: delete the copy constructor again
//    OperatorNode(const OperatorNode &) = delete;
//                          int PropID effect,
 virtual void update_precondition(PropQueue &queue, GraphNode *predecessor) override;
    virtual std::string myname() override {return "OperatorNode";}
};

struct PropositionNode: public GraphNode {
    PropID prop_id;
    // TODO: Make sure in constructor that reached_by does not overflow.
    OpID reached_by : 30;
    /* The following two variables are conceptually bools, but Visual C++ does
       not support packing ints and bools together in a bitfield. */
    unsigned int is_goal : 1;
    unsigned int marked : 1; // used for preferred operators of h^add and h^FF
    int num_precondition_occurrences;
    explicit PropositionNode(PropID prop_id, NodePool *np);
    virtual ~PropositionNode() = default;
    //TODO: delete the copy constructor again
//    PropositionNode(const PropositionNode &) = delete;
    virtual void update_precondition(PropQueue &queue, GraphNode *predecessor) override;
    void update_precondition(PropQueue &queue);
    virtual std::string myname() override {return "PropositionNode";}

};

//static_assert(sizeof(GraphNode) == 28, "GraphNode has wrong size");

struct NodePool {
    std::vector<OperatorNode> operator_nodes;
    std::vector<PropositionNode> propositions;

    GraphNode *get_address(SaveLocation save_location);
};

class RelaxationHeuristic : public Heuristic {
    void build_unary_operators(const OperatorProxy &op);
    void simplify();

    // proposition_offsets[var_no]: first PropID related to variable var_no
    std::vector<PropID> proposition_offsets;
protected:
    NodePool node_pool;
    std::vector<PropID> goal_propositions;

//    array_pool::ArrayPool preconditions_pool;
//    array_pool::ArrayPool precondition_of_pool;

//    array_pool::ArrayPoolSlice get_preconditions(OpID op_id) const {
//        const InnerNode &op = (InnerNode&) inner_nodes[op_id];
//        return preconditions_pool.get_slice(op.preconditions_props, op.num_preconditions);
//    }

    // HACK!
//    std::vector<PropID> get_preconditions_vector(OpID op_id) const {
//        auto view = get_preconditions(op_id);
//        return std::vector<PropID>(view.begin(), view.end());
//        std::vector<PropID> temp;
//        return temp;
//    }

    /*
      TODO: Some of these protected methods are only needed for the
      CEGAR hack in the additive heuristic and should eventually go
      away.
    */


//    PropID get_prop_id(const Proposition &prop) const {
//        PropID prop_id = &prop - propositions.data();
//        assert(utils::in_bounds(prop_id, propositions));
//        return prop_id;
//    }
//
//    OpID get_op_id(const EffectNode &op) const {
//        OpID op_id = &op - effect_nodes.data();
//        assert(utils::in_bounds(op_id, effect_nodes));
//        return op_id;
//    }
    int get_num_cond_effects();
//    PropositionNode * get_prop_node(int var, int value) const;
//    PropositionNode * get_prop_node(const FactProxy &fact) const;
    PropID get_prop_id(int var, int value) const;
    PropID get_prop_id(const FactProxy &fact) const;

//    PropositionNode *get_proposition(PropID prop_id) {
//        return propositions[prop_id];
//    }
//    GraphNode *get_operator(OpID op_id) {
//        return &inner_nodes[op_id];
//    }

//    const Proposition *get_proposition(int var, int value) const;
    int get_proposition_cost(int var, int value) const;
//    Proposition *get_proposition(const FactProxy &fact);
public:
    explicit RelaxationHeuristic(const options::Options &options);
    virtual bool dead_ends_are_reliable() const override;
};
}

#endif
