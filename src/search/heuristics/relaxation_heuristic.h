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

using PropQueue = priority_queues::BucketQueue<PropositionNode*>;

struct OperatorNode;
struct PropositionNode;
struct Operator;
struct Proposition;
struct GraphNode {
    int cost; // Used for h^max cost or h^add cost;
    std::vector<OperatorNode*> precondition_of_op;
    std::vector<PropositionNode*> precondition_of_prop;
    explicit GraphNode();
    explicit GraphNode(std::vector<OperatorNode*> &&precondition_of_op, std::vector<PropositionNode*> &&precondition_of_prob);
    virtual ~GraphNode() = default;
//    virtual void update_precondition(PropQueue &queue, GraphNode *predecessor)=0;
};


struct OperatorNode : public GraphNode {
    const int base_cost;
    const int num_preconditions;
    int operator_no; // -1 for axioms; index into the task's operators otherwise
    int unsatisfied_preconditions;
    std::vector<PropositionNode*> preconditions;
    OperatorNode* parent_node;

    explicit OperatorNode(int base_cost, int num_preconditions, int operator_no);
    virtual ~OperatorNode() = default;
    //TODO: delete the copy constructor again

//    void update_precondition(PropQueue &queue, GraphNode *predecessor) override;

    Operator* corresponding_op;
};

struct PropositionNode: public GraphNode {
    PropID prop_id;
    // TODO: Make sure in constructor that reached_by does not overflow.
    OperatorNode* reached_by;
    /* The following two variables are conceptually bools, but Visual C++ does
       not support packing ints and bools together in a bitfield. */
    int num_precondition_occurrences : 30;
    unsigned int is_goal : 1;
    unsigned int marked : 1; // used for preferred operators of h^add and h^FF
    Proposition* corresponding_prop;
    explicit PropositionNode(PropID prop_id);
    virtual ~PropositionNode() = default;
    //TODO: delete the copy constructor again
//    PropositionNode(const PropositionNode &) = delete;
//    void update_precondition(PropQueue &queue, GraphNode *predecessor) override;
//    void update_precondition(PropQueue &queue);
};

//static_assert(sizeof(GraphNode) == 28, "GraphNode has wrong size");

struct Proposition {
    int cost; // Used for h^max cost or h^add cost;
    PropID prop_id;
    // TODO: Make sure in constructor that reached_by does not overflow.
    Operator* reached_by;
    array_pool::ArrayPoolIndex<Operator*> precondition_of_op_index;
    int precondition_of_op_size : 30;
    unsigned int is_goal : 1;
    unsigned int marked : 1;
    explicit Proposition(PropositionNode *);
};

static_assert(sizeof(Proposition) == 24, "Proposition has wrong size");

struct Operator {
    const int base_cost;
    const int num_preconditions;
    int operator_no; // -1 for axioms; index into the task's operators otherwise
    int unsatisfied_preconditions;
    int cost; // Used for h^max cost or h^add cost;
    array_pool::ArrayPoolIndex<Operator*> precondition_of_op_index;
    int precondition_of_op_size;
    array_pool::ArrayPoolIndex<Proposition*> precondition_of_prop_index;
    int precondition_of_prop_size;
    array_pool::ArrayPoolIndex<Proposition*> precondition_index;
    int precondition_size;

    explicit Operator(OperatorNode *);
};

static_assert(sizeof(Operator) == 44, "Operator has wrong size");

class RelaxationHeuristic : public Heuristic {
    void build_unary_operators(const OperatorProxy &op);
    void simplify();

    std::vector<PropID> proposition_offsets;
protected:
    std::vector<PropositionNode*> propositions_nodes;
    std::vector<OperatorNode*> operator_nodes;

    std::vector<PropID> goal_propositions;

    std::vector<Proposition*> propositions;
    std::vector<Operator*> operators;

    array_pool::ArrayPool<Proposition*> prop_precond_of_pool;
    array_pool::ArrayPool<Operator*> op_precond_of_pool;
    array_pool::ArrayPool<Proposition*> preconds_pool;


    /*
      TODO: Some of these protected methods are only needed for the
      CEGAR hack in the additive heuristic and should eventually go
      away.
    */


    PropID get_prop_id(int var, int value) const;
    PropID get_prop_id(const FactProxy &fact) const;

    Proposition *get_proposition(PropID prop_id) {
        return propositions[prop_id];
    }

    Operator *get_operator(OpID op_id) {
        return operators[op_id];
    }


    int get_proposition_cost(int var, int value) const;

public:
    explicit RelaxationHeuristic(const options::Options &options);
    virtual bool dead_ends_are_reliable() const override;

    virtual ~RelaxationHeuristic();
};
}

#endif
