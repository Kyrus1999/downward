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
struct Node;
struct Proposition;
struct Operator;


using PropQueue = priority_queues::AdaptiveQueue<Proposition*>;

struct GraphNode {
    int cost; // Used for h^max cost or h^add cost;
    std::vector<GraphNode*> precondition_of;
    Node* corresponding_node;
    explicit GraphNode();
    explicit GraphNode(std::vector<GraphNode*> &&precondition_of);
    virtual ~GraphNode() = default;
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

};

struct PropositionNode: public GraphNode {
    PropID prop_id;
    // TODO: Make sure in constructor that reached_by does not overflow.
    OperatorNode* reached_by;
    int num_precondition_occurrences : 30;
    /* The following two variables are conceptually bools, but Visual C++ does
       not support packing ints and bools together in a bitfield. */
    unsigned int is_goal : 1;
    unsigned int marked : 1; // used for preferred operators of h^add and h^FF

    explicit PropositionNode(PropID prop_id);
    virtual ~PropositionNode() = default;
    //TODO: delete the copy constructor again
};

struct Node {
    int cost; // Used for h^max cost or h^add cost;
    array_pool::ArrayPoolIndex<Node*> precondition_of_index;
    int precondition_of_size : 30;
    virtual int get_type() = 0;
    explicit Node(int cost);
    virtual ~Node() = default;
};

    struct Proposition : Node {
        PropID prop_id;
        // TODO: Make sure in constructor that reached_by does not overflow.
        Operator* reached_by;
        unsigned int is_goal : 1;
        unsigned int marked : 1;
        int get_type() override {return 0;}
        explicit Proposition(PropositionNode *);
        virtual ~Proposition() = default;
    };

    //static_assert(sizeof(Proposition) == 24, "Proposition has wrong size");
    struct Operator : Node {
        const int base_cost;
        const int num_preconditions;
        int operator_no; // -1 for axioms; index into the task's operators otherwise
        int unsatisfied_preconditions;
        array_pool::ArrayPoolIndex<Proposition*> precondition_index;
        int precondition_size;
        int get_type() override {return 1;}
        explicit Operator(OperatorNode *);
        virtual ~Operator() = default;
    };
    //static_assert(sizeof(Operator) == 44, "Operator has wrong size");

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

    array_pool::ArrayPool<Node*> precond_of_pool;
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

    int get_proposition_cost(int var, int value) const;

public:
    explicit RelaxationHeuristic(const options::Options &options);
    virtual bool dead_ends_are_reliable() const override;

    virtual ~RelaxationHeuristic();

    void update_precondition(PropQueue &queue, Node *predecessor, Node *current);

    void update_precondition(PropQueue &queue, Proposition *prop);

    void update_precondition(PropQueue &queue, Operator *op);
};
}

#endif
